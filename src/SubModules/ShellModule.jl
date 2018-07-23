module ShellModule
export AbstractShell, Shell, expandShell, SparseMatrixBCSR, SparseMatrixBCSRSymmetric
using ..BaseModule
using ..BasisFunctionsModule
using ..BasisSetModule
using ..GeometryModule
using ..MatrixModule
import Base.display
import ..MatrixModule.computeBlockingPattern
import ..MatrixModule.SparseMatrixBCSR
import ..MatrixModule.SparseMatrixBCSRSymmetric

abstract type AbstractShell end

"""
    Shell(center::Position, lqn::LQuantumNumber, exponents::[Float], coefficients::[Float64]; renorm=true)

A shell of basis functions is the collection of basis functions that share a center and radial definition (e.g.
the px, py and pz orbital that are defined together by their center, exponents and coefficients. Contrary to other codes we also consider a shell to have a unique LQuantumNumber not allowing for shared-sp tricks. Working with shells
instead of the individual basis functions allows for some performance gains in integral evaluation.
If renorm is set to false, then the coefficients are read in "as is". With renorm set to true the coefficients are expected to be given as with a basis set file and normalized accordingly. This is the default in accordance with LibInt2Shell, where this is necessary due to performance reasons. (As such the normalization scaling is taken from the renorm() function in libint2/shell.h.

In a shell, there is no discrimination between MQuantumNumbers. This allows for efficient recursion schemes, but poses a challenge to normalization, because the normalization factor depends explicitly on the MQuantumNumber. Following the choice of libint2 (see: https://github.com/evaleev/libint/blob/80d829b0e9932cd118a85e8aab03a41fe0699a7b/doc/progman/progman.tex#L498) to scale all coefficients according to the axial MQuantumNumbers (e.g. x, yyy, zz, ...). After evaluation of the integrals, those including non-axial (xy, xz, ...) MQuantumNumbers need to be rescaled by the relative factors (`renormNonAxial!`).

[Fermann,Valeev: Fundamentals of Molecular Integrals Evaluation, Eq. 2.25]
"""
type Shell <: AbstractShell
  center::Position
  lqn::LQuantumNumber
  exponents::Vector{Float64}
  coefficients::Vector{Float64}
  function Shell(center,lqn,exponents,coefficients;renorm=true)
    L = lqn.exponent
    if renorm
      nprim = length(exponents)
      # renorm primitive-by-primitive
      for p in 1:nprim
	coefficients[p] *= sqrt( (2^L*(2*exponents[p])^(L+3/2)) / (π^(3/2)*doublefactorial(2*L-1)) )
      end
      # normalize total shell
      norm = 0.
      for p in 1:nprim, q in 1:nprim
        norm += (doublefactorial(2*L-1)*π^(3/2)*coefficients[p]*coefficients[q])/
	              ((2^L)*(exponents[p]+exponents[q])^(L+3/2))
      end
      coefficients /= sqrt(norm)
    end
    new(center,lqn,exponents,coefficients)
  end
end

function computeFactorsNonAxial(lqn::LQuantumNumber)
  L = lqn.exponent
  result = Vector{Float64}()
  for mqn in MQuantumNumbers(lqn)
    na = doublefactorial(2*mqn.x-1)*doublefactorial(2*mqn.y-1)*doublefactorial(2*mqn.z-1)/doublefactorial(2*L-1)
    push!(result,1/√na)
  end
  return result
end

function computeFactorsNonAxial(shells::Vector{Shell})
  naxfacs = Vector{Float64}()
  for sh in shells
    push!(naxfacs,computeFactorsNonAxial(sh.lqn)...)
  end
  return naxfacs
end

function renormNonAxial!(vec::Vector{Float64},shells::Vector{<:AbstractShell})
  naxfacs = computeFactorsNonAxial(shells)
  vec .*= naxfacs
end

function renormNonAxial!(cgbfs::Vector{ContractedGaussianBasisFunction},shells::Vector{<:AbstractShell})
  naxfacs = computeFactorsNonAxial(shells)
  for (cgbf,fac) in zip(cgbfs,naxfacs)
    scale!(cgbf.coefficients,fac)
  end
end

function renormNonAxial!(mat::Matrix{Float64},shells::Vector{<:AbstractShell})
  naxfacs = computeFactorsNonAxial(shells)
  scale!(naxfacs,mat)
  scale!(mat,naxfacs)
end

function renormNonAxial!(mat::Matrix{Float64},shells1::Vector{<:AbstractShell},shells2::Vector{<:AbstractShell})
  naxfacs1 = computeFactorsNonAxial(shells1)
  scale!(naxfacs1,mat)
  naxfacs2 = computeFactorsNonAxial(shells2)
  scale!(mat,naxfacs2)
end

function renormNonAxial!(tensor::Array{Float64,N},shells::Vararg{Vector{<:AbstractShell},N}) where {N}
  naxfacs = computeFactorsNonAxial.(shells)
  for idx in CartesianRange(size(tensor))
    for dim in 1:length(shells)
      tensor[idx] *= naxfacs[dim][idx.I[dim]]
    end
  end
end

function display(sh::Shell)
  println("$(sh.lqn.symbol)-Shell   @    $(sh.center)")
  print("  Exponents:    ")
  for exp in sh.exponents
    @printf("  %f",exp)
  end
  println("")
  print("  Coefficients: ")
  for coeff in sh.coefficients
    @printf("  %f",coeff)
  end
  println("")
end

"""
Compute a list of ContractedGaussianBasisFunctions which are described within the input shell.
"""
function expandShell(sh::Shell)
  result = ContractedGaussianBasisFunction[]
  for mqn in MQuantumNumbers(sh.lqn)
    pgbfs = [PrimitiveGaussianBasisFunction(sh.center,exp,mqn) for exp in sh.exponents]
    cgbf = ContractedGaussianBasisFunction(copy(sh.coefficients),pgbfs)
    push!(result,cgbf)
  end
  renormNonAxial!(result,[sh])
  return result
end

"""
    nbf(shell)
computes the number of basisfunctions (cartesian) described by the shell.
"""
nbf(shell::Shell) = BaseModule.numberMQNsCartesian(shell.lqn)

"""
    computeBlockingPattern(shells)
returns the matrix blocking pattern for given shells, e.g. in case of p orbitals a block consists of [p.x,p.y,p.z]
"""
function computeBlockingPattern(shells::Array{Shell,1})
    n::Array{Int64,1}     = []
    pattern::Array{Tuple{Int64,Int64},1} = []
    sum::Int64							 = 0

    [push!(n,nbf(shells[i])) for i in 1:length(shells)]
    δ::Array{Int64,1} = [n[i+1]-n[i] for i in 1:length(n)-1]
    δ = vcat(0,δ)
    [if(δ[i]<0) δ[i] = 0. end for i in 1:length(n)]
    for i = 1:length(shells)
        sum += 1
		push!(pattern,(sum,sum+n[i]-1))
        sum = pattern[end][2]
    end
	
    return Array{Tuple{Int64,Int64},1}(pattern)
end

"""
    computeIncreasedBlockSizePattern(shells)
returns the matrix blocking pattern with a minimum block size of 100x100
"""
#TODO: Atome nicht in der Mitte zerschneiden
function computeIncreasedBlockSizePattern(shells::Array{Shell,1})
	n::Array{Int64,1} = []
	[push!(n,nbf(shells[i])) for i in 1:length(shells)]
	sum1::Int64 = 0
	sum2::Int64 = 0	 
	pattern::Array{Tuple{Int64,Int64},} = []
	
	for i = 1:length(n)
		sum2 += n[i]
		if sum2 > 100
			push!(pattern,(sum1+1,sum1+sum2))
			sum1 += sum2
			sum2 = 0
		elseif i == length(n)
			push!(pattern,(sum1+1,sum1+sum2))
		end
	end
	
	return pattern
end

"""
    SparseMatrixBCSR(M,shells,Bool)
returns a matrix M in SparseMatrixBCSR format with blocking pattern from shells, if Bool = true 100x100 blocks are used	
"""
function SparseMatrixBCSR(M::Array{Float64,2},shells::Array{Shell,1},b::Bool)
	if b == true
		p = computeIncreasedBlockSizePattern(shells::Array{Shell,1})
		pattern::Array{Tuple{Int64,Int64},1} = []
		for i = 1:length(p)
			push!(pattern,(p[i][1],p[i][2]))
		end
		return SparseMatrixBCSR(M::Array{Float64,2},pattern,pattern)
	else
		return SparseMatrixBCSR(M::Array{Float64,2},computeBlockingPattern(shells::Array{Shell,1}),computeBlockingPattern(shells::Array{Shell,1}))
	end
end
	
"""
    SparseMatrixBCSRSymmetric(M,shells)
returns a matrix M in SparseMatrixBCSRSymmetric format with blocking pattern from shells
"""
function SparseMatrixBCSRSymmetric(M::Array{Float64,2},shells::Array{ShellModule.Shell,1})                                                                              
	SparseMatrixBCSRSymmetric(M::Array{Float64,2},computeBlockingPattern(shells::Array{ShellModule.Shell,1}))
end

end #end of module
