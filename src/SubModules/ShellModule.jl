module ShellModule
export AbstractShell, Shell, expandShell
using ..BaseModule
using ..BasisFunctionsModule
using ..BasisSetModule
using ..GeometryModule
import Base.display

abstract type AbstractShell end

"""
    Shell(center::Position, lqn::LQuantumNumber, exponents::[Float], coefficients::[Float64]; renorm=true)
A shell of basis functions is the collection of basis functions that share a center and radial definition (e.g.
the px, py and pz orbital that are defined together by their center, exponents and coefficients. Contrary to other codes we also consider a shell to have a unique LQuantumNumber not allowing for shared-sp tricks. Working with shells
instead of the individual basis functions allows for some performance gains in integral evaluation.
If renorm is set to false, then the coefficients are read in "as is". With renorm set to true the coefficients are expected to be given as with a basis set file and normalized accordingly. This is the default in accordance with LibInt2Shell, where this is necessary due to performance reasons. (As such the normalization scaling is taken from the renorm() function in libint2/shell.h.
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
      for p in 1:nprim
        coefficients[p] /= sqrt(norm)
      end
    end
    new(center,lqn,exponents,coefficients)
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
    cgbf = ContractedGaussianBasisFunction(sh.coefficients,[])
    for exp in sh.exponents
      push!(cgbf.primitiveBFs,PrimitiveGaussianBasisFunction(sh.center,exp,mqn))
    end
    push!(result,cgbf)
  end
  return result
end

"""
compute the number of basisfunctions (cartesian) described by the shell
"""
function nbf(shell::Shell)
  div((shell.lqn.exponent+1)^2+(shell.lqn.exponent+1),2)
end
end # module
