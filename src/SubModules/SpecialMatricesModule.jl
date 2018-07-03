module SpecialMatricesModule
export computeMatrixKinetic, computeTensorElectronRepulsionIntegrals, computeMatrixNuclearAttraction, computeMatrixOverlap, computeMatrixCoulomb, computeMatrixExchange, computeMatrixFock
using ..BaseModule
using ..AtomModule
using ..IntegralsModule
using ..BasisModule
using ..Geometry
using TensorOperations
using ..LibInt2Module
using ..ShellModule
using ..ShellModule.nbf
import ..IntegralsModule.computeTensorBlockElectronRepulsionIntegrals

#HELPERS
function scatterMatrixBlocks2D(shells,totaldim::Integer,blockevaluator::Function,blocklength::Function)
  result = zeros(Float64,totaldim,totaldim)
  # column-major order (sh1 counts top to bottom, sh2 counts left to right)
  lastcol,lastrow = 0,0
  for sh1 in shells
    block_col_size = blocklength(sh1)
    lastcol += block_col_size
    for sh2 in shells
      block_row_size = blocklength(sh2)
      lastrow += block_row_size
      result[lastcol-block_col_size+1:lastcol, lastrow-block_row_size+1:lastrow] = blockevaluator(sh1,sh2)
    end
    lastrow = 0
  end
  return result
end

function contractTensorBlocksWithMatrix_CoulombLike(shells,density,totaldim::Integer,tensorblockevaluator::Function,blocklength::Function)
  result = zeros(Float64,totaldim,totaldim)
  # column-major order (sh1 counts top to bottom ("col"-direction), sh2 counts left to right ("row"-direction)
  lastcol,lastrow = 0,0
  lastcol_in_P,lastrow_in_P = 0,0
  for sh1 in shells
    block_col_size = blocklength(sh1)
    lastcol += block_col_size
    for sh2 in shells
      block_row_size = blocklength(sh2)
      lastrow += block_row_size
      for sh3 in shells
	P_block_col_size = blocklength(sh3)
	lastcol_in_P += P_block_col_size
	for sh4 in shells
	  P_block_row_size = blocklength(sh4)
	  lastrow_in_P += P_block_row_size
	  P = density[lastcol_in_P-P_block_col_size+1:lastcol_in_P, lastrow_in_P-P_block_row_size+1:lastrow_in_P]
	  tensor = tensorblockevaluator(sh1,sh2,sh3,sh4)
	  #@tensor block[μ,ν] := tensor[μ,ν,λ,σ] * P[λ,σ]
	  block = zeros(block_col_size,block_row_size)
	  for μ in 1:block_col_size
	    for ν in 1:block_row_size
	      block[μ,ν] = dot(vec(tensor[μ,ν,:,:]),vec(P))
	    end
	  end
	  #end tensor stuff
	  result[lastcol-block_col_size+1:lastcol, lastrow-block_row_size+1:lastrow] += block
	end
	lastrow_in_P = 0
      end
      lastcol_in_P = 0
    end
    lastrow = 0
  end
  return result
end

function contractTensorBlocksWithMatrix_ExchangeLike(shells,density,totaldim::Integer,tensorblockevaluator::Function,blocklength::Function)
  result = zeros(Float64,totaldim,totaldim)
  # column-major order (sh1 counts top to bottom ("col"-direction), sh2 counts left to right ("row"-direction)
  lastcol,lastrow = 0,0
  lastcol_in_P,lastrow_in_P = 0,0
  for sh1 in shells
    block_col_size = blocklength(sh1)
    lastcol += block_col_size
    for sh2 in shells
      block_row_size = blocklength(sh2)
      lastrow += block_row_size
      for sh3 in shells
	P_block_col_size = blocklength(sh3)
	lastcol_in_P += P_block_col_size
	for sh4 in shells
	  P_block_row_size = blocklength(sh4)
	  lastrow_in_P += P_block_row_size
	  P = density[lastcol_in_P-P_block_col_size+1:lastcol_in_P, lastrow_in_P-P_block_row_size+1:lastrow_in_P]
	  tensor = tensorblockevaluator(sh1,sh3,sh2,sh4)
	  @tensor block[μ,ν] := tensor[μ,λ,ν,σ] * P[λ,σ]
	  result[lastcol-block_col_size+1:lastcol, lastrow-block_row_size+1:lastrow] += block
	end
	lastrow_in_P = 0
      end
      lastcol_in_P = 0
    end
    lastrow = 0
  end
  return result
end

#OVERLAP
function computeMatrixOverlap(basis::GaussianBasis)
  return [computeIntegralOverlap(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixOverlap(shells::Vector{Shell})
  totaldim = mapreduce(nbf,+,0,shells)
  return scatterMatrixBlocks2D(shells,totaldim,computeMatrixBlockOverlap,nbf)
end

if (LibInt2Shell != Shell)
  function computeMatrixOverlap(shells::Vector{LibInt2Shell})
    totaldim, maxprims, maxlqn = LibInt2Module.computeDimensions(shells)
    engine = LibInt2EngineOverlap(maxprims,maxlqn)
    result = scatterMatrixBlocks2D(shells,totaldim,(sh1,sh2)->computeMatrixBlockOverlap(engine,sh1,sh2),nbf)
    destroy!(engine)
    return result
  end
end

#KINETIC
function computeMatrixKinetic(basis::GaussianBasis)
  return [computeIntegralKinetic(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixKinetic(shells::Vector{Shell})
  totaldim = mapreduce(nbf,+,0,shells)
  return scatterMatrixBlocks2D(shells,totaldim,computeMatrixBlockKinetic,nbf)
end

if (LibInt2Shell != Shell)
  function computeMatrixKinetic(shells::Vector{LibInt2Shell})
    totaldim, maxprims, maxlqn = LibInt2Module.computeDimensions(shells)
    engine = LibInt2EngineKinetic(maxprims,maxlqn)
    result = scatterMatrixBlocks2D(shells,totaldim,(sh1,sh2)->computeMatrixBlockKinetic(engine,sh1,sh2),nbf)
    destroy!(engine)
    return result
  end
end

#NUCLEARATTRACTION
function computeMatrixNuclearAttraction(basis::GaussianBasis,geo::Geometry)
  return [sum([computeIntegralNuclearAttraction(cgb1,cgb2,atom) for atom in geo.atoms]) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixNuclearAttraction(shells::Vector{Shell},atom::Atom)
  totaldim = mapreduce(nbf,+,0,shells)
  return scatterMatrixBlocks2D(shells,totaldim,(sh1,sh2)->computeMatrixBlockNuclearAttraction(sh1,sh2,atom),nbf)
end

function computeMatrixNuclearAttraction(shells::Vector{Shell},geo::Geometry)
  totaldim = mapreduce(nbf,+,0,shells)
  result = zeros(totaldim,totaldim)
  for atom in geo.atoms
    result += scatterMatrixBlocks2D(shells,totaldim,(sh1,sh2)->computeMatrixBlockNuclearAttraction(sh1,sh2,atom),nbf)
  end
  return result
end

if (LibInt2Shell != Shell)
  function computeMatrixNuclearAttraction(shells::Vector{LibInt2Shell},atom::Atom)
    totaldim, maxprims, maxlqn = LibInt2Module.computeDimensions(shells)
    engine = LibInt2EngineNuclearAttraction(maxprims,maxlqn,atom)
    result = scatterMatrixBlocks2D(shells,totaldim,(sh1,sh2)->computeMatrixBlockNuclearAttraction(engine,sh1,sh2),nbf)
    destroy!(engine)
    return result
  end

  function computeMatrixNuclearAttraction(shells::Vector{LibInt2Shell},geo::Geometry)
    totaldim, maxprims, maxlqn = LibInt2Module.computeDimensions(shells)
    engine = LibInt2EngineNuclearAttraction(maxprims,maxlqn,geo)
    result = scatterMatrixBlocks2D(shells,totaldim,(sh1,sh2)->computeMatrixBlockNuclearAttraction(engine,sh1,sh2),nbf)
    destroy!(engine)
    return result
  end
end

#ERIS
function computeTensorElectronRepulsionIntegrals(basis::GaussianBasis)
	return ERIs = [computeElectronRepulsionIntegral(μ,ν,λ,σ) for μ in basis.contractedBFs, ν in basis.contractedBFs, λ in basis.contractedBFs, σ in basis.contractedBFs]
end

function computeTensorBlockElectronRepulsionIntegrals(μ::Shell,ν::Shell,λ::Shell,σ::Shell)
  [computeElectronRepulsionIntegral(μcgbf,νcgbf,λcgbf,σcgbf) for μcgbf in expandShell(μ), νcgbf in expandShell(ν), λcgbf in expandShell(λ), σcgbf in expandShell(σ)]
end

function computeTensorElectronRepulsionIntegrals(shells::Vector{LibInt2Shell})
  totaldim, maxprims, maxlqn = LibInt2Module.computeDimensions(shells)
  engine = LibInt2EngineCoulomb(maxprims,maxlqn)
  result = computeTensorElectronRepulsionIntegrals(shells,computeBlockERI=((μ,ν,λ,σ)->copy(computeTensorBlockElectronRepulsionIntegrals(engine,μ,ν,λ,σ))))
  destroy!(engine)
  return result
end

function computeTensorElectronRepulsionIntegrals(shells::Union{Vector{Shell},Vector{LibInt2Shell}}; computeBlockERI=computeTensorBlockElectronRepulsionIntegrals)
  cat(4,
      [cat(3,
	   [cat(2,
		[cat(1,
		     [computeBlockERI(μsh,νsh,λsh,σsh) for μsh in shells]...) 
		 for νsh in shells]...)
	    for λsh in shells]...)
       for σsh in shells]...)
end

#COULOMB
function computeMatrixCoulomb(basis::GaussianBasis, density::Matrix)
  N = length(basis.contractedBFs)

  J = zeros(N,N)
  μIndex = 0
  for μ in basis.contractedBFs
    μIndex+=1
    for νIndex in 1:μIndex
      ν = basis.contractedBFs[νIndex]
      ERImatrix = [computeElectronRepulsionIntegral(μ,ν,λ,σ) for λ in basis.contractedBFs, σ in basis.contractedBFs]
      J[μIndex,νIndex] = trace(ERImatrix*density)
      J[νIndex,μIndex] = J[μIndex,νIndex]
    end
  end
  return J
end

function computeMatrixCoulomb(electronRepulsionIntegrals::Array{Float64,4}, density::Matrix)
  @tensor J[μIndex,νIndex] := electronRepulsionIntegrals[μIndex,νIndex,λIndex,σIndex] * density[λIndex,σIndex]
end

function computeMatrixCoulomb(shells::Vector{Shell}, density)
  totaldim = mapreduce(nbf,+,0,shells)
  contractTensorBlocksWithMatrix_CoulombLike(shells,density,totaldim,computeTensorBlockElectronRepulsionIntegrals,nbf)
end

if (LibInt2Shell != Shell)
  function computeMatrixCoulomb(shells::Vector{LibInt2Shell}, density::Matrix)
    totaldim, maxprims, maxlqn = LibInt2Module.computeDimensions(shells)
    engine = LibInt2EngineCoulomb(maxprims,maxlqn)
    result = contractTensorBlocksWithMatrix_CoulombLike(shells,density,totaldim,(μ,ν,λ,σ)->computeTensorBlockElectronRepulsionIntegrals(engine,μ,ν,λ,σ),nbf)
    destroy!(engine)
    return result
  end
end

#EXCHANGE
function computeMatrixExchange(basis::GaussianBasis, density::Matrix)
  N = length(basis.contractedBFs)

  K = zeros(N,N)
  μIndex = 0
  for μ in basis.contractedBFs
    μIndex+=1
    for νIndex in 1:μIndex
      ν = basis.contractedBFs[νIndex]
      ERImatrix = [computeElectronRepulsionIntegral(μ,λ,ν,σ) for λ in basis.contractedBFs, σ in basis.contractedBFs]
      K[μIndex,νIndex] = trace(ERImatrix*density)
      K[νIndex,μIndex] = K[μIndex,νIndex]
    end
  end
  return K
end

function computeMatrixExchange(electronRepulsionIntegrals::Array{Float64,4}, density::Matrix)
  @tensor K[μIndex,νIndex] := electronRepulsionIntegrals[μIndex,λIndex,νIndex,σIndex] * density[λIndex,σIndex]
end

function computeMatrixExchange(shells::Vector{Shell}, density)
  totaldim = mapreduce(nbf,+,0,shells)
  contractTensorBlocksWithMatrix_ExchangeLike(shells,density,totaldim,computeTensorBlockElectronRepulsionIntegrals,nbf)
end

if (LibInt2Shell != Shell)
  function computeMatrixExchange(shells::Vector{LibInt2Shell}, density::Matrix)
    totaldim, maxprims, maxlqn = LibInt2Module.computeDimensions(shells)
    engine = LibInt2EngineCoulomb(maxprims,maxlqn)
    result = contractTensorBlocksWithMatrix_ExchangeLike(shells,density,totaldim,(μ,ν,λ,σ)->computeTensorBlockElectronRepulsionIntegrals(engine,μ,ν,λ,σ),nbf)
    destroy!(engine)
    return result
  end
end

#FOCK
function computeMatrixFock(
  matrixKinetic::Matrix,
  matrixNuclearAttraction::Matrix,
  matrixCoulomb::Matrix,
  matrixExchange::Matrix)
  return matrixKinetic + matrixNuclearAttraction + 2*matrixCoulomb-matrixExchange
end

function computeMatrixFock(
  basis::GaussianBasis,
  geometry::Geometry,
  density::Matrix,
	computeMatrixCoulomb=computeMatrixCoulomb,
	computeMatrixExchange=computeMatrixExchange)
  T = computeMatrixKinetic(basis)
  V = computeMatrixNuclearAttraction(basis,geometry)
  J = computeMatrixCoulomb(basis,density)
  K = computeMatrixExchange(basis,density)
  return computeMatrixFock(T,V,J,K)
end

function computeMatrixFock(
  density::Matrix,
  matrixKinetic::Matrix,
  matrixNuclearAttraction::Matrix,
  electronRepulsionIntegralsColoumb::Array{Float64,4},
  electronRepulsionIntegralsExchange::Array{Float64,4}=electronRepulsionIntegralsColoumb)
  T = matrixKinetic
  V = matrixNuclearAttraction
  J = computeMatrixCoulomb(electronRepulsionIntegralsColoumb,density)
  K = computeMatrixExchange(electronRepulsionIntegralsExchange,density)
  return computeMatrixFock(T,V,J,K)
end

end
