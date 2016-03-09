module SpecialMatricesModule
export computeMatrixKinetic, computeMatrixNuclearAttraction, computeMatrixOverlap, computeMatrixCoulomb, computeMatrixExchange, computeMatrixFock
using ..IntegralsModule
using ..BasisModule
using ..Geometry
using TensorOperations

function computeMatrixKinetic(basis::GaussianBasis)
  return [KineticIntegral(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixOverlap(basis::GaussianBasis)
  return [Overlap(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixNuclearAttraction(basis::GaussianBasis,geo::Geometry)
  return [sum([NuclearAttractionIntegral(cgb1,cgb2,atom) for atom in geo.atoms]) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixCoulomb(basis::GaussianBasis, density::Matrix)
  N = length(basis.contractedBFs)

  J = zeros(N,N)
  μIndex = 0
  for (μ in basis.contractedBFs)
    μIndex+=1
    for (νIndex in 1:μIndex)
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

function computeMatrixExchange(basis::GaussianBasis, density::Matrix)
  N = length(basis.contractedBFs)

  K = zeros(N,N)
  μIndex = 0
  for (μ in basis.contractedBFs)
    μIndex+=1
    for (νIndex in 1:μIndex)
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


function computeMatrixFock(
  basis::GaussianBasis,
  geometry::Geometry,
  density::Matrix)
  T = computeMatrixKinetic(basis)
  V = computeMatrixNuclearAttraction(basis,geometry)
  J = computeMatrixCoulomb(basis,density)
  K = computeMatrixExchange(basis,density)
  return T+V+2J-K
end

function computeMatrixFock(
  density::Matrix,
  matrixKinetic::Matrix,
  matrixNuclearAttraction::Matrix,
  electronRepulsionIntegrals::Array{Float64,4})
  T = matrixKinetic
  V = matrixNuclearAttraction
  J = computeMatrixCoulomb(electronRepulsionIntegrals,density)
  K = computeMatrixExchange(electronRepulsionIntegrals,density)
  return T+V+2J-K
end

end
