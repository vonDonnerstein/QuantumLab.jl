module SpecialMatricesModule
export computeMatrixKinetic, computeTensorElectronRepulsionIntegrals, computeMatrixNuclearAttraction, computeMatrixOverlap, computeMatrixCoulomb, computeMatrixExchange, computeMatrixFock, computeMatricesResolutionOfTheIdentityCoulombMetric, computeMatricesResolutionOfTheIdentityOverlapMetric
using ..IntegralsModule
using ..BasisModule
using ..BaseModule
using ..Geometry
using TensorOperations
using ..LibInt2Module

function computeMatrixKinetic(basis::GaussianBasis)
  return [KineticIntegral(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeTensorElectronRepulsionIntegrals(basis::GaussianBasis)
	return ERIs = [computeElectronRepulsionIntegral(μ,ν,λ,σ) for μ in basis.contractedBFs, ν in basis.contractedBFs, λ in basis.contractedBFs, σ in basis.contractedBFs]
end

function computeMatrixOverlap(basis::GaussianBasis)
  return [computeValueOverlap(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixNuclearAttraction(basis::GaussianBasis,geo::Geometry)
  return [sum([NuclearAttractionIntegral(cgb1,cgb2,atom) for atom in geo.atoms]) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatricesResolutionOfTheIdentityOverlapMetric(
  basis::GaussianBasis,
  basis_aux::GaussianBasis)
  # compute B=(μν|P) and Ci=(PQ)^-1 ,Di=(Q|R) such that B Ci D Ci^T B^T ≈ (μν|λσ) (=ERI)
  #
  # example usage:
  # B,Ci,D = computeMatricesResolutionOfTheIdentitOverlapMetric(bas,bas)
  # using TensorOperations
  # @tensor ERI[μ,ν,λ,σ] := B[μ,ν,P] * Ci[P,Q] * D[Q,R] * Ci[R,S] * B[λ,σ,S]
  B = [computeValueThreeCenterOverlap(μ,ν,P) for μ in basis.contractedBFs, ν in basis.contractedBFs, P in basis_aux.contractedBFs] #B=(μν|P)
  C = [computeValueOverlap(cgb1,cgb2) for cgb1 in basis_aux.contractedBFs, cgb2 in basis_aux.contractedBFs] #Ci=(PQ)^-1
  D = [computeIntegral2Center2Electron(P,Q)   for P in basis_aux.contractedBFs, Q in basis_aux.contractedBFs] #D=(P|Q)
	Ci = inv(C)
	DiffCi = Ci - Ci'
	MaxErrorCi = max(abs(DiffCi)...)
	#FOR DEBUG PURRPOSES ONLY
	#= println("Max Element inverse Overlap:  $(max(abs(Ci)...))") =#
	#= println("Max Error inverse Overlap:    $(MaxErrorCi)") =#
	#END OF DEBUGGING BLOCK
  return B, inv(C), D, MaxErrorCi
end

function computeMatricesResolutionOfTheIdentityCoulombMetric(
  basis::GaussianBasis,
  basis_aux::GaussianBasis)
  # compute B=(μν|P) and Ci=(P|Q)^-1, such that B Ci B^T ≈ (μν|λσ) (=ERI)
  #
  # example usage:
  # B,Ci = computeMatricesResolutionOfTheIdentityCoulombMetric(bas,bas)
  # using TensorOperations
  # @tensor ERI[μ,ν,λ,σ] := B[μ,ν,P] * Ci[P,Q] * B[λ,σ,Q]
  B = [computeIntegral3Center2Electron(μ,ν,P) for μ in basis.contractedBFs,          ν in basis.contractedBFs,          P in basis_aux.contractedBFs]
  C = [computeIntegral2Center2Electron(P,Q)   for P in basis_aux.contractedBFs, Q in basis_aux.contractedBFs]
  return B, inv(C)
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
  density::Matrix,
	computeMatrixCoulomb=computeMatrixCoulomb,
	computeMatrixExchange=computeMatrixExchange)
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
  electronRepulsionIntegralsColoumb::Array{Float64,4},
  electronRepulsionIntegralsExchange::Array{Float64,4}=electronRepulsionIntegralsColoumb)
  T = matrixKinetic
  V = matrixNuclearAttraction
  J = computeMatrixCoulomb(electronRepulsionIntegralsColoumb,density)
  K = computeMatrixExchange(electronRepulsionIntegralsExchange,density)
  return T+V+2J-K
end

function computeDimensions(shells::Vector{LibInt2Shell})
  totaldim = 0
  maxprims = 0
  maxlqn = 0
  for sh in shells
    totaldim += div((lqn(sh)+1)^2+(lqn(sh)+1),2)
    maxprims = maxprims < nprims(sh) ? nprims(sh) : maxprims
    maxlqn = maxlqn < lqn(sh) ? lqn(sh) : maxlqn
  end
  return (totaldim, maxprims, maxlqn)
end

function computeMatrixOverlap(shells::Vector{LibInt2Shell})
  totaldim, maxprims, maxlqn = computeDimensions(shells)
  result = zeros(Float64,totaldim,totaldim)
  engine = LibInt2EngineOverlap(maxprims,LQuantumNumber(maxlqn))

  # column-major order (sh1 counts top to bottom, sh2 counts left to right)
  lastcol,lastrow = 0,0
  for sh1 in shells
    block_col_size = div((lqn(sh1)+1)^2+(lqn(sh1)+1),2)
    lastcol += block_col_size
    for sh2 in shells
      block_row_size = div((lqn(sh2)+1)^2+(lqn(sh2)+1),2)
      lastrow += block_row_size

      result[lastcol-block_col_size+1:lastcol, lastrow-block_row_size+1:lastrow] = computeMatrixBlockOverlap(engine,sh1,sh2)
    end
    lastrow = 0
  end

  destroy!(engine)
  return result
end

function computeMatrixKinetic(shells::Vector{LibInt2Shell})
  totaldim, maxprims, maxlqn = computeDimensions(shells)
  result = zeros(Float64,totaldim,totaldim)
  engine = LibInt2EngineKinetic(maxprims,LQuantumNumber(maxlqn))

  # column-major order (sh1 counts top to bottom, sh2 counts left to right)
  lastcol,lastrow = 0,0
  for sh1 in shells
    block_col_size = div((lqn(sh1)+1)^2+(lqn(sh1)+1),2)
    lastcol += block_col_size
    for sh2 in shells
      block_row_size = div((lqn(sh2)+1)^2+(lqn(sh2)+1),2)
      lastrow += block_row_size

      result[lastcol-block_col_size+1:lastcol, lastrow-block_row_size+1:lastrow] = computeMatrixBlockKinetic(engine,sh1,sh2)
    end
    lastrow = 0
  end

  destroy!(engine)
  return result
end

function computeMatrixCoulomb(shells::Vector{LibInt2Shell}, density::Matrix)
  totaldim, maxprims, maxlqn = computeDimensions(shells)
  result = zeros(Float64,totaldim,totaldim)
  engine = LibInt2EngineCoulomb(maxprims,LQuantumNumber(maxlqn))

  # column-major order (sh1 counts top to bottom ("col"-direction), sh2 counts left to right ("row"-direction)
  lastcol,lastrow = 0,0
  lastcol_in_P,lastrow_in_P = 0,0
  for sh1 in shells
    block_col_size = div((lqn(sh1)+1)^2+(lqn(sh1)+1),2)
    lastcol += block_col_size
    for sh2 in shells
      block_row_size = div((lqn(sh2)+1)^2+(lqn(sh2)+1),2)
      lastrow += block_row_size
      for sh3 in shells
	P_block_col_size = div((lqn(sh3)+1)^2+(lqn(sh3)+1),2)
	lastcol_in_P += P_block_col_size
	for sh4 in shells
	  P_block_row_size = div((lqn(sh4)+1)^2+(lqn(sh4)+1),2)
	  lastrow_in_P += P_block_row_size
	  P = density[lastcol_in_P-P_block_col_size+1:lastcol_in_P, lastrow_in_P-P_block_row_size+1:lastrow_in_P]

	  eris = computeElectronRepulsionIntegral(engine,sh1,sh2,sh3,sh4)
	  @tensor block[μ,ν] := eris[μ,ν,λ,σ] * P[λ,σ]
	  result[lastcol-block_col_size+1:lastcol, lastrow-block_row_size+1:lastrow] += block
	end
	lastrow_in_P = 0
      end
      lastcol_in_P = 0
    end
    lastrow = 0
  end

  destroy!(engine)
  return result
end

end
