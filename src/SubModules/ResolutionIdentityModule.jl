module ResolutionIdentityModule
export computeTensorElectronRepulsionIntegralsRICoulomb, computeTensorElectronRepulsionIntegralsRIOverlap, computeMatrixExchangeRIK
using TensorOperations
using ..BaseModule
using ..BasisModule
using ..GeometryModule
using ..IntegralsModule
using ..SpecialMatricesModule

function computeTensorElectronRepulsionIntegralsRICoulomb(basis::GaussianBasis, basis_aux::GaussianBasis)
  B,Ci = computeMatricesResolutionOfTheIdentityCoulombMetric(basis,basis_aux) #B=(μν|P) and Ci=(P|Q)^-1
  @tensor ERIs[μ,ν,λ,σ] := B[μ,ν,P] * Ci[P,Q] * B[λ,σ,Q]
  return ERIs
end

function computeTensorElectronRepulsionIntegralsRIOverlap(basis::GaussianBasis, basis_aux::GaussianBasis)
  B,Ci,D = computeMatricesResolutionOfTheIdentityOverlapMetric(basis,basis_aux) #B=(μν P), Ci=(PQ)^-1 and D=(Q|R)
  @tensor ERIs[μ,ν,λ,σ] := B[μ,ν,P] * Ci[P,Q] * D[Q,R] * Ci[R,S] * B[λ,σ,S]
  return ERIs
end

function computeMatrixExchangeRIK(basis::GaussianBasis, basis_aux::GaussianBasis, density::Matrix)
  println("Using RIK")
  N = length(basis.contractedBFs)
  prefactor = 1
  K = zeros(N,N)
  B,Ci = computeMatricesResolutionOfTheIdentityCoulombMetric(basis,basis_aux) #B=(μν|P) and Ci=(P|Q)^-1
  @tensor ERI[μ,ν,λ,σ] := B[μ,ν,P] * Ci[P,Q] * B[λ,σ,Q]
  @tensor K[μ,λ] := prefactor * ERI[μ,ν,λ,σ] * density[ν,σ]
  return K
end
end #module
