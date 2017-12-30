module ResolutionIdentityModule
export computeTensorElectronRepulsionIntegralsRICoulomb, computeTensorElectronRepulsionIntegralsRIOverlap, computeMatrixExchangeRIK, computeMatrixResolutionOfTheIdentityCoulombMetric, computeMatricesResolutionOfTheIdentityOverlapMetric
using TensorOperations
using ..BaseModule
using ..BasisFunctionsModule
using ..BasisModule
using ..GeometryModule
using ..IntegralsModule

"""
Compute all matrices needed for the Resolution-of-the-Identity (RI) in
overlap metric formulation:

*  B = (μ,ν|P), Sinv = (PQ)^(-1), C = (Q|R)
  (returned in this order)

such that

*  B[μ,ν,:] * Sinv * C * Sinv * B[λ,σ,:] ≈ (μ,ν|λ,σ)
"""
function computeMatricesResolutionOfTheIdentityOverlapMetric(
  basis::GaussianBasis,
  basis_aux::GaussianBasis)
  B = [computeIntegralThreeCenterOverlap(μ,ν,P) for μ in basis.contractedBFs, ν in basis.contractedBFs, P in basis_aux.contractedBFs] #B=(μν|P)
  C = Symmetric([computeIntegralOverlap(cgb1,cgb2) for cgb1 in basis_aux.contractedBFs, cgb2 in basis_aux.contractedBFs]) #Ci=(PQ)^-1
  D = [computeElectronRepulsionIntegral(P,identityCGF,Q,identityCGF)   for P in basis_aux.contractedBFs, Q in basis_aux.contractedBFs] #D=(P|Q)
  return B, inv(C), D
end

"""
Compute the matrix needed for the Resolution-of-the-Identity (RI) in
Coulomb metric formulation:

*  B = (μ,ν|P) * (P|Q)^(-½)

such that

*  B[μ,ν,:] * B[λ,σ,:] ≈ (μ,ν|λ,σ)
"""
function computeMatrixResolutionOfTheIdentityCoulombMetric(
  basis::GaussianBasis,
  basis_aux::GaussianBasis)
  Bpart = [computeElectronRepulsionIntegral(μ,ν,P,identityCGF) for μ in basis.contractedBFs,          ν in basis.contractedBFs,          P in basis_aux.contractedBFs]
  C = Symmetric([computeElectronRepulsionIntegral(P,identityCGF,Q,identityCGF)   for P in basis_aux.contractedBFs, Q in basis_aux.contractedBFs])
  N = length(basis.contractedBFs)
  Naux = length(basis_aux.contractedBFs)
  Bpart = reshape(Bpart,(N*N,Naux))
  return reshape( Bpart*sqrtm(inv(C)) , (N,N,Naux) )
end

computeMatrixResolutionOfTheIdentity = computeMatrixResolutionOfTheIdentityCoulombMetric

function computeTensorElectronRepulsionIntegralsRICoulomb(basis::GaussianBasis, basis_aux::GaussianBasis)
  B = computeMatrixResolutionOfTheIdentityCoulombMetric(basis,basis_aux)
  @tensor ERIs[μ,ν,λ,σ] := B[μ,ν,P] * B[λ,σ,P]
  return ERIs
end

function computeTensorElectronRepulsionIntegralsRIOverlap(basis::GaussianBasis, basis_aux::GaussianBasis)
  B,Ci,D = computeMatricesResolutionOfTheIdentityOverlapMetric(basis,basis_aux) #B=(μν P), Ci=(PQ)^-1 and D=(Q|R)
  N = length(basis.contractedBFs)
  Naux = length(basis_aux.contractedBFs)
  B = reshape(B,(N*N,Naux))
  ERIs = B * Ci * D * Ci * B.'
  return reshape(ERIs,(N,N,N,N))
end

function computeMatrixExchangeRIK(basis::GaussianBasis, basis_aux::GaussianBasis, density::Matrix)
  println("Using RIK")
  N = length(basis.contractedBFs)
  prefactor = 1
  K = zeros(N,N)
  B = computeMatrixResolutionOfTheIdentityCoulombMetric(basis,basis_aux)
  @tensor ERI[μ,ν,λ,σ] := B[μ,ν,P] * B[λ,σ,P]
  @tensor K[μ,λ] := prefactor * ERI[μ,ν,λ,σ] * density[ν,σ]
  return K
end
end #module
