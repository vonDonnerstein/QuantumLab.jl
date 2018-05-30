module HartreeFockModule
export computeEnergyHartreeFock, evaluateSCF
using ..BaseModule
using ..BasisModule
using ..GeometryModule
using ..SpecialMatricesModule
using ..IntegralsModule
using ..ShellModule
using ..LibInt2Module

function computeEnergyHartreeFock(
  density::Matrix,
  kinetic::Matrix,
  nuclearAttraction::Matrix,
  coulomb::Matrix,
  exchange::Matrix,
  interatomicRepulsion::Float64;
  info::Bool=true)

  P = density
  info && println("=================================")
  eT = 2*trace(kinetic*P)
  info && println("Kinetic Energy: $eT")
  eV = 2*trace(nuclearAttraction*P)
  info && println("Nuclear Attraction Energy: $eV")
  info && println("Sum 1-Electron Energy: $(eT+eV)")
  info && println("---------------------------------")
  # G = 2J-K
  eJ = 2*trace(coulomb*P)
  info && println("Coulomb Energy: $eJ")
  eK = -trace(exchange*P)
  info && println("Exchange Energy: $eK")
  info && println("Total Electronic Energy: $(eT+eV+eJ+eK)")
  info && println("Total Energy: $(eT+eV+eJ+eK+interatomicRepulsion)")
  info && println("=================================")

  return eT+eV+eJ+eK
end

function computeEnergyHartreeFock(
  basis::GaussianBasis,
  geometry::Geometry,
  density::Matrix;
  info::Bool=true)

  P = density
  T = computeMatrixKinetic(basis)
  V = computeMatrixNuclearAttraction(basis,geometry)
  J = computeMatrixCoulomb(basis,density)
  K = computeMatrixExchange(basis,density)
  Vnn = computeEnergyInteratomicRepulsion(geometry)
  computeEnergyHartreeFock(P,T,V,J,K,Vnn;info=info)
end

function computeEnergyHartreeFock(
  density::Matrix,
  matrixKinetic::Matrix,
  matrixNuclearAttraction::Matrix,
  electronRepulsionIntegrals::Array{Float64,4},
  energyInteratomicRepulsion::Float64;
  info::Bool=true,
  # inceptions
  electronRepulsionIntegralsExchange::Array{Float64,4}=electronRepulsionIntegrals)

  J = computeMatrixCoulomb(electronRepulsionIntegrals,density)
  K = computeMatrixExchange(electronRepulsionIntegralsExchange,density)
  computeEnergyHartreeFock(density, matrixKinetic, matrixNuclearAttraction, J, K, energyInteratomicRepulsion;info=info)
end

function evaluateSCFStep(
  initialGuessDensity::Matrix,
  fock::Matrix,
  overlap::Matrix,
  electronNumber::Integer)

  eig = eigfact(Symmetric(fock),Symmetric(overlap)) # because Symmetric is faster and allows for sorted eigenvalues
  C = eig.vectors
  moEnergies = eig.values
  Cocc = eig.vectors[:,1:electronNumber]
  P = Cocc * Cocc.'

  return (moEnergies,P)
end

computeMatrixCoefficients(fock::Matrix, overlap::Matrix) = eigvecs(Symmetric(fock),Symmetric(overlap))

function evaluateSCF(
  initialGuessDensity::Matrix,
  overlap::Matrix,
  kinetic::Matrix,
  nuclearAttraction::Matrix,
  coulomb::Function,
  exchange::Function,
  interatomicRepulsion::Float64,
  electronNumber::Integer;
  info::Bool=true,
  detailedinfo::Bool=info,
  energyConvergenceCriterion::Float64 = 1e-8)

  P = initialGuessDensity
  energyConverged = false
  local totalEnergy = computeEnergyHartreeFock(P,kinetic,nuclearAttraction,coulomb(P),exchange(P),interatomicRepulsion;info=false)
  local energies::Array{Float64,1}
  local P::Matrix
  i = 0
  while (!energyConverged)
    i+=1
    J = coulomb(P)
    K = exchange(P)
    F = kinetic+nuclearAttraction+2J-K
    energies,P = evaluateSCFStep(P,F,overlap,electronNumber)
    oldEnergy = totalEnergy
    totalEnergy = computeEnergyHartreeFock(P,kinetic,nuclearAttraction,J,K,interatomicRepulsion,info=detailedinfo)
    info && println("E[$i]: $totalEnergy")
    if (abs(totalEnergy-oldEnergy) < energyConvergenceCriterion)
      info && println("Converged!")
      energyConverged = true
    end
  end
  return energies, totalEnergy, P
end

function evaluateSCF(
  basis::GaussianBasis,
  geometry::Geometry,
  initialGuessDensity::Matrix,
  electronNumber::Integer;
  detailedinfo::Bool=true,
  info::Bool=true,
  energyConvergenceCriterion::Float64 = 1e-8,
  computeTensorElectronRepulsionIntegralsCoulomb = computeTensorElectronRepulsionIntegrals,
  computeTensorElectronRepulsionIntegralsExchange = computeTensorElectronRepulsionIntegrals)

  S = computeMatrixOverlap(basis)
  T = computeMatrixKinetic(basis)
  V = computeMatrixNuclearAttraction(basis,geometry)
  Vnn = computeEnergyInteratomicRepulsion(geometry)
  ERIsCoulomb = computeTensorElectronRepulsionIntegralsCoulomb(basis)
  ERIsExchange = computeTensorElectronRepulsionIntegralsCoulomb==computeTensorElectronRepulsionIntegralsExchange ? ERIsCoulomb : computeTensorElectronRepulsionIntegralsExchange(basis)

  evaluateSCF(initialGuessDensity, S, T, V, P->computeMatrixCoulomb(ERIsCoulomb,P), P->computeMatrixExchange(ERIsExchange,P), Vnn, electronNumber; energyConvergenceCriterion=energyConvergenceCriterion, info=info, detailedinfo=detailedinfo)
end

function evaluateSCF(
  shells::Vector{Shell},
  geometry::Geometry,
  initialGuessDensity::Matrix,
  electronNumber::Integer;
  detailedinfo::Bool=true,
  info::Bool=true,
  energyConvergenceCriterion::Float64 = 1e-8)

  S = computeMatrixOverlap(shells)
  T = computeMatrixKinetic(shells)
  V = computeMatrixNuclearAttraction(shells,geometry)
  Vnn = computeEnergyInteratomicRepulsion(geometry)

  #too slow because trivial ERI recomputation in every step:
  #evaluateSCF(initialGuessDensity, S, T, V, P->computeMatrixCoulomb(shells,P), P->computeMatrixExchange(shells,P), Vnn, electronNumber; energyConvergenceCriterion=energyConvergenceCriterion, info=info, detailedinfo=detailedinfo)
  #this would be faster but doesn't allow any shell-structure reuse as we might later want (with sparsity): 
  bas = GaussianBasis([])
  for sh in shells
    append!(bas.contractedBFs,expandShell(sh))
  end
  ERIs = computeTensorElectronRepulsionIntegrals(bas)
  evaluateSCF(initialGuessDensity, S, T, V, P->computeMatrixCoulomb(ERIs,P), P->computeMatrixExchange(ERIs,P), Vnn, electronNumber; energyConvergenceCriterion=energyConvergenceCriterion)
end

if (LibInt2Shell != Shell)
  function evaluateSCF(
    shells::Vector{LibInt2Shell},
    geometry::Geometry,
    initialGuessDensity::Matrix,
    electronNumber::Integer;
    detailedinfo::Bool=true,
    info::Bool=true,
    energyConvergenceCriterion::Float64 = 1e-8)
  
    S = computeMatrixOverlap(shells)
    T = computeMatrixKinetic(shells)
    V = computeMatrixNuclearAttraction(shells,geometry)
    Vnn = computeEnergyInteratomicRepulsion(geometry)
  
    evaluateSCF(initialGuessDensity, S, T, V, P->computeMatrixCoulomb(shells,P), P->computeMatrixExchange(shells,P), Vnn, electronNumber; energyConvergenceCriterion=energyConvergenceCriterion, info=info, detailedinfo=detailedinfo)
  end
end

end # module
