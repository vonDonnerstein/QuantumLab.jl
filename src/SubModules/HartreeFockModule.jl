module HartreeFockModule
export computeEnergyHartreeFock, evaluateSCF, evaluateHartreeFock
using ..BaseModule
using ..BasisModule
using ..BasisSetModule
using ..GeometryModule
using ..SpecialMatricesModule
using ..IntegralsModule
using ..ShellModule
using ..LibInt2Module
using ..InitialGuessModule

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
  fock::Matrix,
  overlap::Matrix,
  electronNumber::Integer)

  eig = eigfact(Symmetric(fock),Symmetric(overlap)) # because Symmetric is faster and allows for sorted eigenvalues
  moEnergies = eig.values
  Cocc = eig.vectors[:,1:electronNumber]
  P = Cocc * Cocc.'

  return (moEnergies,P)
end

computeMatrixCoefficients(fock::Matrix, overlap::Matrix) = eigvecs(Symmetric(fock),Symmetric(overlap))

"""
    evaluateSCF(initialGuessDensity::Matrix, overlap::Matrix, kinetic::Matrix, nuclearAttraction::Matrix,
                coulomb::Function, exchange::Function,
                interatomicRepulsion::Float64, electronNumber::Integer)

All evaluateSCF calls take the following optional named arguments:
 - `info::Bool=true`,
 - `detailedinfo::Bool=info`,
 - `energyConvergenceCriterion::Float64 = 1e-8`

This is the fundamental routine to evaluate the Self-consistent field cycle of
Hartree-Fock theory and returns the converged (orbital energies, total energy,
density). This function call requires each element in its most primitive form
and is intended for use in unusual contexts or for performance reasons. For
more usable versions see below.
"""
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

  P::Matrix = initialGuessDensity
  totalEnergy = computeEnergyHartreeFock(P,kinetic,nuclearAttraction,coulomb(P),exchange(P),interatomicRepulsion;info=false)
  energies = Vector{Float64}()

  energyConverged = false
  i = 0
  J = coulomb(P)
  K = exchange(P)
  while (!energyConverged)
    i+=1
    F = kinetic+nuclearAttraction+2J-K
    energies,P = evaluateSCFStep(F,overlap,electronNumber)
    J = coulomb(P) # By requiring J,K to contain the new density, the variational principle guarantees no energy below the exact energy is obtained
    K = exchange(P) # think: E2 = P II P = tr(P*G[P]) != tr(P*G[P'])
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

"""
    evaluateSCF(basis, geometry, initialGuess, electronNumber)
    
allows to evaluate the Self-consistent field cycle of Hartree-Fock theory in the most
accessible way. Essentially, this function call combines its arguments in such a way
that it can then call upon the more fundamental `evaluateSCF` that takes only primitive
inputs. To this end, the arguments are very flexible:
- `basis`: One of `Basis`, `BasisSet`, Vector{Shell}, Vector{LibInt2Shell}
- `geo`: One of `Geometry`, `String` (the name of a file containing the geometry)
- `initialGuess`: Either a `Matrix` that is used as initial density matrix, or an `InitialGuess` (ZeroGuess, SADGuess)
- `electronNumber`: number of occupied closed shell orbitals (as the other forms are not yet implemented) or `Ionization` of the system

This evaluateSCF method additionally introduces the following optional named arguments:
- `computeTensorElectronRepulsionIntegralsCoulomb`: f: basis/shells/... -> ERIs (used for Coulomb part (and Exchange part if not overridden))
- `computeTensorElectronRepulsionIntegralsExchange`: f: basis/shells/... -> ERIs (used for Exchange part)
"""
function evaluateSCF(
  basis::Union{GaussianBasis,Vector{Shell},Vector{LibInt2Shell},BasisSet},
  geometry::Union{Geometry,String},
  initialGuessDensity::Union{Matrix,InitialGuess},
  electronNumber::Union{Integer,Ionization}=Ionization(0);
  detailedinfo::Bool=true,
  info::Bool=true,
  energyConvergenceCriterion::Float64 = 1e-8,
  computeTensorElectronRepulsionIntegralsCoulomb = computeTensorElectronRepulsionIntegrals,
  computeTensorElectronRepulsionIntegralsExchange = computeTensorElectronRepulsionIntegrals)

  if isa(geometry,String)
    geometry = Geometry(geometry)
  end

  input = basis
  if isa(input,Vector{LibInt2Shell}) && (LibInt2Shell != Shell)
    # keep working with LibInt2Shell
  elseif isa(input,BasisSet) && (LibInt2Shell != Shell)
    basis = computeBasisShellsLibInt2(input,geometry)
  elseif isa(input,BasisSet) && (LibInt2Shell == Shell)
    basis = computeBasis(input,geometry)
  elseif isa(input,GaussianBasis)
    # At the moment we simply use implicit Shell->Basis conversion which avoids ERI recomputation in every step and is therefore faster for the moment
    # later we will want to switch over to working with Shells directly as shell-structure allows for sparsity
    # evaluateSCF(Pinit, S, T, V, P->computeMatrixCoulomb(shells,P), P->computeMatrixExchange(shells,P), Vnn, electronNumber; energyConvergenceCriterion=energyConvergenceCriterion, info=info, detailedinfo=detailedinfo)
    basis = convert(GaussianBasis,input)
  end

  if isa(electronNumber,Ionization)
    electronNumber::Int64 = computeNumberElectrons(geometry,electronNumber.charge)/2
  end

  S = computeMatrixOverlap(basis)
  T = computeMatrixKinetic(basis)
  V = computeMatrixNuclearAttraction(basis,geometry)
  Pinit = InitialGuessModule.computeDensityMatrixGuess(initialGuessDensity,size(S)[1])
  Vnn = computeEnergyInteratomicRepulsion(geometry)

  if (LibInt2Shell != Shell) && isa(basis,Vector{LibInt2Shell})
    Xarg = basis
    Carg = basis
  else
    ERIsCoulomb = computeTensorElectronRepulsionIntegralsCoulomb(basis)
    ERIsExchange = computeTensorElectronRepulsionIntegralsCoulomb==computeTensorElectronRepulsionIntegralsExchange ? ERIsCoulomb : computeTensorElectronRepulsionIntegralsExchange(basis)
    Xarg = ERIsExchange
    Carg = ERIsCoulomb
  end

  evaluateSCF(Pinit, S, T, V, P->computeMatrixCoulomb(Carg,P), P->computeMatrixExchange(Xarg,P), Vnn, electronNumber; energyConvergenceCriterion=energyConvergenceCriterion, info=info, detailedinfo=detailedinfo)
end

evaluateHartreeFock = evaluateSCF

end # module
