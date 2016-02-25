module HartreeFockModule
export computeEnergyHartreeFock, evaluateSCF
using ..BaseModule
using ..BasisModule
using ..GeometryModule
using ..SpecialMatricesModule

function computeEnergyHartreeFock(
  basis::GaussianBasis,
  geometry::Geometry,
  density::Matrix)
  P = density
  T = computeMatrixKinetic(basis)
  eT = 2*trace(T*P)
  println("=================================")
  println("Kinetic Energy: $eT")
  V = computeMatrixNuclearAttraction(basis,geometry)
  eV = 2*trace(V*P)
  println("Nuclear Attraction Energy: $eV")
  println("Sum 1-Electron Energy: $(eT+eV)")
  println("---------------------------------")
  J = computeMatrixCoulomb(basis,density)
  eJ = 2*trace(J*P)
  println("Coulomb Energy: $eJ")
  K = computeMatrixExchange(basis,density)
  eK = -trace(K*P)
  println("Exchange Energy: $eK")
  println("Total Electronic Energy: $(eT+eV+eJ+eK)")
  println("Total Energy: $(eT+eV+eJ+eK+computeEnergyInteratomicRepulsion(geometry))")
  println("=================================")
  G = 2J-K
  return trace((T+V+1/2*G)*P*2)
end

function evaluateSCFStep(
  basis::GaussianBasis,
  geometry::Geometry,
  initialGuessDensity::Matrix,
  overlap::Matrix,
  electronNumber::Integer)

  N = electronNumber
  S = overlap
  F = computeMatrixFock(basis,geometry,initialGuessDensity)
  eig = eigfact(Symmetric(F),Symmetric(S)) # because Symmetric is faster and allows for sorted eigenvalues
  C = eig.vectors
  moEnergies = eig.values
  Cocc = eig.vectors[:,1:N]
  P = Cocc * Cocc.'

  return (moEnergies,P)
end

function evaluateSCF(
  basis::GaussianBasis,
  geometry::Geometry,
  initialGuessDensity::Matrix,
  electronNumber::Integer;
  energyConvergenceCriterion::Float64 = 1e-8)

  N = electronNumber
  normalize!(basis)
  S = computeMatrixOverlap(basis)
  P = initialGuessDensity

  energyConverged = false
  local totalEnergy = computeEnergyHartreeFock(basis,geometry,P)+computeEnergyInteratomicRepulsion(geometry)
  local energies::Array{Float64,1}
  local P::Matrix
  i = 0
  while (!energyConverged)
    i+=1
    energies,P = evaluateSCFStep(basis,geometry,P,S,N)
    oldEnergy = totalEnergy
    totalEnergy = computeEnergyHartreeFock(basis,geometry,P)+computeEnergyInteratomicRepulsion(geometry)
    println("E[$i]: $totalEnergy")
    if (abs(totalEnergy-oldEnergy) < energyConvergenceCriterion)
      println("Converged!")
      energyConverged = true
    end
  end
  return energies, P
end

end # module
