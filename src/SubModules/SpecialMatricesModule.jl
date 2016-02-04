module SpecialMatricesModule
export computeMatrixKinetic, computeMatrixNuclearAttraction, computeMatrixOverlap 
using ..IntegralsModule
using ..BasisModule
using ..Geometry

function computeMatrixKinetic(basis::GaussianBasis)
  return [KineticIntegral(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixOverlap(basis::GaussianBasis)
  return [Overlap(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixNuclearAttraction(basis::GaussianBasis,geo::Geometry)
  return [sum([NuclearAttractionIntegral(cgb1,cgb2,atom) for atom in geo.atoms]) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

end
