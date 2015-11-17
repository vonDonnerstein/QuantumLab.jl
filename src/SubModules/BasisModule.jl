module BasisModule
export GaussianBasis, PrimitiveGaussianBasisFunction, ContractedGaussianBasisFunction, computeBasis

using ..BaseModule
using ..BasisSetModule
using ..GeometryModule
using ..AtomModule

abstract Basis

type PrimitiveGaussianBasisFunction
  center::Position
  exponent::Real
  mqn::MQuantumNumber
end

type ContractedGaussianBasisFunction
  coefficients::Array{Real,1}
  primitiveBFs::Array{PrimitiveGaussianBasisFunction,1}
end

type GaussianBasis <: Basis
  contractedBFs::Array{ContractedGaussianBasisFunction,1}
end

function computeBasis(basSet::BasisSet,geo::Geometry)
  bas=GaussianBasis([])
  for atom in geo.atoms
    for contractedDefinition in basSet.definitions[atom.element]
      for mqn in MQuantumNumbers(contractedDefinition.lQuantumNumber)
	contractedBF = ContractedGaussianBasisFunction([],[])
	for primitiveDefinition in contractedDefinition.primitives
	  exponent=primitiveDefinition.exponent
	  primitiveBF = PrimitiveGaussianBasisFunction(atom.position,exponent,mqn)
	  append!(contractedBF.coefficients,[primitiveDefinition.prefactor])
	  append!(contractedBF.primitiveBFs,[primitiveBF])
	end
	append!(bas.contractedBFs,[contractedBF])
      end
    end
  end
  return bas
end


end # module
