module BasisModule
export GaussianBasis, computeBasis, normalize!, prettyprint
import Base.display

using ..BaseModule
using ..BasisFunctionsModule
using ..BasisSetModule
using ..GeometryModule
using ..AtomModule
using ..IntegralsModule


abstract Basis

"""
A GaussianBasis is a Basis constructed purely from Gaussian type functions (see PrimitiveGaussianBasisFunction).
The list of basis functions defines the basis. It is most easily constructed from a BasisSet and a Geometry by computeBasis().
"""
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
	  norm = IntegralsModule.Overlap(primitiveBF,primitiveBF)
	  append!(contractedBF.coefficients,[primitiveDefinition.prefactor/sqrt(norm)])
	  append!(contractedBF.primitiveBFs,[primitiveBF])
	end
	append!(bas.contractedBFs,[contractedBF])
      end
    end
  end
  return bas
end

function normalize!(basis::GaussianBasis)
  for (cgb in basis.contractedBFs)
    normalize!(cgb)
  end
end

function normalize!(cgb::ContractedGaussianBasisFunction)
  N = Overlap(cgb,cgb)
  scale!(cgb.coefficients,1/sqrt(N))
end

function display(basis::GaussianBasis)
  dump(typeof(basis))
  indent = "  "

  println(indent * "contractedBFs: ")
  indent = indent * "  "

  for cgbf in basis.contractedBFs
    display(cgbf,indent)
  end
end

end # module
