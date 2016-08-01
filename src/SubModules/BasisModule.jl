module BasisModule
export GaussianBasis, computeBasis, computeBasisShells
import Base.display
import ..IntegralsModule.normalize!

using ..BaseModule
using ..BasisFunctionsModule
using ..BasisSetModule
using ..GeometryModule
using ..AtomModule
using ..IntegralsModule
using ..ShellModule


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
					norm = IntegralsModule.computeIntegralOverlap(primitiveBF,primitiveBF)
					push!(contractedBF.coefficients,primitiveDefinition.prefactor/sqrt(norm))
					push!(contractedBF.primitiveBFs,primitiveBF)
				end
				push!(bas.contractedBFs,contractedBF)
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

function display(basis::GaussianBasis)
  dump(typeof(basis))
  indent = "  "

  println(indent * "contractedBFs: ")
  indent = indent * "  "

  for cgbf in basis.contractedBFs
    display(cgbf,indent)
  end
end

function computeBasisShells(basSet::BasisSet,geo::Geometry)
  shells = Shell[]
  for atom in geo.atoms
    for contractedDefinition in basSet.definitions[atom.element]
      exponents = [prim.exponent for prim in contractedDefinition.primitives]
      coefficients = [prim.prefactor for prim in contractedDefinition.primitives]
      sh = Shell(atom.position,contractedDefinition.lQuantumNumber,exponents,coefficients)
      push!(shells,sh)
    end
  end
  return shells
end

end # module
