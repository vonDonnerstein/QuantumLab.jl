module BasisModule
export GaussianBasis, computeBasis, computeBasisShells

using ..BaseModule
using ..BasisFunctionsModule
using ..BasisSetModule
using ..GeometryModule
using ..AtomModule
using ..IntegralsModule
using ..ShellModule

import Base.display
import Base.normalize!
import Base.convert

abstract type Basis end

"""
A GaussianBasis is a Basis constructed purely from Gaussian type functions (see PrimitiveGaussianBasisFunction).
The list of basis functions defines the basis. It is most easily constructed from a BasisSet and a Geometry by computeBasis().
"""
immutable GaussianBasis <: Basis	# doesn't change much to the user as contractedBFs is mutable anyways
  contractedBFs::Vector{ContractedGaussianBasisFunction}
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

normalize!(basis::GaussianBasis) = normalize!(basis.contractedBFs)

function display(basis::GaussianBasis)
  println(typeof(basis))
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

convert(::Type{GaussianBasis},sh::Shell) = GaussianBasis(expandShell(sh::Shell))
function convert(::Type{GaussianBasis},shells::Vector{Shell})
  bas = GaussianBasis([])
  for sh in shells
    append!(bas.contractedBFs,expandShell(sh))
  end
  return bas
end


end # module
