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
import Base.convert

abstract type Basis end

"""
A GaussianBasis is a Basis constructed purely from Gaussian type functions (see PrimitiveGaussianBasisFunction).
The list of basis functions defines the basis. It is most easily constructed from a BasisSet and a Geometry by computeBasis().
"""
immutable GaussianBasis <: Basis	# doesn't change much to the user as contractedBFs is mutable anyways
  contractedBFs::Vector{ContractedGaussianBasisFunction}
end

function computeBasis_primitive(basSet::BasisSet,geo::Geometry)
  bas=GaussianBasis([])
  for atom in geo.atoms
    contrDefinitions = basSet.definitions[atom.element]
    for cdef in contrDefinitions
      for mqn in MQuantumNumbers(cdef.lQuantumNumber)
	coeffs = [pdef.prefactor for pdef in cdef.primitives]
	pgbfs  = [PrimitiveGaussianBasisFunction(atom.position,pdef.exponent,mqn) for pdef in cdef.primitives]
	push!(bas.contractedBFs,ContractedGaussianBasisFunction(coeffs,pgbfs))
      end
    end
  end
  return bas
end

function renormprimitives!(basis::GaussianBasis)
  for cgbf in basis.contractedBFs
    for p in eachindex(cgbf.primitiveBFs)
      pgbf = cgbf.primitiveBFs[p]
      primnorm = computeIntegralOverlap(pgbf,pgbf)
      cgbf.coefficients[p] /= sqrt(primnorm)
    end
  end
end

function computeBasis(basSet::BasisSet,geo::Geometry,normalize::Bool=true)
  bas = computeBasis_primitive(basSet,geo)
  if normalize
    renormprimitives!(bas)
    normalize!.(bas.contractedBFs)
  end
  return bas
end

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
