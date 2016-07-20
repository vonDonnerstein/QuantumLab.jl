module ShellModule
export Shell, expandShell, computeBasisShells
using ..BaseModule
using ..BasisFunctionsModule
using ..BasisSetModule
using ..GeometryModule
import Base.display 

"""
    Shell(center::Position, lqn::LQuantumNumber, exponents::[Float], coefficients::[Float64])
A shell of basis functions is the collection of basis functions that share a center and radial definition (e.g. 
the px, py and pz orbital that are defined together by their center, exponents and coefficients. Contrary to other codes we also consider a shell to have a unique LQuantumNumber not allowing for shared-sp tricks. Working with shells
instead of the individual basis functions allows for some performance gains in integral evaluation.
"""
type Shell
  center::Position
  lqn::LQuantumNumber
  exponents::Vector{Float64}
  coefficients::Vector{Float64}
end

function display(sh::Shell)
  println("$(sh.lqn.symbol)-Shell   @    $(sh.center)")
  print("  Exponents:    ")
  for exp in sh.exponents
    @printf("  %f",exp) 
  end
  println("")
  print("  Coefficients: ")
  for coeff in sh.coefficients
    @printf("  %f",coeff) 
  end
  println("")
end

"""
Compute a list of ContractedGaussianBasisFunctions which are described within the input shell.
"""
function expandShell(sh::Shell)
  result = ContractedGaussianBasisFunction[]
  for mqn in MQuantumNumbers(sh.lqn)
    cgbf = ContractedGaussianBasisFunction(sh.coefficients,[])
    for exp in sh.exponents
      push!(cgbf.primitiveBFs,PrimitiveGaussianBasisFunction(sh.center,exp,mqn))
    end
    push!(result,cgbf)
  end
  return result
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
