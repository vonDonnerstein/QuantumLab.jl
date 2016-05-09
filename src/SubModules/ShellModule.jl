module ShellModule
export Shell
using ..BaseModule
import Base.display

"""
A shell of basis functions is the collection of basis functions that share a center and radial definition (e.g. 
the px, py and pz orbital that are defined together by their center, exponents and coefficients. Contrary to other codes we also consider a shell to have a unique LQuantumNumber not allowing for shared-sp tricks. Working with shells
instead of the individual basis functions allows for some performance gains in integral evaluation.
"""
type Shell
  lqn::LQuantumNumber
  center::Position
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

end # module
