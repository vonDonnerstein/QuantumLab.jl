module BasisFunctionsModule
export PrimitiveGaussianBasisFunction, ContractedGaussianBasisFunction, identityPGF, identityCGF
import Base.display
import QuantumLab.BaseModule.evaluateFunction
using ..BaseModule

"""
A PrimitiveGaussianBasisFunction is a function of the form x^a y^b z^c exp(-ζ r^2) centered at a position in space. (a,b,c) here are by definition the cartesian MQuantumNumber tuple.
"""
type PrimitiveGaussianBasisFunction
  center::Position
  exponent::Float64
  mqn::MQuantumNumber
end

"""
A ContractedGaussianBasisFunction is defined as a sum of PrimitiveGaussianBasisFunctions with fixed relative coefficients. Different from other Quantum Chemistry codes we do not require
an additional scaling for the total ContractedGaussianBasisFunction - scaling the whole ContractedGaussianBasisFunction happens by modifying the contraction coefficients accordingly.
"""
type ContractedGaussianBasisFunction
  coefficients::Vector{Float64}
  primitiveBFs::Vector{PrimitiveGaussianBasisFunction}
end

function display(cgbf::ContractedGaussianBasisFunction,indent="")
  print(indent); println(typeof(cgbf))
  for idx in 1:length(cgbf.coefficients)
    print(indent * "  $(cgbf.coefficients[idx]) × "); println(typeof(cgbf.primitiveBFs[idx]))
    print(indent * "    center:   "); println(cgbf.primitiveBFs[idx].center)
    print(indent * "    exponent: "); println(cgbf.primitiveBFs[idx].exponent)
    print(indent * "    mqn:      "); println(cgbf.primitiveBFs[idx].mqn)
  end
end

function evaluateFunction(x::Position,pgbf::PrimitiveGaussianBasisFunction)
  α = pgbf.exponent
  O = pgbf.center
  r = distance(x,O)
  return (x.x - O.x)^(pgbf.mqn.x) * (x.y - O.y)^(pgbf.mqn.y) * (x.z - O.z)^(pgbf.mqn.z) * exp(-α*r^2)
end

function evaluateFunction(x::Position,cgbf::ContractedGaussianBasisFunction)
  result = 0.
  for (coeff,pgbf) in zip(cgbf.coefficients,cgbf.primitiveBFs)
    result += coeff * evaluateFunction(x,pgbf)
  end
  return result
end

const identityPGF = PrimitiveGaussianBasisFunction(origin, 0., MQuantumNumber(0,0,0))
const identityCGF = ContractedGaussianBasisFunction([1.],[identityPGF])


end
