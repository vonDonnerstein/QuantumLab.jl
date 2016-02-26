module BasisFunctionsModule
export PrimitiveGaussianBasisFunction, ContractedGaussianBasisFunction, prettyprint
import Base.display

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
  coefficients::Array{Float64,1}
  primitiveBFs::Array{PrimitiveGaussianBasisFunction,1}
end

function display(cgbf::ContractedGaussianBasisFunction,indent="")
  print(indent); dump(typeof(cgbf))
  for idx in 1:length(cgbf.coefficients)
    print(indent * "  $(cgbf.coefficients[idx]) × "); dump(typeof(cgbf.primitiveBFs[idx]))
    print(indent * "    center:   "); println(cgbf.primitiveBFs[idx].center)
    print(indent * "    exponent: "); println(cgbf.primitiveBFs[idx].exponent)
    print(indent * "    mqn:      "); println(cgbf.primitiveBFs[idx].mqn)
  end
end

end
