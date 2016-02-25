module BasisFunctionsModule
export PrimitiveGaussianBasisFunction, ContractedGaussianBasisFunction, prettyprint
import ..BaseModule.prettyprint

using ..BaseModule

type PrimitiveGaussianBasisFunction
  center::Position
  exponent::Float64
  mqn::MQuantumNumber
end

type ContractedGaussianBasisFunction
  coefficients::Array{Float64,1}
  primitiveBFs::Array{PrimitiveGaussianBasisFunction,1}
end

function prettyprint(cgbf::ContractedGaussianBasisFunction,indent="")
  print(indent); dump(typeof(cgbf))
  for idx in 1:length(cgbf.coefficients)
    print(indent * "  $(cgbf.coefficients[idx]) × "); dump(typeof(cgbf.primitiveBFs[idx]))
    print(indent * "    center:   "); println(cgbf.primitiveBFs[idx].center)
    print(indent * "    exponent: "); println(cgbf.primitiveBFs[idx].exponent)
    print(indent * "    mqn:      "); println(cgbf.primitiveBFs[idx].mqn)
  end
end

end
