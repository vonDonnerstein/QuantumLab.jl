module BasisFunctionsModule
export PrimitiveGaussianBasisFunction, ContractedGaussianBasisFunction

using ..BaseModule

type PrimitiveGaussianBasisFunction
  center::Position
  exponent::Real
  mqn::MQuantumNumber
end

type ContractedGaussianBasisFunction
  coefficients::Array{Real,1}
  primitiveBFs::Array{PrimitiveGaussianBasisFunction,1}
end

end
