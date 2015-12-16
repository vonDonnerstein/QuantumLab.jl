module BasisFunctionsModule
export PrimitiveGaussianBasisFunction, ContractedGaussianBasisFunction

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

end
