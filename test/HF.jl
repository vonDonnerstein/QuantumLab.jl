basSet = readBasisSetTX93("STO-3G.tx93")
geo = readGeometryXYZ("h2o.xyz")
bas = computeBasis(basSet,geo)

function computeMatrixCoulomb(basis::GaussianBasis, density::Matrix)
  N = length(basis.contractedBFs)

  J = zeros(N,N)
  μIndex = 0
  νIndex = 0
  for (μ in basis.contractedBFs)
    μIndex+=1
    for (ν in basis.contractedBFs)
      νIndex+=1
      ERImatrix = [computeElectronRepulsionIntegral(μ,ν,λ,σ) for λ in basis.contractedBFs, σ in basis.contractedBFs]
      J[μIndex,νIndex] = trace(ERImatrix*density)
    end
  end
  return J
end
