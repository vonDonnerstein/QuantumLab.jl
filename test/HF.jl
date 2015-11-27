P = [ 0.5000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000 0.0000000;
      0.0000000   1.0593543  -0.2507533   0.0000000   0.0000000   0.0000000 0.0000000; 
      0.0000000  -0.2507533   1.0593543   0.0000000   0.0000000   0.0000000 0.0000000; 
      0.0000000   0.0000000   0.0000000   0.6666667   0.0000000   0.0000000 0.0000000; 
      0.0000000   0.0000000   0.0000000   0.0000000   0.6666667   0.0000000 0.0000000; 
      0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.6666667 0.0000000; 
      0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000 0.5000000 ]

basSet = readBasisSetTX93("STO-3G.tx93")
geo = readGeometryXYZ("h2o.xyz")
bas = computeBasis(basSet,geo)


function computeMatrixCoulomb(basis::GaussianBasis, density::Matrix)
  N = length(basis.contractedBFs)

  J = zeros(N,N)
  μIndex = 0
  for (μ in basis.contractedBFs)
    μIndex+=1
    for (νIndex in 1:μIndex)
      ν = basis.contractedBFs[νIndex]
      ERImatrix = [computeElectronRepulsionIntegral(μ,ν,λ,σ) for λ in basis.contractedBFs, σ in basis.contractedBFs]
      J[μIndex,νIndex] = trace(ERImatrix*density)
      J[νIndex,μIndex] = J[μIndex,νIndex]
    end
  end
  return J
end

function computeMatrixExchange(basis::GaussianBasis, density::Matrix)
  N = length(basis.contractedBFs)

  K = zeros(N,N)
  μIndex = 0
  for (μ in basis.contractedBFs)
    μIndex+=1
    for (νIndex in 1:μIndex)
      ν = basis.contractedBFs[νIndex]
      ERImatrix = [computeElectronRepulsionIntegral(μ,λ,ν,σ) for λ in basis.contractedBFs, σ in basis.contractedBFs]
      K[μIndex,νIndex] = trace(ERImatrix*density)
      K[νIndex,μIndex] = K[μIndex,νIndex]
    end
  end
  return K
end

function computeEnergyInteratomicRepulsion(
  atom1::Atom,
  atom2::Atom)
  q1 = atom1.element.atomicNumber
  q2 = atom2.element.atomicNumber
  A = atom1.position
  B = atom2.position
  return q1*q2/(distance(A,B)^2)
end

function computeEnergyInteratomicRepulsion(geo::Geometry)
  result = 0.::Float64
  for (A in 1:length(geo.atoms))
    for (B in 1:A-1)
      result += computeEnergyInteratomicRepulsion(geo.atoms[A],geo.atoms[B])
    end
  end
  return result
end

function computeEnergyHartreeFock(
  basis::GaussianBasis,
  geometry::Geometry,
  density::Matrix)
  P = density
  T = computeMatrixKinetic(basis)
  V = computeMatrixNuclearAttraction(basis,geometry)
  J = computeMatrixCoulomb(basis,density)
  K = computeMatrixExchange(basis,density)
  G = 2J-K
  return trace((T+V+1/2*G)*P)
end

function computeMatrixFock(
  basis::GaussianBasis,
  geometry::Geometry,
  density::Matrix)
  T = computeMatrixKinetic(basis)
  V = computeMatrixNuclearAttraction(basis,geometry)
  J = computeMatrixCoulomb(basis,density)
  K = computeMatrixExchange(basis,density)
  return T+V+2J-K
end

function evaluateSCFStep(
  basis::GaussianBasis,
  geometry::Geometry,
  initialGuessDensity::Matrix,
  overlap::Matrix,
  electronNumber::Integer)

  N = electronNumber
  S = overlap
  F = computeMatrixFock(basis,geometry,initialGuessDensity)
  eig = eigfact(F,S)
  C = eig.vectors
  moEnergies = eig.values
  Cocc = eig.vectors[:,1:N]
  P = Cocc * Cocc.'

  return (moEnergies,P)
end

function evaluateSCF(
  basis::GaussianBasis,
  geometry::Geometry,
  initialGuessDensity::Matrix,
  electronNumber::Integer)

  N = electronNumber
  normalize!(basis)
  S = computeMatrixOverlap(bas)
  P = initialGuessDensity

  for (i in 1:30)
    energies,P = evaluateSCFStep(basis,geometry,P,S,N)
    println("E[$i]: $(computeEnergyHartreeFock(basis,geometry,P))")
  end
end

      
h2o = Geometry([Atom(Element("O"),Position(0.0, 0.0, 0.04851804)),
		Atom(Element("H"),Position(0.75300223,0.00000000,-0.51923377)),
		Atom(Element("H"),Position(-0.75300223,0.00000000,-0.51923377))])
angstrom2bohr!(h2o)

bfs = computeBasis(basSet,h2o)
normalize!(bfs)
bfs.contractedBFs[3],bfs.contractedBFs[5] = bfs.contractedBFs[5],bfs.contractedBFs[3]


function prettyprint(cgbf::ContractedGaussianBasisFunction,indent="")
  print(indent); dump(typeof(cgbf))
  for idx in 1:length(cgbf.coefficients)
    print(indent * "  $(cgbf.coefficients[idx]) × "); dump(typeof(cgbf.primitiveBFs[idx]))
    print(indent * "    center:   "); println(cgbf.primitiveBFs[idx].center)
    print(indent * "    exponent: "); println(cgbf.primitiveBFs[idx].exponent)
    print(indent * "    mqn:      "); println(cgbf.primitiveBFs[idx].mqn)
  end
end

function prettyprint(basis::GaussianBasis)
  dump(typeof(basis))
  indent = "  "

  println(indent * "contractedBFs: ")
  indent = indent * "  "

  for cgbf in basis.contractedBFs
    prettyprint(cgbf,indent)
  end

end
