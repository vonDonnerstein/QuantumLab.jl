using QuantumLab
using Base.Test

# write your own tests here
@test 1 == 1

# test readGeometryXYZ
@test readGeometryXYZ("h2o.xyz").atoms[2].element.symbol == "O"
@test -1.42 < readGeometryXYZ("h2o.xyz").atoms[3].position.x < -1.41

# test AtomModule
@test Element("C") == Element("C")

# test BaseModule
mqn = MQuantumNumber(0,0,0)
for (mqn in MQuantumNumbers(LQuantumNumber("D"))) end
@test mqn == MQuantumNumber(2,0,0)
@test distance(Position(0,0,1),Position(1,0,0)) == sqrt(2)

# test readBasisSetTX93
@test readBasisSetTX93("STO-3G.tx93").definitions[Element("C")][1].lQuantumNumber.symbol == "S"
@test readBasisSetTX93("STO-3G.tx93").definitions[Element("C")][1].primitives[1].exponent == 71.616837

# test computeBasis
bas = computeBasis(readBasisSetTX93("STO-3G.tx93"),readGeometryXYZ("h2o.xyz"))
@test -0.74 < bas.contractedBFs[3].primitiveBFs[1].center.y < -0.73

# test IntegralsModule
normalize!(bas)
@test_approx_eq 0.0 computeMatrixOverlap(bas)[4,1]
@test_approx_eq_eps 29.00 computeMatrixKinetic(bas)[2,2] 1e-2
@test_approx_eq_eps 0.16175843 IntegralsModule.FIntegral(2,0.3) 1e-8
@test_approx_eq -π IntegralsModule.NuclearAttractionIntegral(PrimitiveGaussianBasisFunction(Position(0,0,0),1.,MQuantumNumber(1,0,0)),PrimitiveGaussianBasisFunction(Position(0,0,0),1.,MQuantumNumber(1,0,0)),Atom(Element("C"),Position(0,0,0)))
