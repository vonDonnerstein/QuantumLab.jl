using QuantumLab
using Base.Test

# write your own tests here
@test 1 == 1

# test readGeometryXYZ
@test readGeometryXYZ("h2o.xyz").atoms[2].element.symbol == "O"
@test -0.752 < readGeometryXYZ("h2o.xyz").atoms[3].position.x < -0.750

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
@test -0.389 < bas.contractedBFs[3].primitiveBFs[1].center.y < -0.387

# test IntegralsModule
normalize!(bas)
@test 0.30387 < computeMatrixOverlap(bas)[4,1] < 0.30388 # reference value not checked
@test 1133.98 < computeMatrixKinetic(bas)[2,2] < 1134 # reference value not checked
@test 0.16175843 < IntegralsModule.FIntegral(2,0.3) < 0.16175844
@test -0.5235988 < IntegralsModule.NuclearAttractionIntegral(PrimitiveGaussianBasisFunction(Position(0,0,0),1,MQuantumNumber(1,0,0)),PrimitiveGaussianBasisFunction(Position(0,0,0),1,MQuantumNumber(1,0,0),Atom(Element("C"),Position(0,0,0))) < -0.5235987
