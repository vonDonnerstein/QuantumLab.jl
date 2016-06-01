using QuantumLab
using Base.Test

# write your own tests here
@test 1 == 1

# test readGeometryXYZ
h2o = readGeometryXYZ("h2o.xyz")
@test h2o.atoms[2].element.symbol == "O"
@test_approx_eq_eps -1.4191843 readGeometryXYZ("h2o.xyz").atoms[3].position.x 1e-7

# test AtomModule
@test Element("C") == Element("C")

# test BaseModule
mqn = MQuantumNumber(0,0,0)
for (mqn in MQuantumNumbers(LQuantumNumber("D"))) end
@test mqn == MQuantumNumber(0,0,2)
@test distance(Position(0,0,1),Position(1,0,0)) == sqrt(2)

# test BasisSetExchange
bseEntries = obtainBasisSetExchangeEntries()
stoEntry   = computeBasisSetExchangeEntry("sto-3g",bseEntries)[3]
downloadBasisSetBasisSetExchange(stoEntry,"STO-3G.tx93")

# test readBasisSetTX93
sto3g = readBasisSetTX93("STO-3G.tx93")
@test sto3g.definitions[Element("C")][1].lQuantumNumber.symbol == "S"
@test sto3g.definitions[Element("C")][1].primitives[1].exponent == 71.616837

run(`rm STO-3G.tx93`)

# test computeBasis
bas = computeBasis(sto3g,h2o)
@test_approx_eq_eps -0.7332137 bas.contractedBFs[3].primitiveBFs[1].center.y 1e-7

# test IntegralsModule
normalize!(bas)
@test_approx_eq_eps 0.3129324434238492 computeMatrixOverlap(bas)[4,1] 1e-8
@test_approx_eq_eps 29.00 computeMatrixKinetic(bas)[2,2] 1e-2
@test_approx_eq_eps 0.16175843 IntegralsModule.FIntegral(2,0.3) 1e-8
@test_approx_eq -π IntegralsModule.NuclearAttractionIntegral(PrimitiveGaussianBasisFunction(Position(0,0,0),1.,MQuantumNumber(1,0,0)),PrimitiveGaussianBasisFunction(Position(0,0,0),1.,MQuantumNumber(1,0,0)),Atom(Element("C"),Position(0,0,0)))

# test InitialGuessModule
@test_approx_eq -0.264689 mean(computeDensityGuessSAD("HF","STO-3G",h2o))[3,2]

# test HartreeFockModule
@test_approx_eq_eps -74.96178985 (computeEnergyHartreeFock(bas,h2o,evaluateSCF(bas,h2o,computeDensityGuessSAD("HF","STO-3G",h2o)[1],5)[3]) + computeEnergyInteratomicRepulsion(h2o)) 1e-7 # checked against FermiONs++

# test ShellModule
shell = LibInt2Shell([0.,0.,0.],0,3,[1.,2.,3.],[.1,.2,.3])

# test LibInt2Module
shells = computeBasisShellsLibInt2(sto3g,h2o)
density = computeDensityGuessSAD("HF","STO-3G",h2o)[1]
@test_approx_eq mean(computeMatrixOverlap(shells)) mean(computeMatrixOverlap(bas))
@test_approx_eq mean(computeMatrixOverlap(shells)[1,4]) mean(computeMatrixOverlap(bas)[1,4])
#@test_approx_eq mean(computeMatrixCoulomb(shells,density)) mean(computeMatrixCoulomb(bas,density))


