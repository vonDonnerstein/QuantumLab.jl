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
@test_approx_eq 0.0 computeMatrixOverlap(bas)[6,1]
@test_approx_eq_eps 0.3129324434238492 computeMatrixOverlap(bas)[4,1] 1e-8
@test_approx_eq_eps 29.00 computeMatrixKinetic(bas)[2,2] 1e-2
@test_approx_eq_eps 0.16175843 IntegralsModule.FIntegral(2,0.3) 1e-8
@test_approx_eq -π IntegralsModule.NuclearAttractionIntegral(PrimitiveGaussianBasisFunction(Position(0,0,0),1.,MQuantumNumber(1,0,0)),PrimitiveGaussianBasisFunction(Position(0,0,0),1.,MQuantumNumber(1,0,0)),Atom(Element("C"),Position(0,0,0)))

# test InitialGuessModule
@test_approx_eq -0.264689 mean(computeDensityGuessSAD("HF","STO-3G",h2o))[3,2]

# test HartreeFockModule
densityguess = mean(computeDensityGuessSAD("HF","STO-3G",h2o))
density = evaluateSCF(bas,h2o,densityguess,5)[3]
@test_approx_eq_eps -74.96178985 (computeEnergyHartreeFock(bas,h2o,density) + computeEnergyInteratomicRepulsion(h2o)) 1e-7 # checked against FermiONs++

# test Shells: ShellModule and LibInt2Module
shell_native  =        Shell(LQuantumNumber("S"),Position(0.,0.,0.),[1.,2.,3.],[.1,.2,.3])
shell_libint2 = LibInt2Shell([0.,0.,0.],0,3,[1.,2.,3.],[.1,.2,.3])
shell_nativefromlibint2 = Shell(shell_libint2)
@test_approx_eq shell_native.coefficients[2] .2
@test_approx_eq_eps shell_nativefromlibint2.coefficients[2] 0.41030724 1e-8

# test LibInt2Module
shells = computeBasisShellsLibInt2(sto3g,h2o)
@test_approx_eq computeMatrixOverlap(shells) computeMatrixOverlap(bas)
@test_approx_eq computeMatrixKinetic(shells) computeMatrixKinetic(bas)
@test_approx_eq_eps computeMatrixCoulomb(shells,density) computeMatrixCoulomb(bas,density) 1e-8

# test LaplaceModule
if(!isdir("hackbusch_pretables"))
  downloadLaplacePointsHackbusch("hackbusch_pretables")
end
R = transformRangeToIdealLaplace(0.5,3.)[2]
lp = transformLaplacePointFromIdealLaplace( findLaplacePointsHackbuschPretableLarger(7,R,"hackbusch_pretables")[1], 0.5)
@test_approx_eq_eps LaplaceModule.computeInverseByLaplaceApproximation(2.3,lp) 1./2.3 1e-7

