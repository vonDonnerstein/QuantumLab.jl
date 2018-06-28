using QuantumLab
using Base.Test

info("$(now())  SETTING UP VARIABLES FOR FURTHER TESTING")
indent="=> "
include("setup.jl")





# test DocumentationModule
info("$(now())  TESTING:   DocumentationModule")
docModTester() = 1
@doc """
some markdown documentation
""" docModTester
@test typeof(@doc(docModTester)) == Base.Markdown.MD
@add_doc GenericCitation("generic citation") docModTester
@test typeof(@doc(docModTester)) == Vector{DocumentationModule.Documentation}
@add_doc BookCitation(["C. Darwin"], "On the Origin of Species", "978-0451529060") docModTester
@test typeof(@doc(docModTester)) == Vector{DocumentationModule.Documentation}
@add_doc JournalCitation(["M. Mustermann"],"J. Stup. Mistakes",1,12,2016) docModTester
@test typeof(@doc(docModTester)) == Vector{DocumentationModule.Documentation}
@add_doc Citation(JournalCitation(["M. Mustermann"],"J. Stup. Mistakes",1,12,2016)) docModTester
@test typeof(@doc(docModTester)) == Vector{DocumentationModule.Documentation}
@test isa([Base.Markdown.parse("Hellau!"),BookCitation(["C. Darwin"], "On the Origin of Species", "978-0451529060")],Vector{DocumentationModule.Documentation})
@test isa([Base.Markdown.parse("Hellau!"),JournalCitation(["M. Mustermann"],"J. Stup. Mistakes",1,12,2016)],Vector{DocumentationModule.Documentation})
@test isa([Base.Markdown.parse("Hellau!"),DocumentationModule.Documentation(BookCitation(["C. Darwin"], "On the Origin of Species", "978-0451529060"))],Vector{DocumentationModule.Documentation})
@test isa([GenericCitation("Hellau!"),DocumentationModule.Documentation(BookCitation(["C. Darwin"], "On the Origin of Species", "978-0451529060"))],Vector{DocumentationModule.Documentation})

info("$(now())  TESTING:   Basic Functionality")
# test AtomModule
@test Element("C") == Element("C")
el1 = Element("C")
el2 = Element("C")
@test isequal(el1,el2)
@test !isequal(el1,Element("Na"))
# test BaseModule
mqn = MQuantumNumber(0,0,0)
for mqn in MQuantumNumbers(LQuantumNumber("D")) end
@test mqn == MQuantumNumber(0,0,2)
@test distance(Position(0,0,1),Position(1,0,0)) == sqrt(2)
@test length(MQuantumNumbers(LQuantumNumber("P"))) == 3
@test 2*Position(1,2,3) == Position(1,2,3)*2
@test isless(LQuantumNumber("S"),LQuantumNumber("P"))
@test 0. ≈ BaseModule.evaluateFunction(Position(1.,2.,3.), x->x.x+x.y-x.z)
@test h2o.atoms[2].element.symbol == "O"
@test -1.4191843 ≈ readGeometryXYZ("h2o.xyz").atoms[3].position.x atol=1e-7
# test readBasisSetTX93
@test sto3g.definitions[Element("C")][1].lQuantumNumber.symbol == "S"
@test sto3g.definitions[Element("C")][1].primitives[1].exponent == 71.616837
# test BasisModule
@test -0.7332137 ≈ bas.contractedBFs[3].primitiveBFs[1].center.y atol=1e-7
@test 0.35175381 ≈ BasisFunctionsModule.evaluateFunction(origin, bas.contractedBFs[3]) atol=1e-8

# test IntegralsModule
info("$(now())  TESTING:   IntegralsModule")
@test IntegralsModule.GaussianIntegral1D_Valeev(1,.4) ≈ 0.
@test 0.0 ≈ matrixOverlap[6,1]
@test 0.3129324434238492 ≈ matrixOverlap[4,1] atol=1e-8
@test 29.00 ≈ matrixKinetic[2,2] atol=1e-2
@test 0.16175843 ≈ IntegralsModule.FIntegral(2,0.3) atol=1e-8
@test -π ≈ IntegralsModule.computeIntegralNuclearAttraction(PrimitiveGaussianBasisFunction(Position(0,0,0),1.,MQuantumNumber(1,0,0)),PrimitiveGaussianBasisFunction(Position(0,0,0),1.,MQuantumNumber(1,0,0)),Atom(Element("C"),Position(0,0,0)))
@test IntegralsModule.GaussianIntegral1D_Mathematica(6,1.2) ≈ IntegralsModule.GaussianIntegral1D_Valeev(6,1.2)
@test IntegralsModule.OverlapFundamental(bas.contractedBFs[1].primitiveBFs[1],bas.contractedBFs[1].primitiveBFs[2]) ≈ computeIntegralOverlap(bas.contractedBFs[1].primitiveBFs[1],bas.contractedBFs[1].primitiveBFs[2])

# test InitialGuessModule
info("$(now())  TESTING:   InitialGuessModule")
@test -0.264689 ≈ mean(matrixSADguess)[3,2] atol=1e-6
@test_throws ErrorException computeDensityGuessSAD("NotImplemented","STO-3G",h2o)



# test HartreeFockModule
info("$(now())  TESTING:   HartreeFock (bas, h2o, density)")
@test -74.96178985 ≈ (evaluateSCF(bas,h2o,density,info=false,detailedinfo=false)[2] + computeEnergyInteratomicRepulsion(h2o)) atol=1e-7 # checked against FermiONs++
#@test_approx_eq -4.473355520007 mean(HartreeFockModule.evaluateSCFStep(bas,h2o,mean(matrixSADguess),matrixOverlap)[1])
#@test_approx_eq_eps -74.96178985 (computeEnergyHartreeFock(bas,h2o,density) + computeEnergyInteratomicRepulsion(h2o)) 1e-7 # checked against FermiONs++

# test Shells: ShellModule and LibInt2Module
info("$(now())  TESTING:   Shells/LibInt2Shells/Basis equivalence")
shell_native  = Shell(Position(0.,0.,0.),LQuantumNumber("S"),[1.,2.,3.],[.1,.2,.3],renorm=false)
shell_libint2 = LibInt2Shell([0.,0.,0.],0,3,[1.,2.,3.],[.1,.2,.3])
shell_nativefromlibint2 = Shell(shell_libint2)
@test shell_nativefromlibint2.coefficients[2] ≈ 0.41030724 atol=1e-8
@test shell_native.coefficients[2] ≈ .2
@test computeMatrixBlockOverlap(shell_nativefromlibint2,shell_nativefromlibint2) ≈ computeMatrixBlockOverlap(shell_libint2,shell_libint2)
@test computeTensorBlockElectronRepulsionIntegrals(shell_nativefromlibint2,shell_nativefromlibint2,shell_nativefromlibint2,shell_nativefromlibint2) ≈ computeTensorBlockElectronRepulsionIntegrals(shell_libint2,shell_libint2,shell_libint2,shell_libint2)
@test LibInt2Module.nbf(shell_libint2) ≈ 1
import QuantumLab.libint2_available
if (libint2_available)
  @test_throws ErrorException LibInt2Shell([0.,0.,0.],100,1,[.1],[.5])
end

# test LibInt2Module
info("$(now())  TESTING:   LibInt2Module")
libInt2Finalize()
libInt2Initialize()
tmp_lib = LibInt2Shell([0.,1.,2.],1,2,[1.,2.],[0.5,0.5];renorm=false)
tmp = Shell(tmp_lib).coefficients
@test tmp[1] ≈ tmp[2]
destroy!(tmp_lib)
S = computeMatrixOverlap(shells)
@test S ≈ computeMatrixOverlap(shells_native)
@test S ≈ matrixOverlap
T = computeMatrixKinetic(shells)
@test T ≈ computeMatrixKinetic(shells_native)
@test T ≈ matrixKinetic
J = computeMatrixCoulomb(shells,density)
@test J ≈ computeMatrixCoulomb(bas,density) atol=1e-8
@test computeMatrixBlockOverlap(shells[1],shells[1]) ≈ S[1:1,1:1]
@test computeMatrixBlockKinetic(shells[1],shells[1]) ≈ T[1:1,1:1]
@test computeTensorBlockElectronRepulsionIntegrals(shells[1],shells[5],shells[4],shells[4])[1,1,2,2] ≈ computeElectronRepulsionIntegral(bas.contractedBFs[1],bas.contractedBFs[7],bas.contractedBFs[5],bas.contractedBFs[5])

# test LaplaceModule
info("$(now())  TESTING:   LaplaceModule")
if(!isfile("hackbusch_pretables"))
  downloadLaplacePointsHackbusch("hackbusch_pretables")
end
R = transformRangeToIdealLaplace(0.5,3.)[2]
lp = transformLaplacePointFromIdealLaplace( findLaplacePointsHackbuschPretableLarger(15,R,"hackbusch_pretables")[1], 0.5)
@test_throws ErrorException findLaplacePointsHackbuschPretableSmaller(15,R,"hackbusch_pretables")
@test LaplaceModule.computeInverseByLaplaceApproximation(2.3,lp) ≈ 1./2.3 atol=1e-7
@test LaplaceModule.computeRPADenominatorByDoubleLaplace(1.2,2.3,lp) ≈ 1./(1.2^2 + 2.3^2) atol=1e-3
println("Expecting warning below:")
@test findLaplacePointsHackbuschPretableLarger(2,100.,"hackbusch_pretables") == findLaplacePointsHackbuschPretableLarger(2,50.,"hackbusch_pretables")
@test transformLaplacePointFromIdealLaplace( findLaplacePointsHackbuschPretableSmaller(15,transformRangeToIdealLaplace(0.5,6.)[2],"hackbusch_pretables")[1], 0.5) == lp

# test RI-Module
info("$(now())  TESTING:   RIModule")
@test 0.3950752513027109 ≈ mean(computeMatrixExchangeRIK(bas,bas,matrixSADguess[1]))
tensorRICoulomb = computeTensorElectronRepulsionIntegralsRICoulomb(bas,bas)
@test 0.0429373705056905 ≈ mean(tensorRICoulomb)
@test 0.0012436634924295 ≈ tensorRICoulomb[1,3,4,5]
tensorRIOverlap = computeTensorElectronRepulsionIntegralsRIOverlap(bas,bas)
@test 0.0144597945326691 ≈ tensorRIOverlap[1,3,4,5]
@test 0.0537036598506078 ≈ mean(tensorRIOverlap)

# test SpecialMatricesModule
info("$(now())  SpecialMatricesModule")
@test -0.785008186026 ≈ mean(computeMatrixFock(bas,h2o,matrixSADguess[1])) atol=1e-10                  #### !!!! ####






### display functions ###
# These are simply here, so that the test coverage isn't limited by the display functions.
# They do not however test for anything really, so do not trust them. If you do add a *real*
# test for a display function, please add it to the corresponding module block above instead
# of here.
info("$(now())  DISPLAY FUNCTIONS")
display(bas)
display(bas.contractedBFs[3])
display(JournalCitation(["M. Mustermann"],"J. Stup. Mistakes",1,12,2016))
display(GenericCitation("M. Mustermann personal note"))
display(BookCitation(["C. Darwin"], "On the Origin of Species", "978-0451529060"))
display([BookCitation(["C. Darwin"], "On the Origin of Species", "978-0451529060"), JournalCitation(["M. Mustermann"],"J. Stup. Mistakes",1,12,2016)])
display([Base.Markdown.parse("**Why, god, why??**"),BookCitation(["C. Darwin"], "On the Origin of Species", "978-0451529060"), JournalCitation(["M. Mustermann"],"J. Stup. Mistakes",1,12,2016)])
display([JournalCitation(["M. Mustermann"],"J. Stup. Mistakes",1,12,2016),JournalCitation(["M. Mustermann"],"J. Stup. Mistakes",1,12,2016)])
display(@doc(IntegralsModule.GaussianIntegral1D_Valeev))
display(shells[1])
display(shell_native)
display(lp)
summarize(BasisModule.Basis)
summarize(GaussianBasis)



info("\n\n     CONGRATULATIONS!! ALL TESTS COMPLETED!")
