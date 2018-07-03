using QuantumLab
using Base.Test

if (!isdefined(:indent))
  indent = ""
end

info("$(now())  $(indent)READING:   h2o.xyz -> h2o::Geometry ")
h2o = readGeometryXYZ("h2o.xyz")

# test BasisSetExchange
# as this takes quite some time we primarily only want to do this during continuous integration (travis-ci) and when we haven-t checked this before
if (!isfile("STO-3G.tx93"))
  info("$(now())  $(indent)OBTAINING:   STO-3G.tx93")
  bseEntries = obtainBasisSetExchangeEntries()
  display(bseEntries)
  stoEntry   = computeBasisSetExchangeEntry("sto-3g",bseEntries)[3]
  downloadBasisSetBasisSetExchange(stoEntry,"STO-3G.tx93")
  @test_throws ErrorException computeBasisSetExchangeEntry("NotDefined",bseEntries)
end

info("$(now())  $(indent)READING:   STO-3G.tx93 -> sto3g::BasisSet")
sto3g = readBasisSetTX93("STO-3G.tx93")

info("$(now())  $(indent)COMPUTING:   sto3g,h2o -> bas")
bas = computeBasis(sto3g,h2o)

info("$(now())  $(indent)COMPUTING:   bas,h2o -> matrixOverlap, matrixKinetic, matrixNuclearAttraction, ERIs")
matrixOverlap = computeMatrixOverlap(bas)
matrixKinetic = computeMatrixKinetic(bas)
matrixNuclearAttraction = computeMatrixNuclearAttraction(bas,h2o)
ERIs = computeTensorElectronRepulsionIntegrals(bas)

info("$(now())  $(indent)COMPUTING:   h2o -> matrixSADguess")
matrixSADguess = computeDensityGuessSAD("HF","STO-3G",h2o)

info("$(now())  $(indent)COMPUTING:   sto3g, h2o -> shells, shells_native")
shells = computeBasisShellsLibInt2(sto3g,h2o)
shells_native = computeBasisShells(sto3g,h2o)

info("$(now())  $(indent)COMPUTING:   (HartreeFock)  shells, h2o, matrixSADguess -> density")
density = evaluateSCF(shells,h2o,mean(matrixSADguess),info=false,detailedinfo=false)[3]

info("$(now())  $(indent)COMPUTING:   density, matrixKinetic, matrixNuclearAttraction, ERIs -> matrixFock")
matrixFock = computeMatrixFock(density,matrixKinetic,matrixNuclearAttraction,ERIs)

densityvirt = inv(matrixOverlap) - density

info("$(now())  $(indent)COMPUTING:   shells, density -> S, T, J, K, Fock")
S = computeMatrixOverlap(shells)
T = computeMatrixKinetic(shells)
J = computeMatrixCoulomb(shells,density)
K = computeMatrixExchange(shells,density)
V = computeMatrixNuclearAttraction(shells,h2o)
fock = computeMatrixFock(T,V,J,K)
moenergies = eigvals(Symmetric(fock),Symmetric(S))
