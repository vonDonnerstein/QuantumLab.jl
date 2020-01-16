using QuantumLab
using Test
using Dates
using LinearAlgebra

INFO(str) = @info "$(Dates.format(now(),dateformat"Y/m/d HH:MM:SS,sss"))  $str"

if (!@isdefined(indent))
  indent = ""
end

INFO("$(indent)READING:   h2o.xyz -> h2o::Geometry ")
h2o = readGeometryXYZ("h2o.xyz")

INFO("$(indent)OBTAINING:   BasisSetExchange -> sto3g::BasisSet")
sto3g = BasisSet("sto-3g")

INFO("$(indent)COMPUTING:   sto3g,h2o -> bas")
bas = computeBasis(sto3g,h2o)

INFO("$(indent)COMPUTING:   bas,h2o -> matrixOverlap, matrixKinetic, matrixNuclearAttraction, ERIs")
matrixOverlap = computeMatrixOverlap(bas)
matrixKinetic = computeMatrixKinetic(bas)
matrixNuclearAttraction = computeMatrixNuclearAttraction(bas,h2o)
ERIs = computeTensorElectronRepulsionIntegrals(bas)

INFO("$(indent)COMPUTING:   h2o -> matrixSADguess")
matrixSADguess = computeDensityGuessSAD("HF","STO-3G",h2o)

INFO("$(indent)COMPUTING:   sto3g, h2o -> shells, shells_native")
shells = computeBasisShellsLibInt2(sto3g,h2o)
shells_native = computeBasisShells(sto3g,h2o)

INFO("$(indent)COMPUTING:   (HartreeFock)  shells, h2o, matrixSADguess -> density")
density = evaluateSCF(shells,h2o,sum(matrixSADguess)/2,info=false,detailedinfo=false)[3]

INFO("$(indent)COMPUTING:   density, matrixKinetic, matrixNuclearAttraction, ERIs -> matrixFock")
matrixFock = computeMatrixFock(density,matrixKinetic,matrixNuclearAttraction,ERIs)

densityvirt = inv(matrixOverlap) - density

INFO("$(indent)COMPUTING:   shells, density -> S, T, J, K, Fock")
S = computeMatrixOverlap(shells)
T = computeMatrixKinetic(shells)
J = computeMatrixCoulomb(shells,density)
K = computeMatrixExchange(shells,density)
V = computeMatrixNuclearAttraction(shells,h2o)
fock = computeMatrixFock(T,V,J,K)
moenergies = eigvals(Symmetric(fock),Symmetric(S))
