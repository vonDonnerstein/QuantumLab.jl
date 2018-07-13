module QuantumLab

lib_path = joinpath(dirname(@__FILE__),"..","deps","usr","lib")
push!(Libdl.DL_LOAD_PATH,lib_path)
try
  Libdl.dlopen("libint2-QuantumLab.so")
  global const libint2_available = true
catch e
  warn("Caught an error when trying to dlopen libint2-QuantumLab.so \n
           You can try to rebuild with Pkg.build(\"QuantumLab\"). \n
	   Alternatively, you can set the symlink by hand if you know where the
	   dynamic library is located. Until then we'll expect you can't access
	   libint2 and continue with a stub to get the code running anyways.
	   This is typical behavior under Windows and leads to a drastic performance decrease.")
  global const libint2_available = false
end



print(" + (ProgressMeter................")
using ProgressMeter
ProgressMeter.printover(STDOUT," + (Reexport...................")
using Reexport
ProgressMeter.printover(STDOUT,"")

print(" + DocumentationModule.........")
include("SubModules/DocumentationModule.jl")
@reexport using .DocumentationModule
println("Done.")

print(" + BaseModule..................")
include("SubModules/BaseModule.jl")
@reexport using .BaseModule
println("Done.")

print(" + MatrixModule................")
include("SubModules/MatrixModule.jl")
@reexport using .MatrixModule
println("Done.")

print(" + AtomModule..................")
include("SubModules/AtomModule.jl")
@reexport using .AtomModule
println("Done.")

print(" + GeometryModule..............")
include("SubModules/GeometryModule.jl")
@reexport using .GeometryModule
println("Done.")

print(" + BasisSetModule..............")
include("SubModules/BasisSetModule.jl")
@reexport using .BasisSetModule
println("Done.")

print(" + BasisFunctionsModule........")
include("SubModules/BasisFunctionsModule.jl")
@reexport using .BasisFunctionsModule
println("Done.")

print(" + ShellModule.................")
include("SubModules/ShellModule.jl")
@reexport using .ShellModule
println("Done.")

print(" + IntegralsModule.............")
include("SubModules/IntegralsModule.jl")
@reexport using .IntegralsModule
println("Done.")

print(" + BasisModule.................")
include("SubModules/BasisModule.jl")
@reexport using .BasisModule
println("Done.")

print(" + LibInt2Module...............")
include("SubModules/LibInt2Module.jl")
@reexport using .LibInt2Module
println("Done.")

print(" + SpecialMatricesModule.......")
include("SubModules/SpecialMatricesModule.jl")
@reexport using .SpecialMatricesModule
println("Done.")

print(" + InitialGuessModule..........")
include("SubModules/InitialGuessModule.jl")
@reexport using .InitialGuessModule
println("Done.")

print(" + HartreeFockModule...........")
include("SubModules/HartreeFockModule.jl")
@reexport using .HartreeFockModule
println("Done.")

print(" + BasisSetExchangeModule......")
include("SubModules/BasisSetExchangeModule.jl")
@reexport using .BasisSetExchangeModule
println("Done.")

print(" + LaplaceModule...............")
include("SubModules/LaplaceModule.jl")
@reexport using .LaplaceModule
println("Done.")

print(" + ResolutionIdentityModule....")
include("SubModules/ResolutionIdentityModule.jl")
@reexport using .ResolutionIdentityModule
println("Done.")



end # module
