module QuantumLab

using Reexport

print(" + BaseModule........")
include("SubModules/BaseModule.jl")
@reexport using .BaseModule
println("Done.")

print(" + AtomModule........")
include("SubModules/AtomModule.jl")
@reexport using .AtomModule 
println("Done.")

print(" + GeometryModule....")
include("SubModules/GeometryModule.jl")
@reexport using .GeometryModule
println("Done.")

print(" + BasisSetModule....")
include("SubModules/BasisSetModule.jl")
@reexport using .BasisSetModule
println("Done.")

print(" + BasisModule.......")
include("SubModules/BasisModule.jl")
@reexport using .BasisModule
println("Done.")

print(" + IntegralsModule...")
include("SubModules/IntegralsModule.jl")
@reexport using .IntegralsModule
println("Done.")


end # module
