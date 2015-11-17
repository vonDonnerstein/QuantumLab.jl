module QuantumLab

using Reexport

include("SubModules/BaseModule.jl")
@reexport using .BaseModule
include("SubModules/AtomModule.jl")
@reexport using .AtomModule 
include("SubModules/GeometryModule.jl")
@reexport using .GeometryModule
include("SubModules/BasisSetModule.jl")
@reexport using .BasisSetModule
include("SubModules/BasisModule.jl")
@reexport using .BasisModule
include("SubModules/IntegralsModule.jl")
@reexport using .IntegralsModule


end # module
