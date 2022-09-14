# Module
module IdealGasLib

# Imports
using Reexport
@reexport using EngThermBase

# Show
include("show.jl")

# Includes - specific heat models for ideal gases
include("spHeat/nobleGas.jl")

end # module
