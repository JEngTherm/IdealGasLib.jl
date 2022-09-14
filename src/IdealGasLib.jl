# Module
module IdealGasLib

# Imports
using Reexport
@reexport using EngThermBase

# Includes - specific heat models for ideal gases
include("spHeat/nobleGas.jl")

end # module
