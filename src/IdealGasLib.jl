# Module
module IdealGasLib

# Imports
using Reexport
@reexport using EngThermBase

# Interface
include("interface.jl")

# Includes - specific heat models for ideal gases
include("spHeat/nobleGas.jl")

end # module
