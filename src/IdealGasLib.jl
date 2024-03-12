#----------------------------------------------------------------------------------------------#
#                                        IdealGasLib.jl                                        #
#----------------------------------------------------------------------------------------------#

# Module
module IdealGasLib

# Imports
using Reexport
@reexport using EngThermBase

# Interface
include("interface.jl")

# Includes - specific heat models for ideal gases
include("heat/nobleGas.jl")

# Includes - the ideal gas EoS model
include("subs/idealGas.jl")

# Tiny fluid libraries
include("../lib/heat/nobleGas.jl")

end # module
