using Test, IdealGasLib

# Module-level tests
include("IdealGasLib.test.jl")

# Included sources test
include("interface.test.jl")

# IdealGasLib specific heat model tests
# include("heat/nobleGas.test.jl")

# idealGasLib substance model tests
# include("subs/idealGas.test.jl")

