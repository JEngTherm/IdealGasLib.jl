#----------------------------------------------------------------------------------------------#
#                                    subs/idealGas-oper.jl                                     #
#----------------------------------------------------------------------------------------------#

#······························································································#
#                                          Inquiring                                           #
#······························································································#

import EngThermBase: precof, exacof

"""
`precof(::Type{𝕋} | x::𝕋) where 𝕋<:Medium{𝕡} where 𝕡 = 𝕡`\n
Returns the precision of the `Medium` subtype or instance as a `DataType`.
"""
precof(::Type{𝕋}) where 𝕋<:Medium{𝕡} where 𝕡 = 𝕡
precof(x::𝕋) where 𝕋<:Medium{𝕡} where 𝕡 = 𝕡

"""
`exacof(::Type{𝕋} | x::𝕋) where 𝕋<:Medium{𝕡} where 𝕡 = 𝕡`\n
Returns the exactness of the `Medium` subtype or instance as a `DataType`.
"""
exacof(::Type{𝕋}) where 𝕋<:Medium{𝕡,𝕩} where {𝕡,𝕩} = 𝕩
exacof(x::𝕋) where 𝕋<:Medium{𝕡,𝕩} where {𝕡,𝕩} = 𝕩


