#----------------------------------------------------------------------------------------------#
#                                  Ideal Gas Substance Model                                   #
#----------------------------------------------------------------------------------------------#

import Base: show

# TODO: remove this after EngThermBase's new release
import EngThermBase: PREC, EXAC, BASE

# Type declaration
struct idealGas{𝗽<:PREC,𝘅<:EXAC,𝗛<:Heat} <: Substance{𝗽,𝘅}
	name::String        # The substance name
    form::String        # The chemical formula
    heat::𝗛             # The heat capacity model
    # Inner copy constructor
    idealGas(x::idealGas{𝗽,𝘅,𝗛}) where {𝗽,𝘅,𝗛} = begin
        new{𝗽,𝘅,𝗛}(x.name, x.form, x.heat)
    end
    # Inner checking & promoting constructor
    idealGas(NAM::AbstractString,
             FOR::AbstractString,
             CPM::𝗛) where {𝗛<:Heat{𝗽,𝘅}} where {𝗽,𝘅} = begin
        new{𝗽,𝘅,𝗛}(NAM, FOR, CPM)
    end
end

# Type exporting
export idealGas

