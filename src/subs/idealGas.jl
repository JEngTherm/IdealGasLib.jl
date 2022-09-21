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

# Type displaying
deco(x::idealGas) = Symbol("ideal gas")

Base.show(io::IO, x::idealGas{𝗽,𝘅,𝗛}) where {𝗽,𝘅,𝗛} = begin
    if DEF[:pprint]
        print(io,
            "$(x.name) $(string(deco(x))) \"$(x.form)\" ",
            "with $(x.heat)"
        )
    else
        Base.show_default(io, x)
    end
end

# Type plain info access functions

"""
`(Tref(GAS::idealGas{𝗽,𝘅})::sysT{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular ideal `GAS`'s reference state temperature.
"""
(Tref(GAS::idealGas{𝗽,𝘅})::sysT{𝗽,𝘅}) where {𝗽,𝘅} = Tref(GAS.heat)      # fallback

"""
`(Pref(GAS::idealGas{𝗽,𝘅})::sysP{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular ideal `GAS`'s reference state pressure.
"""
(Pref(GAS::idealGas{𝗽,𝘅})::sysP{𝗽,𝘅}) where {𝗽,𝘅} = Pref(GAS.heat)      # fallback

"""
`(sref(GAS::idealGas{𝗽,𝘅,𝗯})::sAmt{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯}`\n
Returns a particular ideal `GAS`'s reference state specific entropy in the base `𝗯`.
"""
(sref(GAS::idealGas{𝗽,𝘅,𝗯})::sAmt{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = sref(GAS.heat)      # fallback

