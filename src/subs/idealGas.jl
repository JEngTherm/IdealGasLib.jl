#----------------------------------------------------------------------------------------------#
#                                  Ideal Gas Substance Model                                   #
#----------------------------------------------------------------------------------------------#

import Base: show

# Type declaration
struct idealGas{𝕡,𝕩,ℍ} <: Substance{𝕡,𝕩}
	name::String        # The substance name
    form::String        # The chemical formula
    heat::ℍ             # The heat capacity model
    # Inner copy constructor
    idealGas(x::idealGas{ℍ}) where {ℍ<:Heat{𝕡,𝕩}} where {𝕡,𝕩} = begin
        new{𝕡,𝕩,ℍ}(x.name, x.form, x.heat)
    end
    # Inner checking & promoting constructor
    idealGas(NAM::AbstractString,
             FOR::AbstractString,
             CPM::ℍ) where {ℍ<:Heat{𝕡,𝕩}} where {𝕡,𝕩} = begin
        new{𝕡,𝕩,ℍ}(NAM, FOR, CPM)
    end
end

# Type exporting
export idealGas

# Type displaying
deco(x::idealGas) = Symbol("ideal gas")

Base.show(io::IO, x::idealGas{ℍ}) where {ℍ<:Heat{𝕡,𝕩}} where {𝕡,𝕩} = begin
    if DEF[:pprint]
        print(io,
            "$(x.name) $(string(deco(x))) \"$(x.form)\" ",
            "with $(x.heat)"
        )
    else
        Base.show_default(io, x)
    end
end


#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
#                               Type plain info access functions                               #
#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Parameter-less function calls
for FUN in (:Tref, :Pref, :sref)
    @eval $FUN(x::idealGas) = ($FUN)(x.heat)
end

# Additional parameter function calls
for FUN in (:sref, :rebase)
    @eval $FUN(x::idealGas, args::Any...) = ($FUN)(x.heat, args...)
end

# Thermodynamic function calls
for FUN in (:m_,:R_,:cp,:cv,:ga,:k_,:Δu,:u_,:Δh,:h_,:Δs°,:s°,:ds,:s_,:Pr,:vr,:RT,:Pv,:Z_)
    @eval $FUN(x::idealGas, args::Any...) = ($FUN)(x.heat, args...)
end


#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
#                                   Ideal Gas EoS Functions                                    #
#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

    # !center 64 | frame 64 -f'\#⋅\# ' | center 76
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                        Pressure Functions                        #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Fallback method, with uniform PREC, EXAC:
P_(x::idealGas{𝕡,𝕩}, T::T_amt{𝕡,𝕩}, v::v_amt{𝕡,𝕩,MO})::P_amt{𝕡,𝕩} = RT(x, T, MO) / v

"""
`P_(x::idealGas{𝕡,𝕩}, T::T_amt{𝕢,𝕪}, v::v_amt{𝕣,𝕫})::P_amt{𝕡,𝕩}`\n
Returns the pressure for the ideal gas `x` at specified temperature `T` and specific volume `v`. Contrary to most `julia` methods, the `x::idealGas{𝕡,𝕩}` model sets the return value precision and exactness, `{𝕡,𝕩}` instead of performing data type promotions.
"""
P_(x::idealGas{𝕡,𝕩}, T::T_amt{𝕢,𝕪}, v::v_amt{𝕣,𝕫})::P_amt{𝕡,𝕩} where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    T = T_amt{𝕡,𝕩}(T)
    v = v_(x, v)
    return RT(x, T, MO) / v
end

    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                      Temperature Functions                       #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                         Volume Functions                         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Base standardization methods
v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕡,𝕩,MO})::v_amt{𝕡,𝕩,MO} = v
v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕡,𝕩,MA})::v_amt{𝕡,𝕩,MO} = v * m_(x)


#----------------------------------------------------------------------------------------------#
#                                           Includes                                           #
#----------------------------------------------------------------------------------------------#

include("idealGas-oper.jl")


