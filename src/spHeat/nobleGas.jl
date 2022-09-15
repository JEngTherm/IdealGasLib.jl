#----------------------------------------------------------------------------------------------#
#                        Noble Gas Specific Heat Model for Ideal Gases                         #
#----------------------------------------------------------------------------------------------#

# Type declaration
struct nobleGasHeat{𝗽,𝘅,𝗯<:IntBase} <: ConstHeat{𝗽,𝘅}
    name::String        # Substance name -- that has the (M, c) values
    form::String        # Substance formula as a String
    M::mAmt{𝗽,𝘅,MO}     # The precision- exactness- parametric molar mass
    c::cpAmt{𝗽,𝘅,𝗯}     # The precision- exactness- base- parametric cp
end

# Type exporting
export nobleGasHeat

# Type displaying
deco(x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = Symbol("noble-cp(T)")

Base.show(io::IO, x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
    if DEF[:pprint]
        print(io,
            "$(string(deco(x))) for $(x.form): ",
            "($(x.c)) ($(x.M))"
        )
    else
        Base.show_default(io, x)
    end
end


#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
#                         Basic Ideal Gas Properties from nobleGasHeat                         #
#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                 M: Particular gas molecular mass                 #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: M

# Particular gas molecular mass
"""
`(M(x::nobleGasHeat{𝗽,𝘅})::mAmt{𝗽,𝘅,MO}) where {𝗽,𝘅}`\n
Returns the particular gas molecular mass for the substance with specific heat modeled by `x`
without conversions.
"""
(M(x::nobleGasHeat{𝗽,𝘅})::mAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = x.M


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                    R: Particular gas constant                    #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: R

# Particular gas constant -- function syntax thanks to
# https://stackoverflow.com/a/65890762/4038337
"""
`(R(x::nobleGasHeat{𝗽,𝘅})::RAmt{𝗽,𝘅,MA}) where {𝗽,𝘅}`\n
Returns the particular gas constant for the substance with specific heat modeled by `x`.
"""
(R(x::nobleGasHeat{𝗽,𝘅})::RAmt{𝗽,𝘅,MA}) where {𝗽,𝘅} = R(𝗽, 𝘅) / x.M


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              cp: Particular gas iso-P specific heat              #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: cp

# Particular gas cp values: no conversion
"""
`cp(x::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])`\n
Returns the particular gas constant-pressure specific heat in the default or specified base for
the substance with specific heat modeled by `x`, making base conversion only when necessary.
"""
(cp(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MA})::cpAmt{𝗽,𝘅,MA}) where {𝗽,𝘅} = x.c
(cp(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MO})::cpAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = x.c

# Particular gas cp values: w/ conversion
(cp(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MO})::cpAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = x.c * x.M
(cp(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MA})::cpAmt{𝗽,𝘅,MA}) where {𝗽,𝘅} = x.c / x.M

# Particular gas cp value: default base fallback
cp(x::nobleGasHeat) = cp(x, DEF[:IB]) # fallback


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              cv: Particular gas iso-V specific heat              #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: cv

# Particular gas cv values: no base conversion
"""
`cv(x::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])`\n
Returns the particular gas constant-volume specific heat in the default or specified base for
the substance with specific heat modeled by `x`, making base conversion only when necessary.
"""
(cv(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MA})::cvAmt{𝗽,𝘅,MA}) where {𝗽,𝘅} = cv(x.c - R(x))
(cv(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MO})::cvAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = cv(x.c - R(𝗽, 𝘅))

# Particular gas cv values: w/ base conversion
(cv(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MO})::cvAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = begin
    cv(x.c * x.M - R(𝗽, 𝘅))
end
(cv(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MA})::cvAmt{𝗽,𝘅,MA}) where {𝗽,𝘅} = begin
    cv((x.c - R(𝗽, 𝘅)) / x.M)
end

# Particular gas cv value: default base fallback
cv(x::nobleGasHeat) = cv(x, DEF[:IB]) # fallback


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              γ: Particular gas specific heat ratio               #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: γ

"""
`(γ(x::nobleGasHeat{𝗽,𝘅,𝗯})::γAmt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas specific heat ratio for the substance with specific heat modeled by
`x`, without conversions.
"""
(γ(x::nobleGasHeat{𝗽,𝘅,𝗯})::γAmt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = γ(cp(x, 𝗯)/cv(x, 𝗯))


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #         k: Particular gas isentropic expansion exponent          #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: k

"""
`(k(x::nobleGasHeat{𝗽,𝘅,𝗯})::kAmt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas isentropic expansion exponent for the substance with specific heat
modeled by `x`, without conversions. For ideal gases, \$k = γ\$.
"""
(k(x::nobleGasHeat{𝗽,𝘅,𝗯})::kAmt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = k(γ(x))  # γ fallback



