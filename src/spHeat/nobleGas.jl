#----------------------------------------------------------------------------------------------#
#                        Noble Gas Specific Heat Model for Ideal Gases                         #
#----------------------------------------------------------------------------------------------#

# Type declaration
struct nobleGasHeat{𝗽,𝘅,𝗯<:IntBase} <: ConstHeat{𝗽,𝘅}
    name::String        # Substance name -- that has the (M, c) values
    form::String        # Substance formula as a String
    M::mAmt{𝗽,𝘅,MO}     # The precision- exactness- parametric molar mass
    c::cpAmt{𝗽,𝘅,𝗯}     # The precision- exactness- base- parametric cp
    Tref::sysT{𝗽,𝘅}     # The reference state temperature
    sref::sAmt{𝗽,𝘅,𝗯}   # The reference state specific entropy
    # Inner copy constructor
    nobleGasHeat(x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
        new{𝗽,𝘅,𝗯}(x.name, x.form, x.M, x.c, x.Tref, x.sref)
    end
    # Inner checking & promoting constructor
    nobleGasHeat(NAM::AbstractString,
                 FOR::AbstractString,
                 __M::mAmt{𝗽𝗔,𝘅𝗔,MO},
                 __c::cpAmt{𝗽𝗕,𝘅𝗕,𝗯},
                 T_r::sysT{𝗽𝗖,𝘅𝗖}   = T(promote_type(𝗽𝗔, 𝗽𝗕), promote_type(𝘅𝗔, 𝘅𝗕)),
                 s_r::sAmt{𝗽𝗗,𝘅𝗗,𝗯} = sAmt{promote_type(𝗽𝗔, 𝗽𝗕),promote_type(𝘅𝗔, 𝘅𝗕),𝗯}(
                                           zero(promote_type(𝗽𝗔, 𝗽𝗕)))
                ) where {𝗽𝗔,𝘅𝗔,𝗽𝗕,𝘅𝗕,𝗽𝗖,𝘅𝗖,𝗽𝗗,𝘅𝗗,𝗯} = begin
        # Precision and Exactness promotion
        𝗽 = promote_type(𝗽𝗔, 𝗽𝗕, 𝗽𝗖, 𝗽𝗗)
        𝘅 = promote_type(𝘅𝗔, 𝘅𝗕, 𝘅𝗖, 𝘅𝗗)
        # Checks
        @assert amt(__M).val >  0.0
        @assert amt(__c).val >  0.0
        @assert amt(T_r).val >  0.0
        @assert amt(s_r).val >= 0.0
        @assert NAM > ""
        @assert FOR > ""
        # Returns
        new{𝗽,𝘅,𝗯}(NAM, FOR, mAmt{𝗽,𝘅}(__M), cpAmt{𝗽,𝘅}(__c),
                             sysT{𝗽,𝘅}(T_r),  sAmt{𝗽,𝘅}(s_r))
    end
end

# Type exporting
export nobleGasHeat

# Type displaying
deco(x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = Symbol("noble-cp(T)")

Base.show(io::IO, x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
    if DEF[:pprint]
        print(io,
            "$(string(deco(x))) for $(x.form): ",
            "($(x.c)) ($(x.M)) ($(x.Tref)) ($(x.sref))"
        )
    else
        Base.show_default(io, x)
    end
end

# Type plain info access functions
"""
`name(x::nobleGasHeat)::String`\n
Returns a particular gas's name for the substance with specific heat modeled by `x`.
"""
name(x::nobleGasHeat)::String = x.name

"""
`form(x::nobleGasHeat)::String`\n
Returns a particular gas's chemical formula for the substance with specific heat modeled by `x`.
"""
form(x::nobleGasHeat)::String = x.form

"""
`(Tref(x::nobleGasHeat{𝗽,𝘅})::sysT{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular gas's reference state temperature for the substance with specific heat
modeled by `x`.
"""
(Tref(x::nobleGasHeat{𝗽,𝘅})::sysT{𝗽,𝘅}) where {𝗽,𝘅} = x.Tref

"""
`(sref(x::nobleGasHeat{𝗽,𝘅,𝗯})::sAmt{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯}`\n
Returns a particular gas's reference state specific entropy for the substance with specific heat
modeled by `x`.
"""
(sref(x::nobleGasHeat{𝗽,𝘅,𝗯})::sAmt{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = x.sref


#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
#                         Basic Ideal Gas Properties from nobleGasHeat                         #
#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                 M: Particular gas molecular mass                 #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: m

# Particular gas molecular mass
"""
`(m(x::nobleGasHeat{𝗽,𝘅})::mAmt{𝗽,𝘅,MO}) where {𝗽,𝘅}`\n
Returns the particular gas molecular mass for the substance with specific heat modeled by `x`
without conversions.
"""
(m(x::nobleGasHeat{𝗽,𝘅})::mAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = x.M


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
(cp(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MO})::cpAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = cp(x.c * x.M)
(cp(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MA})::cpAmt{𝗽,𝘅,MA}) where {𝗽,𝘅} = cp(x.c / x.M)

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


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #     Δu: Particular gas variation of specific internal energy     #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δu(x::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::sysT{𝗽,𝘅},
     Tf::sysT{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::ΔeAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific internal energy in the specified or default
base for the substance with specific heat modeled by `x`, for process with initial and final
temperatures of `Ti` and `Tf`, respectively.
"""
(Δu(x::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::sysT{𝗽,𝘅},
    Tf::sysT{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::ΔeAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    Δe(cv(x, B) * (Tf - Ti))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #            u: Particular gas specific internal energy            #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: u

"""
`(u(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::sysT{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::uAmt{𝗽,𝘅,B})
`\n
"""
(u(x::nobleGasHeat{𝗽,𝘅,𝗯},
   theT::sysT{𝗽,𝘅},
   B::Type{<:IntBase}=DEF[:IB])::uAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    u(Δu(x, Tref(x), theT, B))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #        Δh: Particular gas variation of specific enthalpy         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δh(x::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::sysT{𝗽,𝘅},
     Tf::sysT{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::ΔeAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific enthalpy in the specified or default base for
the substance with specific heat modeled by `x`, for process with initial and final temperatures
of `Ti` and `Tf`, respectively.
"""
(Δh(x::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::sysT{𝗽,𝘅},
    Tf::sysT{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::ΔeAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    Δe(cp(x, B) * (Tf - Ti))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #               h: Particular gas specific enthalpy                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: h

"""
`(h(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::sysT{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::hAmt{𝗽,𝘅,B})
`\n
"""
(h(x::nobleGasHeat{𝗽,𝘅,𝗯},
   theT::sysT{𝗽,𝘅},
   B::Type{<:IntBase}=DEF[:IB])::hAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    h(Δh(x, Tref(x), theT, B))
end




# TODO: s°, Δs°
