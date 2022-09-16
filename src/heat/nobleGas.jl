#----------------------------------------------------------------------------------------------#
#                        Noble Gas Specific Heat Model for Ideal Gases                         #
#----------------------------------------------------------------------------------------------#

import Base: show

# Type declaration
struct nobleGasHeat{𝗽,𝘅,𝗯<:IntBase} <: ConstHeat{𝗽,𝘅}
    M::mAmt{𝗽,𝘅,MO}     # The precision- exactness- parametric molar mass
    c::cpAmt{𝗽,𝘅,𝗯}     # The precision- exactness- base- parametric cp
    Tref::sysT{𝗽,𝘅}     # The reference state temperature
    Pref::sysP{𝗽,𝘅}     # The reference state pressure
    sref::sAmt{𝗽,𝘅,𝗯}   # The reference state specific entropy
    # Inner copy constructor
    nobleGasHeat(x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
        new{𝗽,𝘅,𝗯}(x.M, x.c, x.Tref, x.Pref, x.sref)
    end
    # Inner checking & promoting constructor
    nobleGasHeat(__M::mAmt{𝗽𝗔,𝘅𝗔,MO},
                 __c::cpAmt{𝗽𝗕,𝘅𝗕,𝗯},
                 T_r::sysT{𝗽𝗖,𝘅𝗖}   = T(promote_type(𝗽𝗔, 𝗽𝗕), promote_type(𝘅𝗔, 𝘅𝗕)),
                 P_r::sysP{𝗽𝗗,𝘅𝗗}   = P(promote_type(𝗽𝗔, 𝗽𝗕), promote_type(𝘅𝗔, 𝘅𝗕)),
                 s_r::sAmt{𝗽𝗘,𝘅𝗘,𝗯} = sAmt{promote_type(𝗽𝗔, 𝗽𝗕),promote_type(𝘅𝗔, 𝘅𝗕),𝗯}(
                                           zero(promote_type(𝗽𝗔, 𝗽𝗕)))
                ) where {𝗽𝗔,𝘅𝗔,𝗽𝗕,𝘅𝗕,𝗽𝗖,𝘅𝗖,𝗽𝗗,𝘅𝗗,𝗽𝗘,𝘅𝗘,𝗯} = begin
        # Precision and Exactness promotion
        𝗽 = promote_type(𝗽𝗔, 𝗽𝗕, 𝗽𝗖, 𝗽𝗗, 𝗽𝗘)
        𝘅 = promote_type(𝘅𝗔, 𝘅𝗕, 𝘅𝗖, 𝘅𝗗, 𝘅𝗘)
        # Checks
        @assert amt(__M).val >  0.0
        @assert amt(__c).val >  0.0
        @assert amt(T_r).val >  0.0
        @assert amt(P_r).val >  0.0
        @assert amt(s_r).val >= 0.0
        # Returns
        new{𝗽,𝘅,𝗯}(mAmt{𝗽,𝘅}(__M),
                   cpAmt{𝗽,𝘅}(__c),
                   sysT{𝗽,𝘅}(T_r),
                   sysP{𝗽,𝘅}(P_r),
                   sAmt{𝗽,𝘅}(s_r))
    end
end

# Type exporting
export nobleGasHeat

# Type displaying
deco(x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = Symbol("noble-cp(T)")

Base.show(io::IO, x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
    if DEF[:pprint]
        print(io,
            "$(string(deco(x))): ",
            "($(x.c)) ($(x.M)) ($(x.Tref)) ($(x.Pref)) ($(x.sref))"
        )
    else
        Base.show_default(io, x)
    end
end

# Type plain info access functions

"""
`(Tref(x::nobleGasHeat{𝗽,𝘅})::sysT{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular gas's reference state temperature for the substance with specific heat
modeled by `x`.
"""
(Tref(x::nobleGasHeat{𝗽,𝘅})::sysT{𝗽,𝘅}) where {𝗽,𝘅} = x.Tref

"""
`(Pref(x::nobleGasHeat{𝗽,𝘅})::sysP{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular gas's reference state pressure for the substance with specific heat modeled
by `x`.
"""
(Pref(x::nobleGasHeat{𝗽,𝘅})::sysP{𝗽,𝘅}) where {𝗽,𝘅} = x.Pref

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
`(R(x::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])::RAmt{𝗽,𝘅,B}) where {𝗽,𝘅}`\n
Returns the particular gas constant for the substance with specific heat modeled by `x` in the
default or specified base.
"""
(R(x::nobleGasHeat{𝗽,𝘅}, B::Type{MA})::RAmt{𝗽,𝘅,MA}) where {𝗽,𝘅} = R(𝗽, 𝘅) / x.M
(R(x::nobleGasHeat{𝗽,𝘅}, B::Type{MO})::RAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = R(𝗽, 𝘅)

(R(x::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])::RAmt{𝗽,𝘅,B}) where {𝗽,𝘅} = R(x, B)


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
(cv(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MA})::cvAmt{𝗽,𝘅,MA}) where {𝗽,𝘅} = cv(x.c - R(x, MA))
(cv(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MO})::cvAmt{𝗽,𝘅,MO}) where {𝗽,𝘅} = cv(x.c - R(x, MO))

# Particular gas cv values: w/ base conversion
(cv(x::nobleGasHeat{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::cvAmt{𝗽,𝘅,B}) where {𝗽,𝘅} = begin
    cv(cp(x, B) - R(x, B))
end


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

# Alias
Du = Δu


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #            u: Particular gas specific internal energy            #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: u

"""
`(u(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::sysT{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::uAmt{𝗽,𝘅,B})
`\n
Returns the particular gas specific internal energy in the specified or default
base for the substance with specific heat modeled by `x`, for states with temperature `theT`.
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

# Alias
Dh = Δh


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #               h: Particular gas specific enthalpy                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

import EngThermBase: h

"""
`(h(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::sysT{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::hAmt{𝗽,𝘅,B})
`\n
Returns the particular gas specific enthalpy in the specified or default base for the substance
with specific heat modeled by `x`, for states with temperature `theT`.
"""
(h(x::nobleGasHeat{𝗽,𝘅,𝗯},
   theT::sysT{𝗽,𝘅},
   B::Type{<:IntBase}=DEF[:IB])::hAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    h(Δh(x, Tref(x), theT, B) + R(x, B) * Tref(x))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    # Δs°: Particular gas variation of ideal gas partial spec. entropy #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δs°(x::nobleGasHeat{𝗽,𝘅,𝗯},
      Ti::sysT{𝗽,𝘅},
      Tf::sysT{𝗽,𝘅},
      B::Type{<:IntBase} = DEF[:IB])::ΔsAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in ideal gas partial specific entropy in the specified or
default base for the substance with specific heat modeled by `x`, for process with initial and
final temperatures of `Ti` and `Tf`, respectively.
"""
(Δs°(x::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::sysT{𝗽,𝘅},
     Tf::sysT{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::ΔsAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    Δs(cp(x, B) * log(Tf/Ti))
end

# Alias
Ds0 = Δs°


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #      s°: Particular gas specific ideal gas partial entropy       #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(s°(x::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::sysT{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::hAmt{𝗽,𝘅,B})
`\n
Returns the particular gas specific ideal gas partial entropy in the specified or default base
for the substance with specific heat modeled by `x`, for states with temperature `theT`.
"""
(s°(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::sysT{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::hAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    s(Δs°(x, Tref(x), theT, B))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #         Δs: Particular gas variation of specific entropy         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δs(x::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::sysT{𝗽,𝘅},
     Tf::sysT{𝗽,𝘅},
     Pi::sysP{𝗽,𝘅},
     Pf::sysP{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::ΔsAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific entropy in the specified or default base for
the substance with specific heat modeled by `x`, for process with initial and final temperatures
and pressures of `Ti` and `Tf`, and `Pi` and `Pf`, respectively.
"""
(Δs(x::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::sysT{𝗽,𝘅},
    Tf::sysT{𝗽,𝘅},
    Pi::sysP{𝗽,𝘅},
    Pf::sysP{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::ΔsAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    Δs(cp(x, B) * log(Tf/Ti) - R(x, B) * log(Pf/Pi))
end

"""
`(Δs(x::nobleGasHeat{𝗽,𝘅,𝗯𝗔},
    Ti::sysT{𝗽,𝘅},
    Tf::sysT{𝗽,𝘅},
    vi::vAmt{𝗽,𝘅,𝗯𝗕},
    vf::vAmt{𝗽,𝘅,𝗯𝗕},
    B::Type{<:IntBase} = DEF[:IB])::ΔsAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯𝗔,𝗯𝗕}`\n
Returns the particular gas variation in specific entropy in the specified or default base for
the substance with specific heat modeled by `x`, for process with initial and final temperatures
and specific volumes of `Ti` and `Tf`, and `vi` and `vf`, respectively.
"""
(Δs(x::nobleGasHeat{𝗽,𝘅,𝗯𝗔},
    Ti::sysT{𝗽,𝘅},
    Tf::sysT{𝗽,𝘅},
    vi::vAmt{𝗽,𝘅,𝗯𝗕},
    vf::vAmt{𝗽,𝘅,𝗯𝗕},
    B::Type{<:IntBase} = DEF[:IB])::ΔsAmt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯𝗔,𝗯𝗕} = begin
    Δs(cv(x, B) * log(Tf/Ti) + R(x, B) * log(vf/vi))
end

# Alias
Ds = Δs


#----------------------------------------------------------------------------------------------#
#                                        Alias exports                                         #
#----------------------------------------------------------------------------------------------#

export Du, Dh, Ds0, s0, Ds


