#----------------------------------------------------------------------------------------------#
#                                       heat/nobleGas.jl                                       #
#----------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------#
#                        Noble Gas Specific Heat Model for Ideal Gases                         #
#----------------------------------------------------------------------------------------------#

import Base: cp, show # Since :cp is further extended here
import EngThermBase: deco, m_, R_, cv, ga, k_, u_, h_, ds, s_, Pr, vr

# Type declaration
struct nobleGasHeat{𝗽,𝘅,𝗯} <: ConstHeat{𝗽,𝘅,𝗯}
    M::m_amt{𝗽,𝘅,MO}        # The precision- exactness- parametric molar mass
    c::cpamt{𝗽,𝘅,𝗯}         # The precision- exactness- base- parametric cp
    Tref::T_amt{𝗽,𝘅}        # The reference state temperature
    Pref::P_amt{𝗽,𝘅}        # The reference state pressure
    sref::s_amt{𝗽,𝘅,𝗯}      # The reference state specific entropy
    # Inner copy constructor
    nobleGasHeat(𝐻::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
        new{𝗽,𝘅,𝗯}(𝐻.M, 𝐻.c, 𝐻.Tref, 𝐻.Pref, 𝐻.sref)
    end
    # Inner checking & promoting constructor
    nobleGasHeat(__M::m_amt{𝗽𝗔,𝘅𝗔,MO},
                 __c::cpamt{𝗽𝗕,𝘅𝗕,𝗯},
                 T_r::T_amt{𝗽𝗖,𝘅𝗖}   = 𝗧(promote_type(𝗽𝗔, 𝗽𝗕),
                                         promote_type(𝘅𝗔, 𝘅𝗕)),
                 P_r::P_amt{𝗽𝗗,𝘅𝗗}   = 𝗣(promote_type(𝗽𝗔, 𝗽𝗕, 𝗽𝗖),
                                         promote_type(𝘅𝗔, 𝘅𝗕, 𝘅𝗖)),
                 s_r::s_amt{𝗽𝗘,𝘅𝗘,𝗯} = s_amt{promote_type(𝗽𝗔, 𝗽𝗕, 𝗽𝗖, 𝗽𝗗),
                                             promote_type(𝘅𝗔, 𝘅𝗕, 𝘅𝗖, 𝘅𝗗),𝗯}(
                                                zero(promote_type(𝗽𝗔, 𝗽𝗕, 𝗽𝗖, 𝗽𝗗)))
                ) where {𝗽𝗔,𝘅𝗔,𝗽𝗕,𝘅𝗕,𝗽𝗖,𝘅𝗖,𝗽𝗗,𝘅𝗗,𝗽𝗘,𝘅𝗘,𝗯} = begin
        # Precision and Exactness promotion
        𝗽 = promote_type(𝗽𝗔, 𝗽𝗕, 𝗽𝗖, 𝗽𝗗, 𝗽𝗘)
        𝘅 = promote_type(𝘅𝗔, 𝘅𝗕, 𝘅𝗖, 𝘅𝗗, 𝘅𝗘)
        # Checks
        @assert amt(__M).val >= 0.0
        @assert amt(__c).val >= 0.0
        @assert amt(T_r).val >  0.0
        @assert amt(P_r).val >  0.0
        ## @assert amt(s_r).val >= 0.0
        # Returns
        new{𝗽,𝘅,𝗯}(m_amt{𝗽,𝘅}(__M),
                   cpamt{𝗽,𝘅}(__c),
                   T_amt{𝗽,𝘅}(T_r),
                   P_amt{𝗽,𝘅}(P_r),
                   s_amt{𝗽,𝘅}(s_r))
    end
end

# Type exporting
export nobleGasHeat

# Type displaying
deco(𝐻::nobleGasHeat{𝗽,𝘅,MA}) where {𝗽,𝘅} = Symbol("noble-cp(T)")
deco(𝐻::nobleGasHeat{𝗽,𝘅,MO}) where {𝗽,𝘅} = Symbol("noble-c̄p(T)")

Base.show(io::IO, 𝐻::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
    if DEF[:pprint]
        print(io,
            "$(string(deco(𝐻))):\n",
            "   $(𝐻.c)\n    $(𝐻.M)\n    $(𝐻.Tref)\n    $(𝐻.Pref)\n    $(𝐻.sref)"
        )
    else
        Base.show_default(io, 𝐻)
    end
end

# Type plain info access functions

"""
`(Tref(𝐻::nobleGasHeat{𝗽,𝘅})::T_amt{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular gas's reference state temperature for the substance with specific heat
modeled by `𝐻`.
"""
(Tref(𝐻::nobleGasHeat{𝗽,𝘅})::T_amt{𝗽,𝘅}) where {𝗽,𝘅} = 𝐻.Tref

"""
`(Pref(𝐻::nobleGasHeat{𝗽,𝘅})::P_amt{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular gas's reference state pressure for the substance with specific heat modeled
by `𝐻`.
"""
(Pref(𝐻::nobleGasHeat{𝗽,𝘅})::P_amt{𝗽,𝘅}) where {𝗽,𝘅} = 𝐻.Pref

"""
`(sref(𝐻::nobleGasHeat{𝗽,𝘅,𝗯})::s_amt{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯}`\n
Returns a particular gas's reference state specific entropy for the substance with specific heat
modeled by `𝐻`.
"""
(sref(𝐻::nobleGasHeat{𝗽,𝘅,𝗯})::s_amt{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = 𝐻.sref

(sref(𝐻::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MA})::s_amt{𝗽,𝘅,MA}) where {𝗽,𝘅} = 𝐻.sref / 𝐻.M
(sref(𝐻::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MO})::s_amt{𝗽,𝘅,MO}) where {𝗽,𝘅} = 𝐻.sref * 𝐻.M

(sref(𝐻::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MO})::s_amt{𝗽,𝘅,MO}) where {𝗽,𝘅} = 𝐻.sref
(sref(𝐻::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MA})::s_amt{𝗽,𝘅,MA}) where {𝗽,𝘅} = 𝐻.sref


#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
#                         Basic Ideal Gas Properties from nobleGasHeat                         #
#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                 M: Particular gas molecular mass                 #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas molecular mass
"""
`(m_(𝐻::nobleGasHeat{𝗽,𝘅})::m_amt{𝗽,𝘅,MO}) where {𝗽,𝘅}`\n
Returns the particular gas molecular mass for the substance with specific heat modeled by `𝐻`
without conversions.
"""
(m_(𝐻::nobleGasHeat{𝗽,𝘅})::m_amt{𝗽,𝘅,MO}) where {𝗽,𝘅} = 𝐻.M


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                    R: Particular gas constant                    #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas constant -- function syntax thanks to
# https://stackoverflow.com/a/65890762/4038337
"""
`(R_(𝐻::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])::R_amt{𝗽,𝘅,B}) where {𝗽,𝘅}`\n
Returns the particular gas constant for the substance with specific heat modeled by `𝐻` in the
default or specified base.
"""
(R_(𝐻::nobleGasHeat{𝗽,𝘅}, B::Type{MA})::R_amt{𝗽,𝘅,MA}) where {𝗽,𝘅} = R_(𝗽, 𝘅) / 𝐻.M
(R_(𝐻::nobleGasHeat{𝗽,𝘅}, B::Type{MO})::R_amt{𝗽,𝘅,MO}) where {𝗽,𝘅} = R_(𝗽, 𝘅)

# Type stable, fallback version
(R_(𝐻::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])::R_amt{𝗽,𝘅,B}) where {𝗽,𝘅} = R_(𝐻, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              cp: Particular gas iso-P specific heat              #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas cp values: no conversion
"""
`cp(𝐻::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])`\n
Returns the particular gas constant-pressure specific heat in the default or specified base for
the substance with specific heat modeled by `𝐻`, making base conversion only when necessary.
"""
(cp(𝐻::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MA})::cpamt{𝗽,𝘅,MA}) where {𝗽,𝘅} = 𝐻.c
(cp(𝐻::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MO})::cpamt{𝗽,𝘅,MO}) where {𝗽,𝘅} = 𝐻.c

# Particular gas cp values: w/ conversion
(cp(𝐻::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MO})::cpamt{𝗽,𝘅,MO}) where {𝗽,𝘅} = cp(𝐻.c * 𝐻.M)
(cp(𝐻::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MA})::cpamt{𝗽,𝘅,MA}) where {𝗽,𝘅} = cp(𝐻.c / 𝐻.M)

# Type-stable, fallback version
(cp(𝐻::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])::cpamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cp(𝐻, B)

# Temperature specifying methods
(cp(𝐻::nobleGasHeat{𝗽,𝘅},
    T::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::cpamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cp(𝐻, B)
(cp(𝐻::nobleGasHeat{𝗽,𝘅},
    B::Type{<:IntBase},
    T::T_amt{𝗽,𝘅})::cpamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cp(𝐻, B)

# Fallback temperature specifying methods though T-Pairs.
(cp(𝐻::nobleGasHeat{𝗽,𝘅},
    𝒫::hasTPair{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::cpamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cp(𝐻, 𝒫.T, B)
(cp(𝐻::nobleGasHeat{𝗽,𝘅},
    B::Type{<:IntBase},
    𝒫::hasTPair{𝗽,𝘅})::cpamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cp(𝐻, 𝒫.T, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              cv: Particular gas iso-V specific heat              #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas cv values: type-stable, fallback methods
"""
`cv(𝐻::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])`\n
Returns the particular gas constant-volume specific heat in the default or specified base for
the substance with specific heat modeled by `𝐻`, making base conversion only when necessary.
"""
(cv(𝐻::nobleGasHeat{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::cvamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cv(cp(𝐻, B) - R_(𝐻, B))

# Temperature specifying methods
(cv(𝐻::nobleGasHeat{𝗽,𝘅},
    T::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::cvamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cv(𝐻, B)
(cv(𝐻::nobleGasHeat{𝗽,𝘅},
    B::Type{<:IntBase},
    T::T_amt{𝗽,𝘅})::cvamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cv(𝐻, B)

# Fallback temperature specifying methods though T-Pairs.
(cv(𝐻::nobleGasHeat{𝗽,𝘅},
    𝒫::hasTPair{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::cvamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cv(𝐻, 𝒫.T, B)
(cv(𝐻::nobleGasHeat{𝗽,𝘅},
    B::Type{<:IntBase},
    𝒫::hasTPair{𝗽,𝘅})::cvamt{𝗽,𝘅,B}) where {𝗽,𝘅} = cv(𝐻, 𝒫.T, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #             ga: Particular gas specific heat ratio               #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(ga(𝐻::nobleGasHeat{𝗽,𝘅,𝗯})::gaamt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas specific heat ratio for the substance with specific heat modeled by
`𝐻`, without conversions.
"""
(ga(𝐻::nobleGasHeat{𝗽,𝘅,𝗯})::gaamt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = ga(cp(𝐻, 𝗯)/cv(𝐻, 𝗯))

# Temperature specifying method
(ga(𝐻::nobleGasHeat{𝗽,𝘅,𝗯}, T::T_amt{𝗽,𝘅})::gaamt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = ga(𝐻)

# Fallback temperature specifying methods though T-Pairs.
(ga(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝒫::hasTPair{𝗽,𝘅})::gaamt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = ga(𝐻, 𝒫.T)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #         k: Particular gas isentropic expansion exponent          #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(k_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯})::k_amt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas isentropic expansion exponent for the substance with specific heat
modeled by `𝐻`, without conversions. For ideal gases, \$k = ga\$.
"""
(k_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯})::k_amt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = k_(ga(𝐻))  # ga fallback

# Temperature specifying method
(k_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯}, T::T_amt{𝗽,𝘅})::k_amt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = k_(𝐻)

# Fallback temperature specifying methods though T-Pairs.
(k_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝒫::hasTPair{𝗽,𝘅})::k_amt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = k_(𝐻, 𝒫.T)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #     Δu: Particular gas variation of specific internal energy     #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δu(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::T_amt{𝗽,𝘅},
     Tf::T_amt{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific internal energy in the specified or default
base for the substance with specific heat modeled by `𝐻`, for process with initial and final
temperatures of `Ti` and `Tf`, respectively.
"""
(Δu(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    de(cv(𝐻, B) * (Tf - Ti))
end

# Fallback method with hasTPair arguments
(Δu(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝒫i::hasTPair{𝗽,𝘅},
    𝒫f::hasTPair{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = Δu(𝐻, 𝒫i.T, 𝒫f.T, B)

# Alias
du = Δu


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #            u: Particular gas specific internal energy            #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(u_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::u_amt{𝗽,𝘅,B})`\n
Returns the particular gas specific internal energy in the specified or default
base for the substance with specific heat modeled by `𝐻`, for states with temperature `theT`.
"""
(u_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::u_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    u_(Δu(𝐻, Tref(𝐻), theT, B))
end

# Fallback method with hasTPair arguments
(u_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝒫::hasTPair{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::u_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = u_(𝐻, 𝒫.T, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #        Δh: Particular gas variation of specific enthalpy         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δh(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::T_amt{𝗽,𝘅},
     Tf::T_amt{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific enthalpy in the specified or default base for
the substance with specific heat modeled by `𝐻`, for process with initial and final temperatures
of `Ti` and `Tf`, respectively.
"""
(Δh(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    de(cp(𝐻, B) * (Tf - Ti))
end

# Fallback method with hasTPair arguments
(Δh(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝒫i::hasTPair{𝗽,𝘅},
    𝒫f::hasTPair{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = Δh(𝐻, 𝒫i.T, 𝒫f.T, B)

# Alias
dh = Δh


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #               h: Particular gas specific enthalpy                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(h_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::h_amt{𝗽,𝘅,B})`\n
Returns the particular gas specific enthalpy in the specified or default base for the substance
with specific heat modeled by `𝐻`, for states with temperature `theT`.
"""
(h_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::h_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    h_(Δh(𝐻, Tref(𝐻), theT, B) + R_(𝐻, B) * Tref(𝐻))
end

# Fallback method with hasTPair arguments
(h_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝒫::hasTPair{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::h_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = h_(𝐻, 𝒫.T, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    # Δs°: Particular gas variation of ideal gas partial spec. entropy #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δs°(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
      Ti::T_amt{𝗽,𝘅},
      Tf::T_amt{𝗽,𝘅},
      B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in ideal gas partial specific entropy in the specified or
default base for the substance with specific heat modeled by `𝐻`, for process with initial and
final temperatures of `Ti` and `Tf`, respectively.
"""
(Δs°(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::T_amt{𝗽,𝘅},
     Tf::T_amt{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    ds(cp(𝐻, B) * log(Tf/Ti))
end

# Fallback method with hasTPair arguments
(Δs°(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     𝒫i::hasTPair{𝗽,𝘅},
     𝒫f::hasTPair{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = Δs°(𝐻, 𝒫i.T, 𝒫f.T, B)

# Alias
ds0 = Δs°


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #      s°: Particular gas specific ideal gas partial entropy       #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(s°(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B})`\n
Returns the particular gas specific ideal gas partial entropy in the specified or default base
for the substance with specific heat modeled by `𝐻`, for states with temperature `theT`.
"""
(s°(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    s_(Δs°(𝐻, Tref(𝐻), theT, B))
end

# Fallback method with hasTPair arguments
(s°(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝒫::hasTPair{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = s°(𝐻, 𝒫.T, B)

# Alias
s0 = s°


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #         Δs: Particular gas variation of specific entropy         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(ds(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::T_amt{𝗽,𝘅},
     Tf::T_amt{𝗽,𝘅},
     Pi::P_amt{𝗽,𝘅},
     Pf::P_amt{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific entropy in the specified or default base for
the substance with specific heat modeled by `𝐻`, for process with initial and final temperatures
and pressures of `Ti` and `Tf`, and `Pi` and `Pf`, respectively.
"""
(ds(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    Pi::P_amt{𝗽,𝘅},
    Pf::P_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    ds(cp(𝐻, B) * log(Tf/Ti) - R_(𝐻, B) * log(Pf/Pi))
end

(ds(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    Pi::P_amt{𝗽,𝘅},
    Pf::P_amt{𝗽,𝘅},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    ds(𝐻, Ti, Tf, Pi, Pf, B)
end

# Fallback versions with <:EoSPair input types
(ds(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝑖::TPPair{𝗽,𝘅}, # initial (T, P)
    𝑓::TPPair{𝗽,𝘅}, # final (T, P)
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = ds(𝐻, 𝑖.T, 𝑓.T, 𝑖.P, 𝑓.P, B)

"""
`(ds(𝐻::nobleGasHeat{𝗽,𝘅,𝗯𝗔},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    vi::v_amt{𝗽,𝘅,𝗯𝗕},
    vf::v_amt{𝗽,𝘅,𝗯𝗕},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯𝗔,𝗯𝗕}`\n
Returns the particular gas variation in specific entropy in the specified or default base for
the substance with specific heat modeled by `𝐻`, for process with initial and final temperatures
and specific volumes of `Ti` and `Tf`, and `vi` and `vf`, respectively.
"""
(ds(𝐻::nobleGasHeat{𝗽,𝘅,𝗯𝗔},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    vi::v_amt{𝗽,𝘅,𝗯𝗕},
    vf::v_amt{𝗽,𝘅,𝗯𝗕},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯𝗔,𝗯𝗕} = begin
    ds(cv(𝐻, B) * log(Tf/Ti) + R_(𝐻, B) * log(vf/vi))
end

(ds(𝐻::nobleGasHeat{𝗽,𝘅,𝗯𝗔},
    vi::v_amt{𝗽,𝘅,𝗯𝗕},
    vf::v_amt{𝗽,𝘅,𝗯𝗕},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯𝗔,𝗯𝗕} = begin
    ds(𝐻, Ti, Tf, vi, vf, B)    # fallback
end

# Fallback versions with <:EoSPair input types
(ds(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝑖::TvPair{𝗽,𝘅}, # initial (T, v)
    𝑓::TvPair{𝗽,𝘅}, # final (T, v)
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = ds(𝐻, 𝑖.T, 𝑓.T, 𝑖.v, 𝑓.v, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                s: Particular gas specific entropy                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(s_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅},
     theP::P_amt{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B})`\n
Returns the particular gas specific entropy in the specified or default base for the substance
with specific heat modeled by `𝐻`, in the specified thermodynamic state (`theT`, `theP`).
"""
(s_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅},
    theP::P_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    s_(ds(𝐻, Tref(𝐻), theT, Pref(𝐻), theP, B))
end

(s_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    theP::P_amt{𝗽,𝘅},
    theT::T_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    s_(𝐻, theT, theP, B)
end

# Fallback method with TPPair arguments
(s_(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    𝒫::TPPair{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = s_(𝐻, 𝒫.T, 𝒫.P, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #               Pr: Particular gas relative pressure               #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Pr(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅})::Pramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas relative pressure for the substance with specific heat modeled by
`𝐻`, in the specified thermodynamic temperature `theT`.
"""
(Pr(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅})::Pramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = begin
    Pr(exp(s°(𝐻, theT, 𝗯) / R_(𝐻, 𝗯)))
end

# Fallback method with hasTPair arguments
(Pr(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    hasT::hasTPair{𝗽,𝘅})::Pramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = Pr(𝐻, hasT.T)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                vr: Particular gas relative volume                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(vr(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅})::vramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas relative volume for the substance with specific heat modeled by `𝐻`,
in the specified thermodynamic temperature `theT`.
"""
(vr(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅})::vramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = begin
    # The be(𝗽(ℯ)) term is a scale factor to render the numerator dimensionless
    vr(theT * be(𝗽(ℯ)) / Pr(𝐻, theT))
end

# Fallback method with hasTPair arguments
(vr(𝐻::nobleGasHeat{𝗽,𝘅,𝗯},
    hasT::hasTPair{𝗽,𝘅})::vramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = vr(𝐻, hasT.T)



#----------------------------------------------------------------------------------------------#
#                                        Alias exports                                         #
#----------------------------------------------------------------------------------------------#

export du, dh, ds0, s0


#----------------------------------------------------------------------------------------------#
#                                           Includes                                           #
#----------------------------------------------------------------------------------------------#

include("nobleGas-oper.jl")


