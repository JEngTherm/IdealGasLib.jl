#----------------------------------------------------------------------------------------------#
#                                       heat/nobleGas.jl                                       #
#----------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------#
#                        Noble Gas Specific Heat Model for Ideal Gases                         #
#----------------------------------------------------------------------------------------------#

import Base: cp, show # Since :cp is further extended here
import EngThermBase: deco, m_, R_, cv, ga, k_, u_, h_, ds, s_, Pr, vr, Pv, RT, Z_

# Type declaration
struct nobleGasHeat{𝕡,𝕩} <: ConstHeat{𝕡,𝕩}
    M::m_amt{𝕡,𝕩,MO}        # The precision- exactness- parametric molar mass
    c::cpamt{𝕡,𝕩,MO}        # The precision- exactness- base- parametric cp
    Tref::T_amt{𝕡,𝕩}        # The reference state temperature
    Pref::P_amt{𝕡,𝕩}        # The reference state pressure
    sref::s_amt{𝕡,𝕩,MO}     # The reference state specific entropy
    # Inner copy constructor
    nobleGasHeat(𝐻::nobleGasHeat{𝕡,𝕩}) where {𝕡,𝕩} = begin
        new{𝕡,𝕩}(𝐻.M, 𝐻.c, 𝐻.Tref, 𝐻.Pref, 𝐻.sref)
    end
    # Inner checking & promoting constructor
    nobleGasHeat(__M::m_amt{𝕡𝔸,𝕩𝔸,MO},
                 __c::cpamt{𝕡𝔹,𝕩𝔹,MO},
                 T_r::T_amt{𝕡ℂ,𝕩ℂ}    = T_(promote_type(𝕡𝔸, 𝕡𝔹),
                                           promote_type(𝕩𝔸, 𝕩𝔹)),
                 P_r::P_amt{𝕡𝔻,𝕩𝔻}    = P_(promote_type(𝕡𝔸, 𝕡𝔹, 𝕡ℂ),
                                           promote_type(𝕩𝔸, 𝕩𝔹, 𝕩ℂ)),
                 s_r::s_amt{𝕡𝔼,𝕩𝔼,MO} = s_amt{promote_type(𝕡𝔸, 𝕡𝔹, 𝕡ℂ, 𝕡𝔻),
                                              promote_type(𝕩𝔸, 𝕩𝔹, 𝕩ℂ, 𝕩𝔻),MO}(
                                                zero(promote_type(𝕡𝔸, 𝕡𝔹, 𝕡ℂ, 𝕡𝔻)))
                ) where {𝕡𝔸,𝕩𝔸,𝕡𝔹,𝕩𝔹,𝕡ℂ,𝕩ℂ,𝕡𝔻,𝕩𝔻,𝕡𝔼,𝕩𝔼} = begin
        # Precision and Exactness promotion
        𝕡 = promote_type(𝕡𝔸, 𝕡𝔹, 𝕡ℂ, 𝕡𝔻, 𝕡𝔼)
        𝕩 = promote_type(𝕩𝔸, 𝕩𝔹, 𝕩ℂ, 𝕩𝔻, 𝕩𝔼)
        # Checks
        @assert amt(__M).val >= 0.0
        @assert amt(__c).val >= 0.0
        @assert amt(T_r).val >  0.0
        @assert amt(P_r).val >  0.0
        ## @assert amt(s_r).val >= 0.0
        # Returns
        new{𝕡,𝕩}(m_amt{𝕡,𝕩}(__M),
                 cpamt{𝕡,𝕩}(__c),
                 T_amt{𝕡,𝕩}(T_r),
                 P_amt{𝕡,𝕩}(P_r),
                 s_amt{𝕡,𝕩}(s_r))
    end
end

# Type exporting
export nobleGasHeat

# Type displaying
deco(𝐻::nobleGasHeat{𝕡,𝕩}) where {𝕡,𝕩} = Symbol("noble-cp(T)")

Base.show(io::IO, 𝐻::nobleGasHeat{𝕡,𝕩}) where {𝕡,𝕩} = begin
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
`(Tref(𝐻::nobleGasHeat{𝕡,𝕩})::T_amt{𝕡,𝕩}) where {𝕡,𝕩}`\n
Returns a particular gas's reference state temperature for the substance with specific heat
modeled by `𝐻`.
"""
(Tref(𝐻::nobleGasHeat{𝕡,𝕩})::T_amt{𝕡,𝕩}) where {𝕡,𝕩} = 𝐻.Tref

"""
`(Pref(𝐻::nobleGasHeat{𝕡,𝕩})::P_amt{𝕡,𝕩}) where {𝕡,𝕩}`\n
Returns a particular gas's reference state pressure for the substance with specific heat modeled
by `𝐻`.
"""
(Pref(𝐻::nobleGasHeat{𝕡,𝕩})::P_amt{𝕡,𝕩}) where {𝕡,𝕩} = 𝐻.Pref

"""
`(sref(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕩}`\n
Returns a particular gas's reference state specific entropy for the substance with specific heat
modeled by `𝐻`.
"""
(sref(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{MA})::s_amt{𝕡,𝕩,MA}) where {𝕡,𝕩} = 𝐻.sref / 𝐻.M
(sref(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{MO})::s_amt{𝕡,𝕩,MO}) where {𝕡,𝕩} = 𝐻.sref

# Type stable, fallback version
(sref(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    sref(𝐻, B)
end


#······························································································#
#                                           Rebasing                                           #
#······························································································#

"""
`(rebase(𝐻::nobleGasHeat{𝕡,𝕩}, 𝑇::T_amt{𝕡,𝕩}, 𝑃::P_amt{𝕡,𝕩})::nobleGasHeat{𝕡,𝕩}) where {𝕡,𝕩}`\n
Returns a `nobleGasHeat` instance based on `𝐻` with `(Tref, Pref) = (𝑇, 𝑃)`, and with `sref`
adjusted so as to yield same entropy values for the same `(T, P)` states than `𝐻`. Values of
`s°` will also coincide only if `𝐻.Pref == 𝑃`.
"""
(rebase(𝐻::nobleGasHeat{𝕡,𝕩},
        𝑇::T_amt{𝕡,𝕩},
        𝑃::P_amt{𝕡,𝕩})::nobleGasHeat{𝕡,𝕩}) where {𝕡,𝕩} = begin
    nobleGasHeat(𝐻.M, 𝐻.c, 𝑇, 𝑃, s_(𝐻, 𝑇, 𝑃))
end

# Fallback versions
(rebase(𝐻::nobleGasHeat{𝕡,𝕩},
        𝑃::P_amt{𝕡,𝕩},
        𝑇::T_amt{𝕡,𝕩})::nobleGasHeat{𝕡,𝕩}) where {𝕡,𝕩} = rebase(𝐻, 𝑇, 𝑃)

(rebase(𝐻::nobleGasHeat{𝕡,𝕩},
        þ::TPPair{𝕡,𝕩})::nobleGasHeat{𝕡,𝕩}) where {𝕡,𝕩} = rebase(𝐻, þ.T, þ.P)

export rebase


#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
#                         Basic Ideal Gas Properties from nobleGasHeat                         #
#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                 M: Particular gas molecular mass                 #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas molecular mass
"""
`(m_(𝐻::nobleGasHeat{𝕡,𝕩})::m_amt{𝕡,𝕩,MO}) where {𝕡,𝕩}`\n
Returns the particular gas molecular mass for the substance with specific heat modeled by `𝐻`
without conversions.
"""
(m_(𝐻::nobleGasHeat{𝕡,𝕩})::m_amt{𝕡,𝕩,MO}) where {𝕡,𝕩} = 𝐻.M


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                    R: Particular gas constant                    #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas constant -- function syntax thanks to
# https://stackoverflow.com/a/65890762/4038337
"""
`(R_(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])::R_amt{𝕡,𝕩,B}) where {𝕡,𝕩}`\n
Returns the particular gas constant for the substance with specific heat modeled by `𝐻` in the
default or specified base.
"""
(R_(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{MA})::R_amt{𝕡,𝕩,MA}) where {𝕡,𝕩} = R_(𝕡, 𝕩) / 𝐻.M
(R_(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{MO})::R_amt{𝕡,𝕩,MO}) where {𝕡,𝕩} = R_(𝕡, 𝕩)

# Type stable, fallback version
(R_(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])::R_amt{𝕡,𝕩,B}) where {𝕡,𝕩} = R_(𝐻, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              cp: Particular gas iso-P specific heat              #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas cp values: no conversion
"""
`cp(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])`\n
Returns the particular gas constant-pressure specific heat in the default or specified base for
the substance with specific heat modeled by `𝐻`, making base conversion only when necessary.
"""
(cp(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{MO})::cpamt{𝕡,𝕩,MO}) where {𝕡,𝕩} = 𝐻.c
(cp(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{MA})::cpamt{𝕡,𝕩,MA}) where {𝕡,𝕩} = cp(𝐻.c / 𝐻.M)

# Type-stable, fallback version
(cp(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])::cpamt{𝕡,𝕩,B}) where {𝕡,𝕩} = cp(𝐻, B)

# Temperature specifying methods
(cp(𝐻::nobleGasHeat{𝕡,𝕩},
    T::T_amt{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::cpamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = cp(𝐻, B)

# Fallback temperature specifying methods though T-Combos (Pairs/Trios).
(cp(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒯::hasT{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::cpamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = cp(𝐻, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              cv: Particular gas iso-V specific heat              #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas cv values: type-stable, fallback methods
"""
`cv(𝐻::nobleGasHeat{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])`\n
Returns the particular gas constant-volume specific heat in the default or specified base for
the substance with specific heat modeled by `𝐻`, making base conversion only when necessary.
"""
(cv(𝐻::nobleGasHeat{𝕡,𝕩},
    B::Type{<:IntBase} = DEF[:IB])::cvamt{𝕡,𝕩,B}) where {𝕡,𝕩} = cv(cp(𝐻, B) - R_(𝐻, B))

# Temperature specifying methods
(cv(𝐻::nobleGasHeat{𝕡,𝕩},
    T::T_amt{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::cvamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = cv(𝐻, B)

# Fallback temperature specifying methods though T-Pairs.
(cv(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒯::hasT{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::cvamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = cv(𝐻, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #             ga: Particular gas specific heat ratio               #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(ga(𝐻::nobleGasHeat{𝕡,𝕩})::gaamt{𝕡,𝕩}) where {𝕡,𝕩}`\n
Returns the particular gas specific heat ratio for the substance with specific heat modeled by
`𝐻`, without conversions.
"""
(ga(𝐻::nobleGasHeat{𝕡,𝕩})::gaamt{𝕡,𝕩}) where {𝕡,𝕩} = ga(cp(𝐻)/cv(𝐻))

# Temperature specifying method
(ga(𝐻::nobleGasHeat{𝕡,𝕩},
    T::T_amt{𝕢,𝕪})::gaamt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = ga(𝐻)

# Fallback temperature specifying methods though T-Pairs.
(ga(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒯::hasT{𝕢,𝕪})::gaamt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = ga(𝐻)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #         k: Particular gas isentropic expansion exponent          #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(k_(𝐻::nobleGasHeat{𝕡,𝕩})::k_amt{𝕡,𝕩}) where {𝕡,𝕩}`\n
Returns the particular gas isentropic expansion exponent for the substance with specific heat
modeled by `𝐻`, without conversions. For ideal gases, \$k = ga\$.
"""
(k_(𝐻::nobleGasHeat{𝕡,𝕩})::k_amt{𝕡,𝕩}) where {𝕡,𝕩} = k_(ga(𝐻))  # ga fallback

# Temperature specifying method
(k_(𝐻::nobleGasHeat{𝕡,𝕩},
    T::T_amt{𝕢,𝕪})::k_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = k_(𝐻)

# Fallback temperature specifying methods though T-Pairs.
(k_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒯::hasT{𝕢,𝕪})::k_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = k_(𝐻)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #     Δu: Particular gas variation of specific internal energy     #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback method
(Δu(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒾::T_amt{𝕡,𝕩},
    𝒻::T_amt{𝕡,𝕩},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    de(cv(𝐻, B) * (𝒻 - 𝒾))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(Δu(𝐻::nobleGasHeat{𝕡,𝕩},
     𝒾::T_amt{𝕢,𝕪},
     𝒻::T_amt{𝕣,𝕫},
     B::Type{<:IntBase} = DEF[:IB])::deamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫}`\n
Returns the particular gas variation in specific internal energy in the specified or default
base for the substance with specific heat modeled by `𝐻`, for process with initial and final
temperatures of `i` and `f`, respectively. Resulting precision, PREC, and exactness, EXAC, are
model-driven, and not promotion-driven.
"""
(Δu(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒾::T_amt{𝕢,𝕪},
    𝒻::T_amt{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    𝑖 = T_amt{𝕡,𝕩}(𝑖)
    𝑓 = T_amt{𝕡,𝕩}(𝑓)
    return Δu(𝐻, 𝑖, 𝑓, B)
end

# Fallback method with hasTPair arguments
(Δu(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒾::hasT{𝕢,𝕪},
    𝒻::hasT{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = Δu(𝐻, 𝒾.T, 𝒻.T, B)

# Alias
du = Δu


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #            u: Particular gas specific internal energy            #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback method
(u_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕡,𝕩},
    B::Type{<:IntBase}=DEF[:IB])::u_amt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    u_(Δu(𝐻, Tref(𝐻), 𝑇, B))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(u_(𝐻::nobleGasHeat{𝕡,𝕩},
     𝑇::T_amt{𝕢,𝕪},
     B::Type{<:IntBase}=DEF[:IB])::u_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪}`\n
Returns the particular gas specific internal energy in the specified or default base for the
substance with specific heat modeled by `𝐻`, for states with temperature `𝑇`. Resulting
precision, PREC, and exactness, EXAC, are model-driven, and not promotion-driven.
"""
(u_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕢,𝕪},
    B::Type{<:IntBase}=DEF[:IB])::u_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = begin
    𝑇 = T_amt{𝕡,𝕩}(𝑇)
    return u_(𝐻, 𝑇, B)
end

# Fallback method with hasTPair arguments
(u_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒯::hasT{𝕢,𝕪},
    B::Type{<:IntBase}=DEF[:IB])::u_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = u_(𝐻, 𝒯.T, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #        Δh: Particular gas variation of specific enthalpy         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback method
(Δh(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒾::T_amt{𝕡,𝕩},
    𝒻::T_amt{𝕡,𝕩},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    de(cp(𝐻, B) * (𝒻 - 𝒾))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(Δh(𝐻::nobleGasHeat{𝕡,𝕩},
     𝒾::T_amt{𝕢,𝕪},
     𝒻::T_amt{𝕣,𝕫},
     B::Type{<:IntBase} = DEF[:IB])::deamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫}`\n
Returns the particular gas variation in specific enthalpy in the specified or default base for
the substance with specific heat modeled by `𝐻`, for process with initial and final temperatures
of `𝒾` and `𝒻`, respectively. Resulting precision, PREC, and exactness, EXAC, are model-driven,
and not promotion-driven.
"""
(Δh(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒾::T_amt{𝕢,𝕪},
    𝒻::T_amt{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    𝑖 = T_amt{𝕡,𝕩}(𝑖)
    𝑓 = T_amt{𝕡,𝕩}(𝑓)
    return Δh(𝐻, 𝑖, 𝑓, B)
end

# Fallback method with hasTPair arguments
(Δh(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒾::hasT{𝕢,𝕪},
    𝒻::hasT{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = Δh(𝐻, 𝒾.T, 𝒻.T, B)

# Alias
dh = Δh


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #               h: Particular gas specific enthalpy                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback method
(h_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝒯::T_amt{𝕡,𝕩},
    B::Type{<:IntBase}=DEF[:IB])::h_amt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    h_(Δh(𝐻, Tref(𝐻), 𝒯, B) + R_(𝐻, B) * Tref(𝐻))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(h_(𝐻::nobleGasHeat{𝕡,𝕩},
     𝑇::T_amt{𝕢,𝕪},
     B::Type{<:IntBase}=DEF[:IB])::h_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪}`\n
Returns the particular gas specific enthalpy in the specified or default base for the substance
with specific heat modeled by `𝐻`, for states with temperature `𝑇`. Resulting precision, PREC,
and exactness, EXAC, are model-driven, and not promotion-driven.
"""
(h_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕢,𝕪},
    B::Type{<:IntBase}=DEF[:IB])::h_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = begin
    𝑇 = T_amt{𝕡,𝕩}(𝑇)
    return h_(𝐻, 𝑇, B)
end

# Fallback method with hasTPair arguments
(h_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::hasT{𝕢,𝕪},
    B::Type{<:IntBase}=DEF[:IB])::h_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = h_(𝐻, 𝑇.T, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    # Δs°: Particular gas variation of ideal gas partial spec. entropy #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback method
(Δs°(𝐻::nobleGasHeat{𝕡,𝕩},
     𝒾::T_amt{𝕡,𝕩},
     𝒻::T_amt{𝕡,𝕩},
     B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    ds(cp(𝐻, B) * log(𝒻/𝒾))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(Δs°(𝐻::nobleGasHeat{𝕡,𝕩},
      𝒾::T_amt{𝕢,𝕪},
      𝒻::T_amt{𝕣,𝕫},
      B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫}`\n
Returns the particular gas variation in ideal gas partial specific entropy in the specified or
default base for the substance with specific heat modeled by `𝐻`, for process with initial and
final temperatures of `𝒾` and `𝒻`, respectively. Resulting precision, PREC, and exactness, EXAC,
are model-driven, and not promotion-driven.
"""
(Δs°(𝐻::nobleGasHeat{𝕡,𝕩},
     𝒾::T_amt{𝕢,𝕪},
     𝒻::T_amt{𝕣,𝕫},
     B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    𝑖 = T_amt{𝕡,𝕩}(𝑖)
    𝑓 = T_amt{𝕡,𝕩}(𝑓)
    return Δs°(𝐻, 𝑖, 𝑓, B)
end

# Fallback method with hasTPair arguments
(Δs°(𝐻::nobleGasHeat{𝕡,𝕩},
     𝒾::hasT{𝕢,𝕪},
     𝒻::hasT{𝕣,𝕫},
     B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = Δs°(𝐻, 𝒾.T, 𝒻.T, B)

# Alias
ds0 = Δs°


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #      s°: Particular gas specific ideal gas partial entropy       #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback method
(s°(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕡,𝕩},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    s_(Δs°(𝐻, Tref(𝐻), 𝑇, B) + sref(𝐻, B))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(s°(𝐻::nobleGasHeat{𝕡,𝕩},
     𝑇::T_amt{𝕢,𝕪},
     B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪,B}`\n
Returns the particular gas specific ideal gas partial entropy in the specified or default base
for the substance with specific heat modeled by `𝐻`, for states with temperature `𝑇`. Resulting
precision, PREC, and exactness, EXAC, are model-driven, and not promotion-driven.
"""
(s°(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕢,𝕪},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪,B} = begin
    𝑇 = T_amt{𝕡,𝕩}(𝑇)
    return s°(𝐻, 𝑇, B)
end

# Fallback method with hasTPair arguments
(s°(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::hasT{𝕢,𝕪},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪,B} = s°(𝐻, 𝑇.T, B)

# Alias
s0 = s°


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #         Δs: Particular gas variation of specific entropy         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback methods
(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    Ti::T_amt{𝕡,𝕩},
    Tf::T_amt{𝕡,𝕩},
    Pi::P_amt{𝕡,𝕩},
    Pf::P_amt{𝕡,𝕩},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    ds(cp(𝐻, B) * log(Tf/Ti) - R_(𝐻, B) * log(Pf/Pi))
end
(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    Ti::T_amt{𝕡,𝕩},
    Tf::T_amt{𝕡,𝕩},
    vi::v_amt{𝕡,𝕩,𝕓},
    vf::v_amt{𝕡,𝕩,𝕓},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕡,𝕩,𝕓} = begin
    ds(cv(𝐻, B) * log(Tf/Ti) + R_(𝐻, B) * log(vf/vi))
end


# Model-driven PREC and EXAC converting (not promoting)
"""
`(ds(𝐻::nobleGasHeat{𝕡,𝕩},
     Ti::T_amt{𝕢,𝕪},
     Tf::T_amt{𝕣,𝕫},
     Pi::P_amt{𝕟,𝕧},
     Pf::P_amt{𝕠,𝕨},
     B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕟,𝕠,𝕡,𝕢,𝕣,𝕧,𝕨,𝕩,𝕪,𝕫,B}`\n
Returns the particular gas variation in specific entropy in the specified or default base for
the substance with specific heat modeled by `𝐻`, for process with initial and final temperatures
and pressures of `Ti` and `Tf`, and `Pi` and `Pf`, respectively. Resulting precision, PREC, and
exactness, EXAC, are model-driven, and not promotion-driven.
"""
(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    Ti::T_amt{𝕢,𝕪},
    Tf::T_amt{𝕣,𝕫},
    Pi::P_amt{𝕟,𝕧},
    Pf::P_amt{𝕠,𝕨},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕟,𝕠,𝕡,𝕢,𝕣,𝕧,𝕨,𝕩,𝕪,𝕫,B} = begin
    Ti = T_amt{𝕡,𝕩}(Ti)
    Tf = T_amt{𝕡,𝕩}(Tf)
    Pi = P_amt{𝕡,𝕩}(Pi)
    Pf = P_amt{𝕡,𝕩}(Pf)
    return ds(𝐻, Ti, Tf, Pi, Pf, B)
end
(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    Pi::P_amt{𝕟,𝕧},
    Pf::P_amt{𝕠,𝕨},
    Ti::T_amt{𝕢,𝕪},
    Tf::T_amt{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕟,𝕠,𝕡,𝕢,𝕣,𝕧,𝕨,𝕩,𝕪,𝕫,B} = begin
    ds(𝐻, Ti, Tf, Pi, Pf, B)
end

# Fallback versions with <:EoSPair input types
(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑖::TPPair{𝕢,𝕪}, # initial (T, P)
    𝑓::TPPair{𝕣,𝕫}, # final (T, P)
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫,B} = begin
    ds(𝐻, 𝑖.T, 𝑓.T, 𝑖.P, 𝑓.P, B)
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    Ti::T_amt{𝕢,𝕪},
    Tf::T_amt{𝕣,𝕫},
    vi::v_amt{𝕟,𝕧,𝕓},
    vf::v_amt{𝕠,𝕨,𝕓},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕓,𝕟,𝕠,𝕡,𝕢,𝕣,𝕧,𝕨,𝕩,𝕪,𝕫,B}`\n
Returns the particular gas variation in specific entropy in the specified or default base for
the substance with specific heat modeled by `𝐻`, for process with initial and final temperatures
and specific volumes of `Ti` and `Tf`, and `vi` and `vf`, respectively. Resulting precision,
PREC, and exactness, EXAC, are model-driven, and not promotion-driven.
"""
(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    Ti::T_amt{𝕢,𝕪},
    Tf::T_amt{𝕣,𝕫},
    vi::v_amt{𝕟,𝕧,𝕓},
    vf::v_amt{𝕠,𝕨,𝕓},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕓,𝕟,𝕠,𝕡,𝕢,𝕣,𝕧,𝕨,𝕩,𝕪,𝕫,B} = begin
    Ti = T_amt{𝕡,𝕩}(Ti)
    Tf = T_amt{𝕡,𝕩}(Tf)
    vi = v_amt{𝕡,𝕩,𝕓}(vi)
    vf = v_amt{𝕡,𝕩,𝕓}(vf)
    return ds(𝐻, Ti, Tf, vi, vf, B)
end

(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    vi::v_amt{𝕟,𝕧,𝕓},
    vf::v_amt{𝕠,𝕨,𝕓},
    Ti::T_amt{𝕢,𝕪},
    Tf::T_amt{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕓,𝕟,𝕠,𝕡,𝕢,𝕣,𝕧,𝕨,𝕩,𝕪,𝕫,B} = begin
    return ds(𝐻, Ti, Tf, vi, vf, B)    # fallback
end

# Fallback versions with <:ChFPair input types
(ds(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑖::TvPair{𝕟,𝕧,𝕓}, # initial (T, v)
    𝑓::TvPair{𝕠,𝕨,𝕓}, # final (T, v)
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝕡,𝕩,B}) where {𝕟,𝕠,𝕡,𝕧,𝕨,𝕩,𝕓} = begin
    return ds(𝐻, 𝑖.T, 𝑓.T, 𝑖.v, 𝑓.v, B)
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                s: Particular gas specific entropy                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback methods
(s_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕡,𝕩},
    𝑃::P_amt{𝕡,𝕩},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕩} = begin
    s_(ds(𝐻, Tref(𝐻), 𝑇, Pref(𝐻), 𝑃, B) + sref(𝐻, B))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(s_(𝐻::nobleGasHeat{𝕡,𝕩},
     𝑇::T_amt{𝕢,𝕪},
     𝑃::P_amt{𝕣,𝕫},
     B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B})`\n
Returns the particular gas specific entropy in the specified or default base for the substance
with specific heat modeled by `𝐻`, in the specified thermodynamic state (`𝑇`, `𝑃`). Resulting
precision, PREC, and exactness, EXAC, are model-driven, and not promotion-driven.
"""
(s_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕢,𝕪},
    𝑃::P_amt{𝕣,𝕫},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫,B} = begin
    𝑇 = T_amt{𝕡,𝕩}(𝑇)
    𝑃 = P_amt{𝕡,𝕩}(𝑃)
    return s_(𝐻, 𝑇, 𝑃, B)
end

(s_(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑃::P_amt{𝕣,𝕫},
    𝑇::T_amt{𝕢,𝕪},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫,B} = begin
    return s_(𝐻, 𝑇, 𝑃, B)
end

# Fallback method with TPPair arguments
(s_(𝐻::nobleGasHeat{𝕡,𝕩},
    þ::TPPair{𝕢,𝕪},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = s_(𝐻, þ.T, þ.P, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #               Pr: Particular gas relative pressure               #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback methods
(Pr(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕡,𝕩})::Pramt{𝕡,𝕩}) where {𝕡,𝕩} = begin
    Pr(exp(s°(𝐻, 𝑇, 𝕓) / R_(𝐻, 𝕓)))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(Pr(𝐻::nobleGasHeat{𝕡,𝕩},
     𝑇::T_amt{𝕢,𝕪})::Pramt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪}`\n
Returns the particular gas relative pressure for the substance with specific heat modeled by
`𝐻`, in the specified thermodynamic temperature `𝑇`. Resulting precision, PREC, and exactness,
EXAC, are model-driven, and not promotion-driven.
"""
(Pr(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕢,𝕪})::Pramt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = begin
    𝑇 = T_amt{𝕡,𝕩}(𝑇)
    return Pr(𝐻, 𝑇)
end

# Fallback method with hasTPair arguments
(Pr(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::hasT{𝕢,𝕪})::Pramt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = Pr(𝐻, 𝑇.T)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                vr: Particular gas relative volume                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Type-homogeneous fallback methods
(vr(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕡,𝕩})::vramt{𝕡,𝕩}) where {𝕡,𝕩} = begin
    # The be(𝕡(ℯ)) term is a scale factor to render the numerator dimensionless
    vr(𝑇 * be(𝕡(ℯ)) / Pr(𝐻, 𝑇))
end

# Model-driven PREC and EXAC converting (not promoting)
"""
`(vr(𝐻::nobleGasHeat{𝕡,𝕩},
     𝑇::T_amt{𝕢,𝕪})::vramt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪}`\n
Returns the particular gas relative volume for the substance with specific heat modeled by `𝐻`,
in the specified thermodynamic temperature `𝑇`. Resulting precision, PREC, and exactness, EXAC,
are model-driven, and not promotion-driven.
"""
# Type-homogeneous fallback methods
(vr(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::T_amt{𝕢,𝕪})::vramt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = begin
    𝑇 = T_amt{𝕡,𝕩}(𝑇)
    return vr(𝐻, 𝑇)
end

# Fallback method with hasTPair arguments
(vr(𝐻::nobleGasHeat{𝕡,𝕩},
    𝑇::hasT{𝕢,𝕪})::vramt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = vr(𝐻, 𝑇.T)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                         RT: RT products                          #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas constant -- function syntax thanks to
# https://stackoverflow.com/a/65890762/4038337
"""
`(RT(𝐻::nobleGasHeat{𝕡,𝕩}, T::T_amt{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])::RTamt{𝕡,𝕩,B}) where {𝕡,𝕩}`\n
Returns the particular gas (RT) product based on the provided temperature and optional base.
"""
(RT(𝐻::nobleGasHeat{𝕡,𝕩},
    T::T_amt{𝕡,𝕩},
    B::Type{<:IntBase} = DEF[:IB])::RTamt{𝕡,𝕩,B}) where {𝕡,𝕩} = R_(𝐻, B) * T

(RT(𝐻::nobleGasHeat{𝕡,𝕩},
    B::Type{<:IntBase},
    T::T_amt{𝕡,𝕩})::RTamt{𝕡,𝕩,B}) where {𝕡,𝕩} = RT(𝐻, T, B)


    # !center 64 | frame 64 -f'\#⋅\# ' | center 76
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                         Pv: Pv products                          #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Pv(𝐻::nobleGasHeat{𝕡,𝕩}, T::T_amt{𝕡,𝕩}, B::Type{<:IntBase} = DEF[:IB])::Pvamt{𝕡,𝕩,B}) where {𝕡,𝕩}`\n
Returns the particular gas (Pv) product based on the provided temperature and optional base.
"""
(Pv(𝐻::nobleGasHeat{𝕡,𝕩},
    T::T_amt{𝕡,𝕩},
    B::Type{<:IntBase} = DEF[:IB])::Pvamt{𝕡,𝕩,B}) where {𝕡,𝕩} = Pv(RT(𝐻, T, B))

(Pv(𝐻::nobleGasHeat{𝕡,𝕩},
    B::Type{<:IntBase},
    T::T_amt{𝕡,𝕩})::Pvamt{𝕡,𝕩,B}) where {𝕡,𝕩} = Pv(RT(𝐻, T, B))


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                       Ideal Gas functions                        #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# ideal gas Z, by definition
(Z_(𝐻::nobleGasHeat{𝕡,EX})::Z_amt{𝕡,EX}) where {𝕡} = Z_(one(𝕡))
(Z_(𝐻::nobleGasHeat{𝕡,MM})::Z_amt{𝕡,MM}) where {𝕡} = Z_(one(𝕡) ± zero(𝕡))

# calculated gas Z, by definition of Z
"""
`(Z_(𝐻::nobleGasHeat{𝕡,𝕩}, T::T_amt{𝕡,𝕩})::Z_amt{𝕡,𝕩}) where {𝕡,𝕩}`\n
Returns the (ideal gas) generalized compressibility factor from it's \$Pv/RT\$ definition.
"""
(Z_(𝐻::nobleGasHeat{𝕡,𝕩}, T::T_amt{𝕡,𝕩})::Z_amt{𝕡,𝕩}) where {𝕡,𝕩} = Pv(𝐻, T) / RT(𝐻, T)


#----------------------------------------------------------------------------------------------#
#                                        Alias exports                                         #
#----------------------------------------------------------------------------------------------#

export du, dh, ds0, s0


#----------------------------------------------------------------------------------------------#
#                                           Includes                                           #
#----------------------------------------------------------------------------------------------#

include("nobleGas-oper.jl")


