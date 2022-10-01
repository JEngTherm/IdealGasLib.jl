#----------------------------------------------------------------------------------------------#
#                        Noble Gas Specific Heat Model for Ideal Gases                         #
#----------------------------------------------------------------------------------------------#

import Base: show

# Type declaration
struct nobleGasHeat{𝗽,𝘅,𝗯<:IntBase} <: ConstHeat{𝗽,𝘅}
    M::m_amt{𝗽,𝘅,MO}     # The precision- exactness- parametric molar mass
    c::cpamt{𝗽,𝘅,𝗯}     # The precision- exactness- base- parametric cp
    Tref::T_amt{𝗽,𝘅}     # The reference state temperature
    Pref::P_amt{𝗽,𝘅}     # The reference state pressure
    sref::s_amt{𝗽,𝘅,𝗯}   # The reference state specific entropy
    # Inner copy constructor
    nobleGasHeat(x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
        new{𝗽,𝘅,𝗯}(x.M, x.c, x.Tref, x.Pref, x.sref)
    end
    # Inner checking & promoting constructor
    nobleGasHeat(__M::m_amt{𝗽𝗔,𝘅𝗔,MO},
                 __c::cpamt{𝗽𝗕,𝘅𝗕,𝗯},
                 T_r::T_amt{𝗽𝗖,𝘅𝗖}   = 𝗧(promote_type(𝗽𝗔, 𝗽𝗕), promote_type(𝘅𝗔, 𝘅𝗕)),
                 P_r::P_amt{𝗽𝗗,𝘅𝗗}   = 𝗣(promote_type(𝗽𝗔, 𝗽𝗕), promote_type(𝘅𝗔, 𝘅𝗕)),
                 s_r::s_amt{𝗽𝗘,𝘅𝗘,𝗯} = s_amt{promote_type(𝗽𝗔, 𝗽𝗕),promote_type(𝘅𝗔, 𝘅𝗕),𝗯}(
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
`(Tref(x::nobleGasHeat{𝗽,𝘅})::T_amt{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular gas's reference state temperature for the substance with specific heat
modeled by `x`.
"""
(Tref(x::nobleGasHeat{𝗽,𝘅})::T_amt{𝗽,𝘅}) where {𝗽,𝘅} = x.Tref

"""
`(Pref(x::nobleGasHeat{𝗽,𝘅})::P_amt{𝗽,𝘅}) where {𝗽,𝘅}`\n
Returns a particular gas's reference state pressure for the substance with specific heat modeled
by `x`.
"""
(Pref(x::nobleGasHeat{𝗽,𝘅})::P_amt{𝗽,𝘅}) where {𝗽,𝘅} = x.Pref

"""
`(sref(x::nobleGasHeat{𝗽,𝘅,𝗯})::s_amt{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯}`\n
Returns a particular gas's reference state specific entropy for the substance with specific heat
modeled by `x`.
"""
(sref(x::nobleGasHeat{𝗽,𝘅,𝗯})::s_amt{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = x.sref

(sref(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MA})::s_amt{𝗽,𝘅,MA}) where {𝗽,𝘅} = x.sref / x.M
(sref(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MO})::s_amt{𝗽,𝘅,MO}) where {𝗽,𝘅} = x.sref * x.M

(sref(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MO})::s_amt{𝗽,𝘅,MO}) where {𝗽,𝘅} = x.sref
(sref(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MA})::s_amt{𝗽,𝘅,MA}) where {𝗽,𝘅} = x.sref


#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
#                         Basic Ideal Gas Properties from nobleGasHeat                         #
#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                 M: Particular gas molecular mass                 #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas molecular mass
"""
`(m_(x::nobleGasHeat{𝗽,𝘅})::m_amt{𝗽,𝘅,MO}) where {𝗽,𝘅}`\n
Returns the particular gas molecular mass for the substance with specific heat modeled by `x`
without conversions.
"""
(m_(x::nobleGasHeat{𝗽,𝘅})::m_amt{𝗽,𝘅,MO}) where {𝗽,𝘅} = x.M


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                    R: Particular gas constant                    #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas constant -- function syntax thanks to
# https://stackoverflow.com/a/65890762/4038337
"""
`(R_(x::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])::R_amt{𝗽,𝘅,B}) where {𝗽,𝘅}`\n
Returns the particular gas constant for the substance with specific heat modeled by `x` in the
default or specified base.
"""
(R_(x::nobleGasHeat{𝗽,𝘅}, B::Type{MA})::R_amt{𝗽,𝘅,MA}) where {𝗽,𝘅} = R_(𝗽, 𝘅) / x.M
(R_(x::nobleGasHeat{𝗽,𝘅}, B::Type{MO})::R_amt{𝗽,𝘅,MO}) where {𝗽,𝘅} = R_(𝗽, 𝘅)

(R_(x::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])::R_amt{𝗽,𝘅,B}) where {𝗽,𝘅} = R_(x, B)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              cp: Particular gas iso-P specific heat              #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas cp values: no conversion
"""
`cp(x::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])`\n
Returns the particular gas constant-pressure specific heat in the default or specified base for
the substance with specific heat modeled by `x`, making base conversion only when necessary.
"""
(cp(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MA})::cpamt{𝗽,𝘅,MA}) where {𝗽,𝘅} = x.c
(cp(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MO})::cpamt{𝗽,𝘅,MO}) where {𝗽,𝘅} = x.c

# Particular gas cp values: w/ conversion
(cp(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MO})::cpamt{𝗽,𝘅,MO}) where {𝗽,𝘅} = cp(x.c * x.M)
(cp(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MA})::cpamt{𝗽,𝘅,MA}) where {𝗽,𝘅} = cp(x.c / x.M)

# Particular gas cp value: default base fallback
cp(x::nobleGasHeat) = cp(x, DEF[:IB]) # fallback


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #              cv: Particular gas iso-V specific heat              #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Particular gas cv values: no base conversion
"""
`cv(x::nobleGasHeat{𝗽,𝘅}, B::Type{<:IntBase} = DEF[:IB])`\n
Returns the particular gas constant-volume specific heat in the default or specified base for
the substance with specific heat modeled by `x`, making base conversion only when necessary.
"""
(cv(x::nobleGasHeat{𝗽,𝘅,MA}, B::Type{MA})::cvamt{𝗽,𝘅,MA}) where {𝗽,𝘅} = cv(x.c - R(x, MA))
(cv(x::nobleGasHeat{𝗽,𝘅,MO}, B::Type{MO})::cvamt{𝗽,𝘅,MO}) where {𝗽,𝘅} = cv(x.c - R(x, MO))

# Particular gas cv values: w/ base conversion
(cv(x::nobleGasHeat{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::cvamt{𝗽,𝘅,B}) where {𝗽,𝘅} = begin
    cv(cp(x, B) - R_(x, B))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #             ga: Particular gas specific heat ratio               #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(ga(x::nobleGasHeat{𝗽,𝘅,𝗯})::gaamt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas specific heat ratio for the substance with specific heat modeled by
`x`, without conversions.
"""
(ga(x::nobleGasHeat{𝗽,𝘅,𝗯})::gaamt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = ga(cp(x, 𝗯)/cv(x, 𝗯))


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #         k: Particular gas isentropic expansion exponent          #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(k_(x::nobleGasHeat{𝗽,𝘅,𝗯})::k_amt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas isentropic expansion exponent for the substance with specific heat
modeled by `x`, without conversions. For ideal gases, \$k = ga\$.
"""
(k_(x::nobleGasHeat{𝗽,𝘅,𝗯})::k_amt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = k_(ga(x))  # ga fallback


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #     Δu: Particular gas variation of specific internal energy     #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δu(x::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::T_amt{𝗽,𝘅},
     Tf::T_amt{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific internal energy in the specified or default
base for the substance with specific heat modeled by `x`, for process with initial and final
temperatures of `Ti` and `Tf`, respectively.
"""
(Δu(x::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    de(cv(x, B) * (Tf - Ti))
end

# Alias
du = Δu


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #            u: Particular gas specific internal energy            #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(u_(x::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::u_amt{𝗽,𝘅,B})`\n
Returns the particular gas specific internal energy in the specified or default
base for the substance with specific heat modeled by `x`, for states with temperature `theT`.
"""
(u_(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::u_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    u_(Δu(x, Tref(x), theT, B))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #        Δh: Particular gas variation of specific enthalpy         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δh(x::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::T_amt{𝗽,𝘅},
     Tf::T_amt{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific enthalpy in the specified or default base for
the substance with specific heat modeled by `x`, for process with initial and final temperatures
of `Ti` and `Tf`, respectively.
"""
(Δh(x::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::deamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    de(cp(x, B) * (Tf - Ti))
end

# Alias
dh = Δh


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #               h: Particular gas specific enthalpy                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(h_(x::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::h_amt{𝗽,𝘅,B})`\n
Returns the particular gas specific enthalpy in the specified or default base for the substance
with specific heat modeled by `x`, for states with temperature `theT`.
"""
(h_(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::h_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    h_(Δh(x, Tref(x), theT, B) + R_(x, B) * Tref(x))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    # Δs°: Particular gas variation of ideal gas partial spec. entropy #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Δs°(x::nobleGasHeat{𝗽,𝘅,𝗯},
      Ti::T_amt{𝗽,𝘅},
      Tf::T_amt{𝗽,𝘅},
      B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in ideal gas partial specific entropy in the specified or
default base for the substance with specific heat modeled by `x`, for process with initial and
final temperatures of `Ti` and `Tf`, respectively.
"""
(Δs°(x::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::T_amt{𝗽,𝘅},
     Tf::T_amt{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    ds(cp(x, B) * log(Tf/Ti))
end

# Alias
ds0 = Δs°


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #      s°: Particular gas specific ideal gas partial entropy       #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(s°(x::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B})`\n
Returns the particular gas specific ideal gas partial entropy in the specified or default base
for the substance with specific heat modeled by `x`, for states with temperature `theT`.
"""
(s°(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    s_(Δs°(x, Tref(x), theT, B))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #         Δs: Particular gas variation of specific entropy         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(ds(x::nobleGasHeat{𝗽,𝘅,𝗯},
     Ti::T_amt{𝗽,𝘅},
     Tf::T_amt{𝗽,𝘅},
     Pi::P_amt{𝗽,𝘅},
     Pf::P_amt{𝗽,𝘅},
     B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas variation in specific entropy in the specified or default base for
the substance with specific heat modeled by `x`, for process with initial and final temperatures
and pressures of `Ti` and `Tf`, and `Pi` and `Pf`, respectively.
"""
(ds(x::nobleGasHeat{𝗽,𝘅,𝗯},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    Pi::P_amt{𝗽,𝘅},
    Pf::P_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    ds(cp(x, B) * log(Tf/Ti) - R_(x, B) * log(Pf/Pi))
end

(ds(x::nobleGasHeat{𝗽,𝘅,𝗯},
    Pi::P_amt{𝗽,𝘅},
    Pf::P_amt{𝗽,𝘅},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    ds(x, Ti, Tf, Pi, Pf, B)
end

"""
`(ds(x::nobleGasHeat{𝗽,𝘅,𝗯𝗔},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    vi::v_amt{𝗽,𝘅,𝗯𝗕},
    vf::v_amt{𝗽,𝘅,𝗯𝗕},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯𝗔,𝗯𝗕}`\n
Returns the particular gas variation in specific entropy in the specified or default base for
the substance with specific heat modeled by `x`, for process with initial and final temperatures
and specific volumes of `Ti` and `Tf`, and `vi` and `vf`, respectively.
"""
(ds(x::nobleGasHeat{𝗽,𝘅,𝗯𝗔},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    vi::v_amt{𝗽,𝘅,𝗯𝗕},
    vf::v_amt{𝗽,𝘅,𝗯𝗕},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯𝗔,𝗯𝗕} = begin
    ds(cv(x, B) * log(Tf/Ti) + R_(x, B) * log(vf/vi))
end

(ds(x::nobleGasHeat{𝗽,𝘅,𝗯𝗔},
    vi::v_amt{𝗽,𝘅,𝗯𝗕},
    vf::v_amt{𝗽,𝘅,𝗯𝗕},
    Ti::T_amt{𝗽,𝘅},
    Tf::T_amt{𝗽,𝘅},
    B::Type{<:IntBase} = DEF[:IB])::dsamt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯𝗔,𝗯𝗕} = begin
    ds(x, Ti, Tf, vi, vf, B)    # fallsback
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                s: Particular gas specific entropy                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(s_(x::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅},
     theP::P_amt{𝗽,𝘅},
     B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B})`\n
Returns the particular gas specific entropy in the specified or default base for the substance
with specific heat modeled by `x`, in the specified thermodynamic state (`theT`, `theP`).
"""
(s_(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅},
    theP::P_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    s_(ds(x, Tref(x), theT, Pref(x), theP, B))
end

(s_(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theP::P_amt{𝗽,𝘅},
    theT::T_amt{𝗽,𝘅},
    B::Type{<:IntBase}=DEF[:IB])::s_amt{𝗽,𝘅,B}) where {𝗽,𝘅,𝗯} = begin
    s_(x, theT, theP, B)
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #               Pr: Particular gas relative pressure               #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(Pr(x::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅})::Pramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas relative pressure for the substance with specific heat modeled by
`x`, in the specified thermodynamic temperature `theT`.
"""
(Pr(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅})::Pramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = begin
    Pr(exp(s°(x, theT, 𝗯) / R_(x, 𝗯)))
end


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                vr: Particular gas relative volume                #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

"""
`(vr(x::nobleGasHeat{𝗽,𝘅,𝗯},
     theT::T_amt{𝗽,𝘅})::vramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯}`\n
Returns the particular gas relative volume for the substance with specific heat modeled by `x`,
in the specified thermodynamic temperature `theT`.
"""
(vr(x::nobleGasHeat{𝗽,𝘅,𝗯},
    theT::T_amt{𝗽,𝘅})::vramt{𝗽,𝘅}) where {𝗽,𝘅,𝗯} = begin
    # The be(ℯ) term is a scale factor to render the numerator dimensionless
    vr(theT * be(ℯ) / Pr(x, theT))
end



#----------------------------------------------------------------------------------------------#
#                                        Alias exports                                         #
#----------------------------------------------------------------------------------------------#

export du, dh, ds0, s0


