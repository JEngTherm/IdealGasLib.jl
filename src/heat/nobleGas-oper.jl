#----------------------------------------------------------------------------------------------#
#                                    heat/nobleGas-oper.jl                                     #
#----------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------#
#                Operations with Noble Gas Specific Heat Model for Ideal Gases                 #
#----------------------------------------------------------------------------------------------#

import Base: +, -, *, /, zero

# Adittion, promoting while keeping the first operand's base
+(𝐴::nobleGasHeat{𝗽𝗔,𝘅𝗔,𝗯𝗔},
  𝐵::nobleGasHeat{𝗽𝗕,𝘅𝗕,𝗯𝗕}) where {𝗽𝗔,𝘅𝗔,𝗽𝗕,𝘅𝗕,𝗯𝗔,𝗯𝗕} = begin
    # Precision and Exactness promotion
    𝗽 = promote_type(𝗽𝗔, 𝗽𝗕)
    𝘅 = promote_type(𝘅𝗔, 𝘅𝗕)
    @assert T_amt{𝗽}(𝐴.Tref) == T_amt{𝗽}(𝐵.Tref)
    @assert P_amt{𝗽}(𝐴.Pref) == P_amt{𝗽}(𝐵.Pref)
    nobleGasHeat(
        m_amt{𝗽}(𝐴.M) + m_amt{𝗽}(𝐵.M),
        cpamt{𝗽}(𝐴.c) + cpamt{𝗽}(cp(𝐵, 𝗯𝗔)),
        T_amt{𝗽}(𝐴.Tref),
        P_amt{𝗽}(𝐴.Pref),
        s_amt{𝗽}(𝐴.sref) + s_amt{𝗽}(sref(𝐵, 𝗯𝗔))
    )
end

# Subtraction, promoting while keeping the first operand's base
-(𝐴::nobleGasHeat{𝗽𝗔,𝘅𝗔,𝗯𝗔},
  𝐵::nobleGasHeat{𝗽𝗕,𝘅𝗕,𝗯𝗕}) where {𝗽𝗔,𝘅𝗔,𝗽𝗕,𝘅𝗕,𝗯𝗔,𝗯𝗕} = begin
    # Precision and Exactness promotion
    𝗽 = promote_type(𝗽𝗔, 𝗽𝗕)
    𝘅 = promote_type(𝘅𝗔, 𝘅𝗕)
    @assert T_amt{𝗽}(𝐴.Tref) == T_amt{𝗽}(𝐵.Tref)
    @assert P_amt{𝗽}(𝐴.Pref) == P_amt{𝗽}(𝐵.Pref)
    # Inner constructor checks for {M, c, Tref, Pref} out of physical bounds
    nobleGasHeat(
        m_amt{𝗽}(𝐴.M) - m_amt{𝗽}(𝐵.M),
        cpamt{𝗽}(𝐴.c) - cpamt{𝗽}(cp(𝐵, 𝗯𝗔)),
        T_amt{𝗽}(𝐴.Tref),
        P_amt{𝗽}(𝐴.Pref),
        s_amt{𝗽}(𝐴.sref) - s_amt{𝗽}(sref(𝐵, 𝗯𝗔))
    )
end

# Scalar multiplication, promoting while keeping the first operand's base
*(𝐴::nobleGasHeat{𝗽𝗔,𝘅,𝗯}, N::EngThermBase.plnF{𝗽𝗡}) where {𝗽𝗔,𝘅,𝗯,𝗽𝗡<:PREC} = begin
    # Precision and Exactness promotion
    𝗽 = promote_type(𝗽𝗔, 𝗽𝗡)
    nobleGasHeat(
        m_amt{𝗽}(𝐴.M) * N,
        cpamt{𝗽}(𝐴.c) * N,
        T_amt{𝗽}(𝐴.Tref),
        P_amt{𝗽}(𝐴.Pref),
        s_amt{𝗽}(𝐴.sref),
    )
end
# Fallback version
*(N::EngThermBase.plnF{𝗽𝗡}, 𝐴::nobleGasHeat{𝗽𝗔,𝘅,𝗯}) where {𝗽𝗔,𝘅,𝗯,𝗽𝗡<:PREC} = 𝐴 * N

# Division by scalar, promoting
/(𝐴::nobleGasHeat{𝗽𝗔,𝘅,𝗯}, N::EngThermBase.plnF{𝗽𝗡}) where {𝗽𝗔,𝘅,𝗯,𝗽𝗡<:PREC} = begin
    # Precision and Exactness promotion
    𝗽 = promote_type(𝗽𝗔, 𝗽𝗡)
    nobleGasHeat(
        m_amt{𝗽}(𝐴.M) / N,
        cpamt{𝗽}(𝐴.c) / N,
        T_amt{𝗽}(𝐴.Tref),
        P_amt{𝗽}(𝐴.Pref),
        s_amt{𝗽}(𝐴.sref),
    )
end

# Base.zero
(zero(𝐴::nobleGasHeat{𝗽,𝘅,𝗯})::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
    nobleGasHeat(
        zero(𝐴.M),
        zero(𝐴.c),
        𝐴.Tref,
        𝐴.Pref,
        zero(𝐴.sref)
    )
end


#······························································································#
#                                          Inquiring                                           #
#······························································································#

import EngThermBase: precof, exacof, baseof

"""
`precof(::Type{𝗧} | x::𝗧) where 𝗧<:Heat{𝗽} where 𝗽 = 𝗽`\n
Returns the precision of the `Heat` subtype or instance as a `DataType`.
"""
precof(::Type{𝗧}) where 𝗧<:Heat{𝗽} where 𝗽 = 𝗽
precof(x::𝗧) where 𝗧<:Heat{𝗽} where 𝗽 = 𝗽

"""
`exacof(::Type{𝗧} | x::𝗧) where 𝗧<:Heat{𝗽} where 𝗽 = 𝗽`\n
Returns the exactness of the `Heat` subtype or instance as a `DataType`.
"""
exacof(::Type{𝗧}) where 𝗧<:Heat{𝗽,𝘅} where {𝗽,𝘅} = 𝘅
exacof(x::𝗧) where 𝗧<:Heat{𝗽,𝘅} where {𝗽,𝘅} = 𝘅

"""
`baseof(::Type{𝗧} | x::𝗧) where 𝗧<:BasedAmt{𝗽,𝘅,𝗯} where {𝗽,𝘅,𝗯} = 𝗯`\n
Returns the thermodynamic base of the `Heat` subtype or instance as a `DataType`.
"""
baseof(::Type{𝗧}) where 𝗧<:Heat{𝗽,𝘅,𝗯} where {𝗽,𝘅,𝗯} = 𝗯
baseof(x::𝗧) where 𝗧<:Heat{𝗽,𝘅,𝗯} where {𝗽,𝘅,𝗯} = 𝗯


#······························································································#
#                                            Mixing                                            #
#······························································································#

# Mixing, promoting while keeping the first heat operand's base
(mx(ys::NTuple{𝗡,EngThermBase.plnF{𝗽𝗬}},
    hs::Tuple{
        nobleGasHeat{𝗽,𝘅,𝗯},
        Vararg{nobleGasHeat,𝗡}
    })::nobleGasHeat{promote_type(𝗽𝗬, 𝗽),𝘅,𝗯}) where {𝗡,𝗽𝗬,𝗽,𝘅,𝗯} = begin
    Σy = sum(ys)
    yr = one(𝗽𝗬) - Σy
    return (hcat(ys..., yr) * vcat(hs...))[1]
end

export mx


