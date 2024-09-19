#----------------------------------------------------------------------------------------------#
#                                  Ideal Gas Substance Model                                   #
#----------------------------------------------------------------------------------------------#

import Base: show
import EngThermBase: P_, T_, v_

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
(P_(x::idealGas{𝕡,𝕩},
    T::T_amt{𝕡,𝕩},
    v::v_amt{𝕡,𝕩,MO})::P_amt{𝕡,𝕩}) where {𝕡,𝕩} = Pv(x, T, MO) / v

"""
`(P_(x::idealGas{𝕡,𝕩}, T::T_amt{𝕢,𝕪}, v::v_amt{𝕣,𝕫})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫}`\n
Returns the pressure for the ideal gas `x` at specified temperature `T` and specific volume `v`.
Contrary to most `julia` methods, the `x::idealGas{𝕡,𝕩}` model sets the return value precision
and exactness, `{𝕡,𝕩}` instead of performing data type promotions.
"""
(P_(x::idealGas{𝕡,𝕩}, T::T_amt{𝕢,𝕪}, v::v_amt{𝕣,𝕫})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    T = T_amt{𝕡,𝕩}(T)
    v = v_(x, v)
    return P_(x, T, v)      # fallback
end

# Out-of order methods
(P_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕣,𝕫}, T::T_amt{𝕢,𝕪})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    P_(x, T, v)
end

# Other signatures
(P_(x::idealGas{𝕡,𝕩},
    v::v_amt{𝕣,𝕫}, 𝒯::hasT{𝕢,𝕪})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = P_(x, v, 𝒯.T)
(P_(x::idealGas{𝕡,𝕩},
    𝒯::hasT{𝕢,𝕪}, v::v_amt{𝕣,𝕫})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = P_(x, v, 𝒯.T)

(P_(x::idealGas{𝕡,𝕩},
    𝓋::hasv{𝕣,𝕫}, T::T_amt{𝕢,𝕪})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = P_(x, 𝓋.v, T)
(P_(x::idealGas{𝕡,𝕩},
    T::T_amt{𝕢,𝕪}, 𝓋::hasv{𝕣,𝕫})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = P_(x, 𝓋.v, T)

(P_(x::idealGas{𝕡,𝕩},
    𝓋::hasv{𝕣,𝕫}, 𝒯::hasT{𝕢,𝕪})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = P_(x, 𝓋.v, 𝒯.T)
(P_(x::idealGas{𝕡,𝕩},
    𝒯::hasT{𝕢,𝕪}, 𝓋::hasv{𝕣,𝕫})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = P_(x, 𝓋.v, 𝒯.T)

(P_(x::idealGas{𝕡,𝕩},
    þ::TvPair{𝕢,𝕪})::P_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = P_(x, þ.v, þ.T)


    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                      Temperature Functions                       #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Fallback method, with uniform PREC, EXAC:
(T_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕡,𝕩},
    v::v_amt{𝕡,𝕩,MO})::T_amt{𝕡,𝕩}) where {𝕡,𝕩} = P * v / R_(x, MO)

"""
`(T_(x::idealGas{𝕡,𝕩}, P::P_amt{𝕢,𝕪}, v::v_amt{𝕣,𝕫})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫}`\n
Returns the temperature for the ideal gas `x` at specified pressure `P` and specific volume `v`.
Contrary to most `julia` methods, the `x::idealGas{𝕡,𝕩}` model sets the return value precision
and exactness, `{𝕡,𝕩}` instead of performing data type promotions.
"""
(T_(x::idealGas{𝕡,𝕩}, P::P_amt{𝕢,𝕪}, v::v_amt{𝕣,𝕫})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    P = P_amt{𝕡,𝕩}(P)
    v = v_(x, v)
    return T_(x, P, v)      # fallback
end

# Out-of order methods
(T_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕣,𝕫}, P::P_amt{𝕢,𝕪})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    T_(x, P, v)
end

# Other signatures
(T_(x::idealGas{𝕡,𝕩},
    𝒫::hasP{𝕢,𝕪}, v::v_amt{𝕣,𝕫})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = T_(x, 𝒫.P, v)
(T_(x::idealGas{𝕡,𝕩},
    v::v_amt{𝕣,𝕫}, 𝒫::hasP{𝕢,𝕪})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = T_(x, 𝒫.P, v)

(T_(x::idealGas{𝕡,𝕩},
    𝓋::hasv{𝕢,𝕪}, P::P_amt{𝕣,𝕫})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = T_(x, P, 𝓋.v)
(T_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕣,𝕫}, 𝓋::hasv{𝕢,𝕪})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = T_(x, P, 𝓋.v)

(T_(x::idealGas{𝕡,𝕩},
    𝒫::hasP{𝕣,𝕫}, 𝓋::hasv{𝕢,𝕪})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = T_(x, 𝒫.P, 𝓋.v)
(T_(x::idealGas{𝕡,𝕩},
    𝓋::hasv{𝕢,𝕪}, 𝒫::hasP{𝕣,𝕫})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = T_(x, 𝒫.P, 𝓋.v)

(T_(x::idealGas{𝕡,𝕩},
    þ::PvPair{𝕢,𝕪})::T_amt{𝕡,𝕩}) where {𝕡,𝕢,𝕩,𝕪} = T_(x, þ.P, þ.v)



    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
    #                         Volume Functions                         #
    #⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#

# Base standardization methods
# Fallback methods, with uniform PREC, EXAC:
(v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕡,𝕩,MO})::v_amt{𝕡,𝕩,MO}) where {𝕡,𝕩} = v
(v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕡,𝕩,MA})::v_amt{𝕡,𝕩,MO}) where {𝕡,𝕩} = v * m_(x)

# Base-explicit methods
(v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕡,𝕩,BA}, B::Type{MO})::v_amt{𝕡,𝕩,MO})
    where {𝕡,𝕩,BA<:IntBase} = v_(x, v)
(v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕡,𝕩,BA}, B::Type{MA})::v_amt{𝕡,𝕩,MA})
    where {𝕡,𝕩,BA<:IntBase} = v_(x, v) / m_(x)

"""
`(v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕢,𝕪,BA})::v_amt{𝕡,𝕩,MO}) where {𝕡,𝕢,𝕩,𝕪,BA<:IntBase}`\n
Returns the `x::idealGas{𝕡,𝕩}` specific volume as `v_amt{𝕡,𝕩,MO}`, thus adopting the model's
precision and exactness rather than doing promotions.
"""
(v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕢,𝕪,BA})::v_amt{𝕡,𝕩,MO}) where {𝕡,𝕢,𝕩,𝕪,BA<:IntBase} = begin
    if 𝕪 == MM
        if 𝕩 == MM
            valv = Measurement{𝕡}(bare(v))  # Transports uncertainty
        else
            valv = 𝕡(pod(v))                # Ignores uncertainty
        end
    else
        if 𝕩 == MM
            valv = Measurement{𝕡}(bare(v))  # Initializes uncertainty = zero(𝕡)
        else
            valv = 𝕡(pod(v))                # No uncertainty
        end
    end
    return v_(x, v_amt{𝕡,𝕩,BA}(valv))       # fallback
end

# Base-explicit methods
(v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕢,𝕪,BA}, B::Type{MO})::v_amt{𝕡,𝕩,MO})
    where {𝕡,𝕢,𝕩,𝕪,BA<:IntBase} = v_(x, v)
(v_(x::idealGas{𝕡,𝕩}, v::v_amt{𝕢,𝕪,BA}, B::Type{MA})::v_amt{𝕡,𝕩,MA})
    where {𝕡,𝕢,𝕩,𝕪,BA<:IntBase} = v_(x, v) / m_(x)

# Ideal Gas calculation methods
# Fallback method, with uniform PREC, EXAC:
(v_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕡,𝕩},
    T::T_amt{𝕡,𝕩},
    B::Type{MO})::v_amt{𝕡,𝕩,MO}) where {𝕡,𝕩} = RT(x, T, MO) / P

(v_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕡,𝕩},
    T::T_amt{𝕡,𝕩},
    B::Type{MA})::v_amt{𝕡,𝕩,MA}) where {𝕡,𝕩} = RT(x, T, MA) / P

(v_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕡,𝕩},
    T::T_amt{𝕡,𝕩},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕩} = v_(x, P, T, B)    # fallback

# Different precision and/or exactness
(v_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕢,𝕪},
    T::T_amt{𝕣,𝕫},
    B::Type{MO})::v_amt{𝕡,𝕩,MO}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    P = P_amt{𝕡,𝕩}(P)
    T = T_amt{𝕡,𝕩}(T)
    return v_(x, P, T, MO)
end

(v_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕢,𝕪},
    T::T_amt{𝕣,𝕫},
    B::Type{MA})::v_amt{𝕡,𝕩,MA}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    P = P_amt{𝕡,𝕩}(P)
    T = T_amt{𝕡,𝕩}(T)
    return v_(x, P, T, MA)
end

"""
`v_(x::idealGas{𝕡,𝕩}, P::P_amt{𝕢,𝕪}, T::T_amt{𝕣,𝕫}, B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}`\n
Returns the specific volume (at base `B`) for the ideal gas `x` at specified pressure `P` and
temperature `T`.  Contrary to most `julia` methods, the `x::idealGas{𝕡,𝕩}` model sets the return
value precision and exactness, `{𝕡,𝕩}` instead of performing data type promotions. If ommitted,
the base `B` defaults to `DEF[:IB]` (from `EngThermBase`).
"""
(v_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕢,𝕪},
    T::T_amt{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = begin
    P = P_amt{𝕡,𝕩}(P)
    T = T_amt{𝕡,𝕩}(T)
    return v_(x, P, T, B)
end

# Out-of order methods
(v_(x::idealGas{𝕡,𝕩},
    T::T_amt{𝕣,𝕫},
    P::P_amt{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = v_(x, P, T, B)

# Other signatures
(v_(x::idealGas{𝕡,𝕩},
    𝒫::hasP{𝕢,𝕪},
    T::T_amt{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = v_(x, 𝒫.P, T, B)
(v_(x::idealGas{𝕡,𝕩},
    T::T_amt{𝕣,𝕫},
    𝒫::hasP{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = v_(x, 𝒫.P, T, B)

(v_(x::idealGas{𝕡,𝕩},
    P::P_amt{𝕢,𝕪},
    𝒯::hasT{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = v_(x, P, 𝒯.T, B)
(v_(x::idealGas{𝕡,𝕩},
    𝒯::hasT{𝕣,𝕫},
    P::P_amt{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = v_(x, P, 𝒯.T, B)

(v_(x::idealGas{𝕡,𝕩},
    𝒫::hasP{𝕢,𝕪},
    𝒯::hasT{𝕣,𝕫},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = v_(x, 𝒫.P, 𝒯.T, B)
(v_(x::idealGas{𝕡,𝕩},
    𝒯::hasT{𝕣,𝕫},
    𝒫::hasP{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕣,𝕩,𝕪,𝕫} = v_(x, 𝒫.P, 𝒯.T, B)

(v_(x::idealGas{𝕡,𝕩},
    þ::TPPair{𝕢,𝕪},
    B::Type{<:IntBase} = DEF[:IB])::v_amt{𝕡,𝕩,B}) where {𝕡,𝕢,𝕩,𝕪} = v_(x, þ.P, þ.T, B)


#----------------------------------------------------------------------------------------------#
#                                           Includes                                           #
#----------------------------------------------------------------------------------------------#

include("idealGas-oper.jl")


