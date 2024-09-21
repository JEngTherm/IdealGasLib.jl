#----------------------------------------------------------------------------------------------#
#                                       heat/series-T.jl                                       #
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
#                             Simple Parametric Power Series Type                              #
#----------------------------------------------------------------------------------------------#

import Base: show #, cp # Since :cp is further extended here
import EngThermBase: deco#, m_, R_, cv, ga, k_, u_, h_, ds, s_, Pr, vr, Pv, RT, Z_

# Type declaration
struct simplePowerSeries{𝕡,𝕩,𝕀,𝕆} <: Any where {𝕡<:PREC,
                                                𝕩<:EXAC,
                                                𝕀<:WProperty{𝕡,𝕩},
                                                𝕆<:BProperty{𝕡,𝕩,MO}}
    coef::Vector{ø_amt{𝕡,𝕩}}    # The power series term coefficients
    powr::Vector{ø_amt{𝕡,𝕩}}    # The power series term powers
    Linf::𝕀                     # The inferior limit for input variable
    Lsup::𝕀                     # The superior limit for input variable
    type::Type{𝕆}               # The output type
    # Inner checking & promoting constructor
    simplePowerSeries(COEF::Vector{ø_amt{𝕡𝔸,𝕩𝔸}},
                      POWR::Vector{ø_amt{𝕡𝔹,𝕩𝔹}},
                      LINF::WProperty{𝕡ℂ,𝕩ℂ},
                      LSUP::WProperty{𝕡𝔻,𝕩𝔻},
                      TYPE::Type{<:BProperty} = Type{cpamt}
                     ) where {𝕡𝔸,𝕩𝔸,𝕡𝔹,𝕩𝔹,𝕡ℂ,𝕩ℂ,𝕡𝔻,𝕩𝔻} = begin
        # Precision and Exactness promotion
        𝕡 = promote_type(𝕡𝔸, 𝕡𝔹, 𝕡ℂ, 𝕡𝔻)
        𝕩 = promote_type(𝕩𝔸, 𝕩𝔹, 𝕩ℂ, 𝕩𝔻)
        # Checks
        @assert length(COEF)        == length(POWR)
        @assert length(COEF)        >  0
        @assert typeof(LINF).name   == typeof(LSUP).name
        @assert 𝕡(bare(LINF))       >= zero(𝕡)
        @assert 𝕡(bare(LINF))       <  𝕡(bare(LSUP))
        # Parameter construction
        𝕀 = typeof(LINF).name.wrapper{𝕡,𝕩}
        while typeof(TYPE) == UnionAll
            TYPE = TYPE.body
        end
        𝕆 = TYPE.name.wrapper{𝕡,𝕩,MO}
        # Returns
        new{𝕡,𝕩,𝕀,𝕆}(ø_amt{𝕡,𝕩}.(COEF),
                     ø_amt{𝕡,𝕩}.(POWR),
                     𝕀(LINF),
                     𝕀(LSUP),
                     𝕆)
    end
end

# Type exporting
export simplePowerSeries

# Type displaying
deco(p::simplePowerSeries{𝕡,𝕩,𝕀,𝕆}) where {𝕡,𝕩,𝕀,𝕆} = begin
    Symbol("$(deco(one(𝕆)))$(EngThermBase.pDeco(𝕡))($(deco(one(𝕀))))")
end

Base.show(io::IO, p::simplePowerSeries{𝕡,𝕩,𝕀,𝕆}) where {𝕡,𝕩,𝕀,𝕆} = begin
    if DEF[:pprint]
        print(io, "$(length(p.coef))-term $(string(deco(p)))")
    else
        print(io, "$(string(deco(p)))")
    end
end

# Functor usage
(p::simplePowerSeries{𝕡,𝕩,𝕀,𝕆})(x::𝕀) where {𝕡,𝕩,𝕀,𝕆} = begin
    @assert p.Linf <= x <= p.Lsup
    𝑥 = bare(x)
    return 𝕆(bare(sum(p.coef .* [𝑥^𝑝 for 𝑝 in p.powr])))
end

(p::simplePowerSeries{𝕡,𝕩,𝕀,𝕆})(x::Vector{𝕀}) where {𝕡,𝕩,𝕀,𝕆} = begin
    return [ p(𝑥) for 𝑥 in x ]
end


#----------------------------------------------------------------------------------------------#
#         Univariate Specific Heat Model for Ideal Gases as a Temperature Power Series         #
#----------------------------------------------------------------------------------------------#

