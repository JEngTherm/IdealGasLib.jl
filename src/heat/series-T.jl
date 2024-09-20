#----------------------------------------------------------------------------------------------#
#                                       heat/series-T.jl                                       #
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
#                             Simple Parametric Power Series Type                              #
#----------------------------------------------------------------------------------------------#

struct simplePowerSeries{𝕡,𝕩,𝕀,𝕆} <: Any where {𝕡<:PREC,
                                                𝕩<:EXAC,
                                                𝕀<:WProperty{𝕡,𝕩},
                                                𝕆<:BProperty{𝕡,𝕩,MO}}
    coef::Vector{ø_amt{𝕡,𝕩}}    # The power series term coefficients
    powr::Vector{ø_amt{𝕡,𝕩}}    # The power series term powers
    Linf::𝕀                     # The inferior limit for input variable
    Lsup::𝕀                     # The superior limit for input variable
    unit::𝕆                     # The properly typed UNIT output
    # Inner checking & promoting constructor
    simplePowerSeries(COEF::Vector{ø_amt{𝕡𝔸,𝕩𝔸}},
                      POWR::Vector{ø_amt{𝕡𝔹,𝕩𝔹}},
                      LINF::WProperty{𝕡ℂ,𝕩ℂ},
                      LSUP::WProperty{𝕡𝔻,𝕩𝔻},
                      UNIT::BProperty{𝕡𝔼,𝕩𝔼,MO} = cpamt{promote_type(𝕡𝔸,𝕡𝔹,𝕡ℂ,𝕡𝔻),
                                                        promote_type(𝕩𝔸,𝕩𝔹,𝕩ℂ,𝕩𝔻),
                                                        MO}(one(promote_type(𝕡𝔸,𝕡𝔹,𝕡ℂ,𝕡𝔻)))
                     ) where {𝕡𝔸,𝕩𝔸,𝕡𝔹,𝕩𝔹,𝕡ℂ,𝕩ℂ,𝕡𝔻,𝕩𝔻,𝕡𝔼,𝕩𝔼} = begin
        # Precision and Exactness promotion
        𝕡 = promote_type(𝕡𝔸, 𝕡𝔹, 𝕡ℂ, 𝕡𝔻, 𝕡𝔼)
        𝕩 = promote_type(𝕩𝔸, 𝕩𝔹, 𝕩ℂ, 𝕩𝔻, 𝕩𝔼)
        # Checks
        @assert length(COEF)        == length(POWR)
        @assert length(COEF)        >  0
        @assert typeof(LINF).name   == typeof(LSUP).name
        @assert 𝕡(bare(LINF))       >= zero(𝕡)
        @assert 𝕡(bare(LINF))       <  𝕡(bare(LSUP))
        @assert 𝕡(pod(UNIT))        == one(𝕡)
        # Parameter construction
        𝕀 = typeof(LINF).name.wrapper{𝕡,𝕩}
        𝕆 = typeof(UNIT).name.wrapper{𝕡,𝕩}
        # Returns
        new{𝕡,𝕩,𝕀,𝕆}(ø_amt{𝕡,𝕩}.(COEF),
                     ø_amt{𝕡,𝕩}.(POWR),
                     𝕀(LINF),
                     𝕀(LSUP),
                     𝕆(UNIT))
    end
end


#----------------------------------------------------------------------------------------------#
#         Univariate Specific Heat Model for Ideal Gases as a Temperature Power Series         #
#----------------------------------------------------------------------------------------------#

