#----------------------------------------------------------------------------------------------#
#                                  Ideal Gas Substance Model                                   #
#----------------------------------------------------------------------------------------------#

import Base: show

# Type declaration
struct idealGas{𝕡<:PREC,𝕩<:EXAC,ℍ<:Heat} <: Substance{𝕡,𝕩}
	name::String        # The substance name
    form::String        # The chemical formula
    heat::ℍ             # The heat capacity model
    # Inner copy constructor
    idealGas(x::idealGas{𝕡,𝕩,ℍ}) where {𝕡,𝕩,ℍ} = begin
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

Base.show(io::IO, x::idealGas{𝕡,𝕩,ℍ}) where {𝕡,𝕩,ℍ} = begin
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
for FUN in (:m_,:R_,:cp,:cv,:ga,:k_,:Δu,:u_,:Δh,:h_,:Δs°,:s°,:ds,:s_,:Pr,:vr,:RT,:Pv)
    @eval $FUN(x::idealGas, args::Any...) = ($FUN)(x.heat, args...)
end


#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#
#                                     Ideal Gas functions                                      #
#⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅⋅#


