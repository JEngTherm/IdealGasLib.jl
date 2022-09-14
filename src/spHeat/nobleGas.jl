#----------------------------------------------------------------------------------------------#
#                        Noble Gas Specific Heat Model for Ideal Gases                         #
#----------------------------------------------------------------------------------------------#

# Type declaration
struct nobleGasHeat{𝗽,𝘅,𝗯<:IntBase} <: ConstHeat{𝗽,𝘅}
    name::String        # Substance name -- that has the (M, c) values
    form::String        # Substance formula as a String
    M::mAmt{𝗽,𝘅,MO}     # The precision- exactness- parametric molar mass
    c::cpAmt{𝗽,𝘅,𝗯}     # The precision- exactness- base- parametric cp
end

# Type exporting
export nobleGasHeat

# Type displaying
deco(x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = Symbol("noble-cp(T)")

Base.show(io::IO, x::nobleGasHeat{𝗽,𝘅,𝗯}) where {𝗽,𝘅,𝗯} = begin
    if DEF[:pprint]
        print(io,
            "$(string(deco(x))) for $(x.form): ",
            "($(x.c)) ($(x.M))"
        )
    else
        Base.show_default(io, x)
    end
end

