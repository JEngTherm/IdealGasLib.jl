#----------------------------------------------------------------------------------------------#
#                                    subs/idealGas-oper.jl                                     #
#----------------------------------------------------------------------------------------------#

# Metaprogramming on the underlying specific heat model
for FUN in (:+,:-)
    fun = String(FUN)
    @eval $FUN(x::idealGas, y::idealGas, args::Any...) = 
        idealGas(x.name*String($fun)*y.name,
                 x.form*String($fun)*y.form,
                 ($FUN)(x.heat, y.heat, args...))
end

for FUN in (:*,:/)
    fun = String(FUN)
    @eval $FUN(x::idealGas, args::Any...) = 
        idealGas(x.name*String($fun)*args[1],
                 x.form*String($fun)*args[1],
                 ($FUN)(x.heat, args...))
end


