# !center 92 | frame 92
#----------------------------------------------------------------------------------------------#
#                                      interface.test.jl                                       #
#----------------------------------------------------------------------------------------------#

@testset "interface.test.jl                                                       " begin
    for funcSymb in (:name, :form, :Tref, :Pref, :sref, :Δu, :Δh, :Δs°, :s°)
        @test isdefined(IdealGasLib, funcSymb)  # Whether defined in Module
        @test isdefined(Main, funcSymb)         # Whether exported
    end
end


