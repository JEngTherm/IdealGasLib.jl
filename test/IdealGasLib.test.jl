# !center 92 | frame 92
#----------------------------------------------------------------------------------------------#
#                                     IdealGasLib.test.jl                                      #
#----------------------------------------------------------------------------------------------#

@testset "IdealGasLib.test.jl                                                     " begin
    @test IdealGasLib isa Module
    @test IdealGasLib.EngThermBase isa Module
end


