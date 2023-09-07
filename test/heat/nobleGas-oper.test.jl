# !center 92 | frame 92
#----------------------------------------------------------------------------------------------#
#                                  heat/nobleGas-oper.test.jl                                  #
#----------------------------------------------------------------------------------------------#

@testset "heat/nobleGas-oper.test.jl - Comparison tests                           " begin
    hHe = nobleGasHeat(m_(4.003, MO), cp((5/2)*R_()))
    cHe = nobleGasHeat(hHe)
    begin
        @test cHe == hHe
    end
end

