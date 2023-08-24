# !center 92 | frame 92
#----------------------------------------------------------------------------------------------#
#                                    heat/nobleGas.test.jl                                     #
#----------------------------------------------------------------------------------------------#

@testset "heat/nobleGas.test.jl - Inner checking constructors                     " begin
    INPP = Dict(
        :name   => "Helium",
        :form   => "He",
        :M      => m_(4.003, MO),
        :c      => cp((5/2)*R_()),
        :Tref   => T_(),
        :Pref   => P_(),
        :sref   => s_(0, MO),
    )
    PREC = (Float16, Float32, Float64, BigFloat)
    EXAC = (1.0, 1.0 ± 0.08)
    begin
        for 𝕡 in PREC
            for 𝕩 in EXAC
                hHe = nobleGasHeat(
                    𝕡(INPP[:M]) * 𝕩,
                    𝕡(INPP[:c]) * 𝕩,
                )
                @test hHe isa nobleGasHeat{𝕡, 𝕩, MO}
            end
        end
    end
end

