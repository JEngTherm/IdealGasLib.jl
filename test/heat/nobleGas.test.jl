# !center 92 | frame 92
#----------------------------------------------------------------------------------------------#
#                                    heat/nobleGas.test.jl                                     #
#----------------------------------------------------------------------------------------------#

@testset "heat/nobleGas.test.jl - Inner copy constructor                          " begin
    hHe = nobleGasHeat(m_(4.003, MO), cp((5/2)*R_()))
    cHe = nobleGasHeat(hHe)
    begin
        @test typeof(cHe) == typeof(hHe)
        @test cHe.M == hHe.M
        @test cHe.c == hHe.c
        @test cHe.Tref == hHe.Tref
        @test cHe.Pref == hHe.Pref
        @test cHe.sref == hHe.sref
    end
end

@testset "heat/nobleGas.test.jl - Inner checking constructor                      " begin
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
            for 𝕩 in (EX, MM)
                f = 𝕩 == EX ? one(𝕡) : one(𝕡) ± one(𝕡) / 20
                hHe = nobleGasHeat(
                    m_amt{𝕡,𝕩,MO}(𝕡(pod(INPP[:M])) * f),
                    cpamt{𝕡,𝕩,MO}(𝕡(pod(INPP[:c])) * f),
                    T_(𝕡) * f,
                    P_(𝕡) * f,
                    s_(zero(𝕡) * f, MO),
                )
                @test hHe isa nobleGasHeat{𝕡,𝕩,MO}          # expected type {params}
                @test deco(hHe) == Symbol("noble-cp(T)")    # deco independence of {params}
                @test Tref(hHe) == T_(𝕡) * f
                @test Pref(hHe) == P_(𝕡) * f
                @test sref(hHe) == s_(zero(𝕡) * f, MO)
            end
        end
    end
end

