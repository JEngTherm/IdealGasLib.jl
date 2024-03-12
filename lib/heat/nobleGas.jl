#----------------------------------------------------------------------------------------------#
#                            Noble Gas Heat Capacity Model Library                             #
#----------------------------------------------------------------------------------------------#

import EngThermBase: atoM_64, atoM_32, atoM_16

function atoMass(𝕡::Type{<:PREC})
    (𝕡 in (Float64, BigFloat)) && (return atoM_64)
    (𝕡 == Float32) && (return atoM_32)
    (𝕡 == Float16) && (return atoM_16)
end

# Library signature
HEAT = Dict{Symbol, Dict{DataType, Dict{DataType, nobleGasHeat}}}()

# Populate the library
for GAS in (:He, :Ne, :Ar, :Kr, :Xe, :Rn)
    HEAT[GAS] = Dict{DataType, Dict{DataType, nobleGasHeat}}()
    for PRE in (Float64, Float32, Float16)
        HEAT[GAS][PRE] = Dict{DataType, nobleGasHeat}()
        for EXA in (MM, EX)
            HEAT[GAS][PRE][EXA] = nobleGasHeat(
                m_amt{PRE,EXA}(
                    m_(molParse(String(GAS)), atoMass(PRE))),
                cp((5//2)R_(PRE, EXA)),
            )
        end
    end
end

# export
export HEAT


