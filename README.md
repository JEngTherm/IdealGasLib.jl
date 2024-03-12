![IdealGasLib](https://raw.githubusercontent.com/JEngTherm/IdealGasLib.jl/main/docs/src/assets/logo-036.png)

# IdealGasLib.jl

Ideal Gas Library for Engineering Thermodynamics in Julia

# Description

The `IdealGasLib.jl` package is a Ideal Gas Library for Engineering Thermodynamics in Julia
based on the `EngThermBase.jl` infrastructure.

The ideal gas type, `idealGas` is `PREC`-ision, `EXAC`-tness, and heat-capacity model
parameterized, so that in order to instantiate an ideal gas object, a heat capacity model object
must be instantiated first.

`EngThermBase` heat capacity model hooks include:

```julia
julia> using EngThermBase, TypeTree

julia> print.(TypeTree.tt(Heat));
Heat
 ├─ BivarHeat
 ├─ ConstHeat
 ├─ GenerHeat
 └─ UnvarHeat
```

The heat models are thus classified by the amount of variable parameters, from zero (for the
`ConstHeat` model), to one (for the `UnvarHeat` model), to two (for the `BivarHeat` model), to
any (for the `GenerHeat` model).

The available `IdealGasLib` heat capacity model is of an all-temperature, constant specific
heats one suitable for ideal gases, therefore, when the library is loaded, the `Heat` model tree
becomes:

```julia
julia> using IdealGasLib

julia> print.(TypeTree.tt(Heat));
Heat
 ├─ BivarHeat
 ├─ ConstHeat
 │   └─ nobleGasHeat
 ├─ GenerHeat
 └─ UnvarHeat
```

Thus demosntrating that it hooks the `nobleGasHeat` model under `ConstHeat` `Heat` model, which
is a subtype of `MODEL`:

```julia
julia> supertype(Heat)
MODELS{𝗽, 𝘅} where {𝗽, 𝘅}
```

Moreover, the `IdealGasLib` also hooks under `EngThermBase.Substance` models:

```julia
julia> using EngThermBase, TypeTree

julia> print.(TypeTree.tt(Medium));
Medium
 └─ Substance

julia> using IdealGasLib

julia> print.(TypeTree.tt(Medium));
Medium
 └─ Substance
     └─ idealGas
```

## Ideal Noble Gas Heat Model

The ideal noble gas heat model is a constant-specific-heat for all temperatures. There are no
explicit bounds in temperature included in this model's implementation.

The model is instantiated with a mandatory (i) molar molecular mass, of type `m_amt{𝕡,𝕩,MO}
where {𝕡,𝕩}`, a mandatory (ii) molar specific heat at constant pressure, of type `cpamt{𝕡,𝕩,MO}
where {𝕡,𝕩}`, and optional [iii] reference temperature, `T_amt{𝕡,𝕩} where {𝕡,𝕩}`, and [iv]
reference pressure, `P_amt{𝕡,𝕩} where {𝕡,𝕩}`, and [v] molar specific entropy at the reference
state, `s_amt{𝕡,𝕩,MO} where {𝕡,𝕩}`.

The following is a `nobleGasHeat` instantiation for water vapor, based on tabulated values of
`M`, and `cp`, in a mass base, in which the value of `M` is explicitly and manually used prior
to the instantiation, so as to change the specific heat base, from `MA` to `MO`:

```julia
julia> wMass = m_(18.015, MO)
M₆₄: 18.015 kg/kmol

julia> wSpHt = cp(1.8723, MA) * wMass
c̄p₆₄: 33.729 kJ/K/kmol

julia> wHeat = nobleGasHeat(wMass, wSpHt)
noble-cp(T):
   c̄p₆₄: 33.729 kJ/K/kmol
    M₆₄: 18.015 kg/kmol
    T₆₄: 298.15 K
    P₆₄: 101.35 kPa
    s̄₆₄: 0.0000 kJ/K/kmol
```

Owing to `EngThermBase` functionality, the molecular mass of any molecule on Earth can be
accurately calculated from standard values of elemental isotope mass and abundance fractions,
usig the following syntax:

```julia
julia> m_(molParse("H2O"))	# Water, with default element atomic masses
M₃₂: (18.015 ± 0.00033106) kg/kmol

julia> typeof(ans)
m_amt{Float32, MM, MO}

julia> m_(molParse("C8H18"), EngThermBase.atoM_64)	# Octane, with 64-bit element atomic masses
M₆₄: (114.23 ± 0.0065229) kg/kmol

julia> typeof(ans)
m_amt{Float64, MM, MO}

julia> m_amt{Float32,EX}(m_(molParse("C2H5(OH)")))	# Ethanol, converted to 32-bit, EXact base
M₃₂: 46.068 kg/kmol

julia> typeof(ans)
m_amt{Float32, EX, MO}
```

### Using `nobleGasHeat` objects

Except for the $P-T-v$ behavior, which is left for `idealGas` objects to deal vith, most if not
all other thermodynamic behavior, i.e., the "caloric" side, meaning energy and entropy
quantities are (almost) implemented based on the heat model, thus:

```julia
julia> [ @eval ( $FUN(wHeat, T_(1000)) ) for FUN in (:cp, :cv, :ga, :u_, :h_, :Pr, :vr) ]
7-element Vector{AMOUNTS{Float64, EX}}:
 cp₆₄: 1.8723 kJ/K/kg
 cv₆₄: 1.4108 kJ/K/kg
 γ₆₄: 1.3271 –
 u₆₄: 990.15 kJ/kg
 h₆₄: 1451.7 kJ/kg
 Pr₆₄: 135.54 –
 vr₆₄: 20.055 –

julia> [ @eval ( $FUN(wHeat, T_(1000), P_(100)) ) for FUN in (:s_, ) ]
1-element Vector{s_amt{Float64, EX, MA}}:
 s₆₄: 2.2720 kJ/K/kg

julia> T1 = T_(1000)
T₆₄: 1000.0 K

julia> P1 = P_(100)
P₆₄: 100.00 kPa

julia> u_(wHeat, T1) - T1 * s_(wHeat, T1, P1)
a₆₄: -1281.8 kJ/kg

julia> ans / T1
j₆₄: 1.2818 kJ/K/kg

julia> h_(wHeat, T1) - T1 * s_(wHeat, T1, P1)
g₆₄: -820.29 kJ/kg

julia> ans / T1
y₆₄: 0.82029 kJ/K/kg
```

The example above illustrates that, despite Helmholtz, Massieu, Gibbs, and Plank functions not
being directly implemented just yet, their values can be obtained by applying their definitions,
as shown.

## Ideal Gas Model

Once underlying heat models are understood, `idealGas`es are simple and easy to instantiate, as
they only require a name and a chemical formula, beyond the heat model. This makes it possible
to have a substance modeled in different ways, as it is common in engineering thermodynamics.

It is possible to pick a heat model from a tiny library, implemented in this package, as
follows:

```julia
julia> using IdealGasLib

julia> He = idealGas("Helium", "He", HEAT[:He][Float32][MM])
idealGas{Float32, MM, nobleGasHeat{Float32, MM}}("Helium", "He", noble-cp(T):
   c̄p₃₂: (20.786 ± 0.000037500) kJ/K/kmol
    M₃₂: (4.0026 ± 0.0000020000) kg/kmol
    T₃₂: (298.15 ± 0.000000000000056843) K
    P₃₂: (101.35 ± 0.000000000000014211) kPa
    s̄₃₂: (0.0000 ± 0.0000) kJ/K/kmol)
```

Here, the `HEAT` variable is a dictionary of noble gas models parameterized by dictionary keys,
so that `HEAT[:He][Float32][MM]` points to the 32-bit precision, `MM` (i.e., measurements)
exactness noble gas heat model of constant specific heat for Helium. The output explicitly
states the 32-bit precision, as well as each quantity's uncertainties, after the "±" sign.

Once the ideal gas model is instantiated, $P-T-v$ calculations can _also_ be performed, since
all other calculations that can be made with the underlying heat model, can also be performed
with the ideal gas:

```julia
julia> P_(He, T_(300), v_(1, MA))
P₃₂: (623.18 ± 0.0011666) kPa

julia> T_(He, P_(500), v_(1, MA))
T₃₂: (240.70 ± 0.00045059) K

julia> v_(He, T_(300), P_(500))
v₃₂: (1.2464 ± 0.0000023332) m³/kg

julia> v_(He, T_(300), P_(500), MO)
v̄₃₂: (4.9887 ± 0.0000090000) m³/kmol

```

The functions are quite versatile, accepting arguments in a different order, provided that the
`idealGas` model is always the *first* argument, while the (optional) `BASE`, if required, is
always the *last* argument, so the following also work:

```julia
julia> v_(He, P_(500), T_(300), MO)
v̄₃₂: (4.9887 ± 0.0000090000) m³/kmol

julia> TPstate = TPPair(P_(500), T_(300))
TPPair{Float64, EX}(T₆₄: 300.00 K, P₆₄: 500.00 kPa)

julia> v_(He, TPstate, MO)
v̄₃₂: (4.9887 ± 0.0000090000) m³/kmol

```

Moreover, functions implemented for the underlying heat model also have `idealGas` versions:

```julia
julia> [ @eval ( $FUN(He, T_(1000)) ) for FUN in (:cp, :cv, :ga, :u_, :h_, :Pr, :vr) ]
7-element Vector{AMOUNTS{Float32, MM}}:
 cp₃₂: (5.1932 ± 0.0000097216) kJ/K/kg
 cv₃₂: (3.1159 ± 0.0000058330) kJ/K/kg
 γ₃₂: (1.6667 ± 0.00000000000022352) –
 u₃₂: (2186.9 ± 0.0040939) kJ/kg
 h₃₂: (4264.2 ± 0.0079825) kJ/kg
 Pr₃₂: (20.602 ± 0.0000000061399) –
 vr₃₂: (131.94 ± 0.000011573) –

```

## Author

Prof. C. Naaktgeboren, PhD. [Lattes](http://lattes.cnpq.br/8621139258082919).

Federal University of Technology, Paraná
[(site)](http://portal.utfpr.edu.br/english), Guarapuava Campus.

`NaaktgeborenC <dot!> PhD {at!} gmail [dot!] com`

## License

This project is [licensed](https://github.com/JEngTherm/EngThermBase.jl/blob/master/LICENSE)
under the MIT license.

## Citations

How to cite this project:

```bibtex
@Misc{2023-NaaktgeborenC-IdealGasLib,
  author       = {C. Naaktgeboren},
  title        = {{IdealGasLib.jl} -- Ideal Gas Library for Engineering Thermodynamics in Julia},
  howpublished = {Online},
  month        = {September},
  year         = {2023},
  journal      = {GitHub repository},
  publisher    = {GitHub},
  url          = {https://github.com/JEngTherm/IdealGasLib.jl},
  note         = {release 0.2.1 of 24-03-12},
}
```



