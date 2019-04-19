"""
```
DuncanNolanRBC{T} <: AbstractModel{T}
```

The `DuncanNolanRBC` type defines the structure of a simple RBC model
used in Duncan Nolan 2019.

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model
  parameters.

* `steady_state::Vector{AbstractParameter}`: Model steady-state values, computed
  as a function of elements of `parameters`.

* `keys::OrderedDict{Symbol,Int}`: Maps human-readable names for all model
  parameters and steady-states to their indices in `parameters` and
  `steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readable names to row and
column indices in the matrix representations of of the measurement equation and
equilibrium conditions.

* `endogenous_states::OrderedDict{Symbol,Int}`: Maps each state to a column in
  the measurement and equilibrium condition matrices.

* `exogenous_shocks::OrderedDict{Symbol,Int}`: Maps each shock to a column in
  the measurement and equilibrium condition matrices.

* `expected_shocks::OrderedDict{Symbol,Int}`: Maps each expected shock to a
  column in the measurement and equilibrium condition matrices.

* `equilibrium_conditions::OrderedDict{Symbol,Int}`: Maps each equlibrium
  condition to a row in the model's equilibrium condition matrices.

* `endogenous_states_augmented::OrderedDict{Symbol,Int}`: Maps lagged states to
  their columns in the measurement and equilibrium condition equations. These
  are added after `gensys` solves the model.

* `observables::OrderedDict{Symbol,Int}`: Maps each observable to a row in the
  model's measurement equation matrices.

* `pseudo_observables::OrderedDict{Symbol,Int}`: Maps each pseudo-observable to
  a row in the model's pseudo-measurement equation matrices.

#### Model Specifications and Settings

* `spec::String`: The model specification identifier, \"duncan_nolan_rbc\", cached
  here for filepath computation.

* `subspec::String`: The model subspecification number, indicating that some
  parameters from the original model spec (\"ss0\") are initialized
  differently. Cached here for filepath computation.

* `settings::Dict{Symbol,Setting}`: Settings/flags that affect computation
  without changing the economic or mathematical setup of the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. Can be is seeded to ensure
  reproducibility in algorithms that involve randomness (such as
  Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If `true`,
  settings from `m.test_settings` are used in place of those in `m.settings`.

* `observable_mappings::OrderedDict{Symbol,Observable}`: A dictionary that
  stores data sources, series mnemonics, and transformations to/from model units.
  DSGE.jl will fetch data from the Federal Reserve Bank of
  St. Louis's FRED database; all other data must be downloaded by the
  user. See `load_data` and `Observable` for further details.

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
type DuncanNolanRBC{T} <: AbstractModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                       # model steady-state values
    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
                                                           # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::OrderedDict{Symbol,Int}               #
    equilibrium_conditions::OrderedDict{Symbol,Int}        #
    endogenous_states_augmented::OrderedDict{Symbol,Int}   #
    observables::OrderedDict{Symbol,Int}                   #
    pseudo_observables::OrderedDict{Symbol,Int}            #

    spec::String                                           # Model specification number (eg "m990")
    subspec::String                                        # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                         # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                    # Settings/flags for testing mode
    rng::MersenneTwister                                   # Random number generator
    testing::Bool                                          # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
    pseudo_observable_mappings::OrderedDict{Symbol, PseudoObservable}
end

description(m::DuncanNolanRBC) = "Julia implementation of RBC model DuncanNolanRBC, $(m.subspec)"

"""
`init_model_indices!(m::DuncanNolanRBC)`

Arguments:
`m:: DuncanNolanRBC`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::DuncanNolanRBC)
    # Endogenous states
    endogenous_states = collect([
        :y_t,
        :c_t,
        :nw_t,   # Start with observables in the same order
        :y_t1,
        :c_t1,
        :nw_t1,
        :n_t,
        :w_t,
        :r_t,
        :k_t,
        :i_t,
        :g_t,
        :z_t,
        :Ec_t1,
        :Er_t1])

    # Exogenous shocks
    exogenous_shocks = collect([
        :z_sh, :g_sh])

    # Expectations shocks
    expected_shocks = collect([
        :Ec_sh, :Er_sh])

    # Equilibrium conditions
    equilibrium_conditions = collect([
        :eq_euler,
        :eq_labsup,
        :eq_prod,
        :eq_ad,
        :eq_cap,
        :eq_fpcap,
        :eq_fplab,
        :eq_nw,
        :eq_nw_t1,
        :eq_y_t1,
        :eq_c_t1,
        :eq_g,
        :eq_z,
        :eq_Ec,
        :eq_Er])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(expected_shocks);             m.expected_shocks[k]             = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
    for (i,k) in enumerate(pseudo_observables);          m.pseudo_observables[k]          = i end
end


function DuncanNolanRBC(subspec::String="ss0";
                       custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                       testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = DuncanNolanRBC{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(), OrderedDict{Symbol,Int}(),

            # model indices
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Set settings
    settings_duncan_nolan_rbc!(m)
    default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Set observable and pseudo-observable transformations
    init_observable_mappings!(m)
    init_pseudo_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)

    init_model_indices!(m)
    init_subspec!(m)
    steadystate!(m)

    return m
end

"""
```
init_parameters!(m::DuncanNolanRBC)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::DuncanNolanRBC)
    # Initialize parameters

    m <= parameter(:alpha, 0.25,fixed=true,# (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), GammaAlt(0.25, 0.05), fixed=false,
                   description="alpha: The capital share of output.",
                   tex_label="\\alpha")

    m <= parameter(:gamma, 1.9937,fixed=true,# (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), GammaAlt(2., 0.5), fixed=false,
                   description="gamma: The inverse of the intemporal elasticity of substitution.",
                   tex_label="\\gamma")

    m <= parameter(:psi, 2.0,fixed=true,# (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), GammaAlt(2., 0.5), fixed=false,
                   description="psi: The inverse of the Frisch elasticity of labor supply.",
                   tex_label="\\psi")

    m <= parameter(:delta, 0.025,fixed=true,# (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), Uniform(0.01, 0.04), fixed=false,
                   description="delta: Depreciation rate.",
                   tex_label="\\delta")

    m <= parameter(:beta, 0.995,fixed=true,# (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), Uniform(.99, 0.997), fixed=false,
                   description="beta: Discount factor.",
                   tex_label="\\beta")

    m <= parameter(:CoY, 0.64,fixed=true,# (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), Normal(.64, 0.02), fixed=false,
                   description="C over Y: Consumption share of output.",
                   tex_label="CoY")

    m <= parameter(:IoY, 0.17,fixed=true,# (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), Normal(.17, 0.02), fixed=false,
                   description="I over Y: Investment share of output.",
                   tex_label="IoY")

    m <= parameter(:γ_Q, 1.5, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), Normal(0.40, 0.20), fixed=false,

                   description="γ_Q: Steady state growth rate of technology.",
                   tex_label="\\gamma_Q")

    m <= parameter(:ρ_g, 0.3777, (1e-20, 1-1e-7), (1e-20, 1-1e-7), DSGE.SquareRoot(), Uniform(0,1), fixed=false,
                   description="ρ_g: AR(1) coefficient on g_t government spending as a fraction of output.",
                   tex_label="\\rho_g")

    m <= parameter(:ρ_z, 0.9579, (1e-20, 1-1e-7), (1e-20, 1-1e-7), DSGE.SquareRoot(), Uniform(0,1), fixed=false,
                   description="ρ_z: AR(1) coefficient on shocks to the technology growth rate.",
                   tex_label="\\rho_z")

    m <= parameter(:σ_g, 1.4594, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), DSGE.RootInverseGamma(4, 1.), fixed=false,
                   description="σ_g: Standard deviation of shocks to the government spending process.",
                   tex_label="\\sigma_g")

    m <= parameter(:σ_z, 0.9247, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), DSGE.RootInverseGamma(4, 0.5), fixed=false,
                   description="σ_z: Standard deviation of shocks to the technology growth rate process.",
                   tex_label="\\sigma_z")

    m <= parameter(:e_y, 0.20*0.579923, fixed=true,
                   description="e_y: Measurement error on GDP growth.",
                   tex_label="e_y")

    m <= parameter(:e_c, 0.20*0.579923, fixed=true,
                   description="e_c: Measurement error on consumption growth.",
                   tex_label="e_c")

    m <= parameter(:e_nw, 0.20*0.579923, fixed=true,
                  description="e_{nw}: Measurement error on labor income growth.",
                  tex_label="e_{nw}")
end

"""
```
steadystate!(m::DuncanNolanRBC)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::DuncanNolanRBC)
    return m
end

function settings_duncan_nolan_rbc!(m::DuncanNolanRBC)
    default_settings!(m)

    # Data
    m <= Setting(:data_id, 0, "Dataset identifier")
    m <= Setting(:cond_full_names, [:obs_gdp, :obs_nominalrate],
        "Observables used in conditional forecasts")
    m <= Setting(:cond_semi_names, [:obs_nominalrate],
        "Observables used in semiconditional forecasts")

    # Metropolis-Hastings
    m <= Setting(:mh_cc, 0.27,
                 "Jump size for Metropolis-Hastings (after initialization)")

    # Estimation
    m <= Setting(:reoptimize, true)
    m <= Setting(:recalculate_hessian, true)

    # Forecast
    m <= Setting(:use_population_forecast, false,
                 "Whether to use population forecasts as data")
    m <= Setting(:forecast_zlb_value, 0.13,
        "Value of the zero lower bound in forecast periods, if we choose to enforce it")
end

function shock_groupings(m::DuncanNolanRBC)
    gov = ShockGroup("g", [:g_sh], RGB(0.70, 0.13, 0.13)) # firebrick
    tfp = ShockGroup("z", [:z_sh], RGB(1.0, 0.55, 0.0)) # darkorange
    det = ShockGroup("dt", [:dettrend], :gray40)

    return [gov, tfp, det]
end
