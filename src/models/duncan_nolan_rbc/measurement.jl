"""
```
measurement{T<:AbstractFloat}(m::DuncanNolanRBC{T}, TTT::Matrix{T},
                              RRR::Matrix{T}, CCC::Vector{T})
```

Assign measurement equation

```
y_t = ZZ*s_t + DD + u_t
```

where

```
Var(ϵ_t) = QQ
Var(u_t) = EE
Cov(ϵ_t, u_t) = 0
```
"""
function measurement{T<:AbstractFloat}(m::DuncanNolanRBC{T},
                                       TTT::Matrix{T},
                                       RRR::Matrix{T},
                                       CCC::Vector{T})
    endo     = m.endogenous_states
    endo_new = m.endogenous_states_augmented
    exo      = m.exogenous_shocks
    obs      = m.observables

    _n_observables = n_observables(m)
    _n_states = n_states_augmented(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    ZZ = zeros(_n_observables, _n_states)
    DD = zeros(_n_observables)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    ## Output growth
    ZZ[obs[:obs_gdp], endo[:y_t]]  = 1.0
    ZZ[obs[:obs_gdp], endo[:y_t1]] = -1.0
    DD[obs[:obs_gdp]]              = m[:γ_Q]

    ## Consumption growth
    ZZ[obs[:obs_con], endo[:ctot_t]]  = 1.0
    ZZ[obs[:obs_con], endo[:ctot_t1]] = -1.0
    DD[obs[:obs_con]]                 = m[:γ_Q]

    ## Labor income growth
    ZZ[obs[:obs_wagsal], endo[:nw_t]]  = 1.0
    ZZ[obs[:obs_wagsal], endo[:nw_t1]] = -1.0
    DD[obs[:obs_wagsal]]               = m[:γ_Q]

    # Measurement error
    EE[obs[:obs_gdp], endo[:y_t]]      = m[:e_y]^2
    EE[obs[:obs_con], endo[:ctot_t]]   = m[:e_ctot]^2
    EE[obs[:obs_wagsal], endo[:nw_t]]  = m[:e_nw]^2

    # Variance of innovations
    QQ[exo[:z_sh], exo[:z_sh]]   = (m[:σ_z])^2
    QQ[exo[:g_sh], exo[:g_sh]]   = (m[:σ_g])^2
    QQ[exo[:xi_sh],exo[:xi_sh]]  = (m[:σ_xi])^2

    return Measurement(ZZ, DD, QQ, EE)
end
