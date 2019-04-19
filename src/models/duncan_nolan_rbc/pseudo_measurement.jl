"""
```
pseudo_measurement{T<:AbstractFloat}(m::DuncanNolanRBC{T},
    TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T})
```

Assign pseudo-measurement equation (a linear combination of states):

```
x_t = ZZ_pseudo*s_t + DD_pseudo
```
"""
function pseudo_measurement{T<:AbstractFloat}(m::DuncanNolanRBC{T},
                                              TTT::Matrix{T},
                                              RRR::Matrix{T},
                                              CCC::Vector{T})
    endo   = m.endogenous_states
    pseudo = m.pseudo_observables

    _n_states = n_states_augmented(m)
    _n_pseudo = n_pseudo_observables(m)

    # Initialize pseudo ZZ and DD matrices
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)

    ##########################################################
    ## PSEUDO-OBSERVABLE EQUATIONS
    ##########################################################

    ## Output
    ZZ_pseudo[pseudo[:y_t],endo[:y_t]] = 1.

    ## Consumption
    ZZ_pseudo[pseudo[:ctot_t],endo[:ctot_t]] = 1.

    ## Wages and Salaries
    ZZ_pseudo[pseudo[:nw_t],endo[:nw_t]] = 1.

    return PseudoMeasurement(ZZ_pseudo, DD_pseudo)
end
