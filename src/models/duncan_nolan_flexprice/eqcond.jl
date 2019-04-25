"""
```
eqcond(m::Schorf)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
Using the mappings of states/equations to integers defined in schorf.jl, coefficients are
specified in their proper positions.

### Outputs

* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::DuncanNolanFlexPrice)
    endo = m.endogenous_states
    exo  = m.exogenous_shocks
    ex   = m.expected_shocks
    eq   = m.equilibrium_conditions

    Γ0 = zeros(n_states(m), n_states(m))
    Γ1 = zeros(n_states(m), n_states(m)) 
    C  = zeros(n_states(m))
    Ψ  = zeros(n_states(m), n_shocks_exogenous(m))
    Π  = zeros(n_states(m), n_shocks_expectational(m))

    ### ENDOGENOUS STATES ###

    ### 1. Consumption Euler Equation

    Γ0[eq[:eq_euler], endo[:c_t]]   = -1
    Γ0[eq[:eq_euler], endo[:Ec_t1]] =  1
    Γ0[eq[:eq_euler], endo[:Er_t1]] = -1/m[:gamma]

    ### 2. Labour supply

    Γ0[eq[:eq_labsup], endo[:n_t]] = -m[:psi]
    Γ0[eq[:eq_labsup], endo[:w_t]] =  1
    Γ0[eq[:eq_labsup], endo[:c_t]] = -m[:gamma]

    ### 3. Production

    Γ0[eq[:eq_prod], endo[:y_t]] =  1
    Γ0[eq[:eq_prod], endo[:z_t]] = -1
    Γ0[eq[:eq_prod], endo[:n_t]] = -(1-m[:alpha])
    Γ1[eq[:eq_prod], endo[:k_t]] =  m[:alpha]

    ### 4. Consumption aggregator

    Γ0[eq[:eq_con], endo[:ctot_t]] =  -1
    Γ0[eq[:eq_con], endo[:c_t]]    =  m[:CoCtot]
    Γ0[eq[:eq_con], endo[:ce_t]]   = (1 - m[:CoCtot])

    ### 5. Aggregate demand

    Γ0[eq[:eq_ad], endo[:y_t]]    = -1
    Γ0[eq[:eq_ad], endo[:ctot_t]] =  m[:CoY]
    Γ0[eq[:eq_ad], endo[:i_t]]    =  m[:IoY]
    Γ0[eq[:eq_ad], endo[:g_t]]    = (1 - m[:CoY] - m[:IoY])

    ### 6. Capital Accumulation

    Γ0[eq[:eq_cap], endo[:k_t]]    = 1
    Γ1[eq[:eq_cap], endo[:k_t]]    = (1 - m[:delta])
    Γ0[eq[:eq_cap], endo[:i_t]]    = -m[:delta]

    ### 7. Entrepreneurs: Consumption savings

    Γ0[eq[:eq_entcon], endo[:ce_t]]      =  1
    Γ1[eq[:eq_entcon], endo[:omegae_t]]  =  1
    Γ0[eq[:eq_entcon], endo[:re_t]]      = -1

    ### 8. Entrepreneurs: Leverage

    Γ0[eq[:eq_entlev], endo[:lev_t]]     = -1
    Γ0[eq[:eq_entlev], endo[:y_t]]       =  1
    Γ1[eq[:eq_entlev], endo[:omegae_t]]  =  1
    Γ0[eq[:eq_entlev], endo[:r_t]]       = -1

    ### 9. Entrepreneurs: Wedge

    Γ0[eq[:eq_entweg], endo[:lev_t]] = -1
    Γ0[eq[:eq_entweg], endo[:tau_t]] =  m[:L]
    Γ0[eq[:eq_entweg], endo[:xi_t]]  = -(1+m[:erp])

    ### 10. Entrepreneurs: Wealth evolution

    Γ0[eq[:eq_entwel], endo[:omegae_t]] =  1
    Γ0[eq[:eq_entwel], endo[:re_t]]     = -1
    Γ1[eq[:eq_entwel], endo[:omegae_t]] =  1

    ### 11. Factor Prices: Capital

    Γ0[eq[:eq_fpcap], endo[:r_t]]   = -1/(1-m[:beta]*(1-m[:delta]))
    Γ0[eq[:eq_fpcap], endo[:y_t]]   =  1
    Γ1[eq[:eq_fpcap], endo[:k_t]]   =  1
    Γ0[eq[:eq_fpcap], endo[:tau_t]] = -1

    ### 12. Factor Prices: ERP

    Γ0[eq[:eq_fperp], endo[:re_t]]   = -1
    Γ0[eq[:eq_fperp], endo[:r_t]]    =  1
    Γ0[eq[:eq_fperp], endo[:lev_t]]  =  m[:erp]
    Γ0[eq[:eq_fperp], endo[:tau_t]]  =  m[:L]


    ### 13. Factor Prices: Labor

    Γ0[eq[:eq_fplab], endo[:w_t]]   = -1
    Γ0[eq[:eq_fplab], endo[:y_t]]   =  1
    Γ0[eq[:eq_fplab], endo[:n_t]]   = -1
    Γ0[eq[:eq_fplab], endo[:tau_t]] = -1

    ### 14. Wages and salaries

    Γ0[eq[:eq_nw], endo[:nw_t]] = -1
    Γ0[eq[:eq_nw], endo[:w_t]]  =  1
    Γ0[eq[:eq_nw], endo[:n_t]]  =  1

    ### 15. Wages and salaries lag

    Γ0[eq[:eq_nw_t1], endo[:nw_t1]] = 1
    Γ1[eq[:eq_nw_t1], endo[:nw_t]]  = 1

    ### 16. Output lag

    Γ0[eq[:eq_y_t1], endo[:y_t1]] = 1
    Γ1[eq[:eq_y_t1], endo[:y_t]]  = 1

    ### 17. Con lag

    Γ0[eq[:eq_ctot_t1], endo[:ctot_t1]] = 1
    Γ1[eq[:eq_ctot_t1], endo[:ctot_t]]  = 1

    ### 18. Government spending

    Γ0[eq[:eq_g], endo[:g_t]] = 1
    Γ1[eq[:eq_g], endo[:g_t]] = m[:ρ_g]
    Ψ[eq[:eq_g],  exo[:g_sh]] = 1

    ### 19. Technology

    Γ0[eq[:eq_z], endo[:z_t]] = 1
    Γ1[eq[:eq_z], endo[:z_t]] = m[:ρ_z]
    Ψ[eq[:eq_z],  exo[:z_sh]] = 1

    ### 20. Uncertainty

    Γ0[eq[:eq_xi], endo[:xi_t]] = 1
    Γ1[eq[:eq_xi], endo[:xi_t]] = m[:ρ_xi]
    Ψ[eq[:eq_xi],  exo[:xi_sh]] = 1

    ### 21. Expected consumption

    Γ0[eq[:eq_Ec], endo[:c_t]]   = 1
    Γ1[eq[:eq_Ec], endo[:Ec_t1]] = 1
    Π[eq[:eq_Ec],  ex[:Ec_sh]]   = 1

    ### 22. Expected real interest rate

    Γ0[eq[:eq_Er], endo[:r_t]]   = 1
    Γ1[eq[:eq_Er], endo[:Er_t1]] = 1
    Π[eq[:eq_Er],  ex[:Er_sh]]   = 1

    return Γ0, Γ1, C, Ψ, Π
end
