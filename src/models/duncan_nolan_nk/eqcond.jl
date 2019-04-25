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
function eqcond(m::DuncanNolanNK)
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

    Γ0[eq[:eq_euler], endo[:c_t]]   =  1
    Γ0[eq[:eq_euler], endo[:R_t]]   =  1/m[:τ]
    Γ0[eq[:eq_euler], endo[:Ec_t1]] = -1
    Γ0[eq[:eq_euler], endo[:Eπ_t1]] = -1/m[:τ]

    ### 2. Labour supply

    Γ0[eq[:eq_labsup], endo[:n_t]] = -m[:psi]
    Γ0[eq[:eq_labsup], endo[:w_t]] =  1
    Γ0[eq[:eq_labsup], endo[:c_t]] = -m[:τ]

    ### 3. Production

    Γ0[eq[:eq_prod], endo[:y_t]] = -1
    Γ0[eq[:eq_prod], endo[:z_t]] =  1
    Γ0[eq[:eq_prod], endo[:n_t]] =  (1-m[:alpha])

    ### 4. Entrepreneurs: Consumption savings

    Γ0[eq[:eq_entcon], endo[:ce_t]]      = -1
    Γ0[eq[:eq_entcon], endo[:omegae_t]]  =  1
    Γ0[eq[:eq_entcon], endo[:re_t]]      =  1

    ### 5. Entrepreneurs: Leverage

    Γ0[eq[:eq_entlev], endo[:lev_t]]     = -1
    Γ0[eq[:eq_entlev], endo[:y_t]]       =  1
    Γ0[eq[:eq_entlev], endo[:r_t]]       = -1
    Γ0[eq[:eq_entlev], endo[:ce_t]]      = -1
    Γ0[eq[:eq_entlev], endo[:re_t]]      =  1

    ### 6. Entrepreneurs: Wedge

    Γ0[eq[:eq_entweg], endo[:lev_t]] = -1
    Γ0[eq[:eq_entweg], endo[:tau_t]] =  m[:L]
    Γ0[eq[:eq_entweg], endo[:xi_t]]  = -(1+m[:erp])

    ### 7. Entrepreneurs: Wealth evolution

    Γ0[eq[:eq_entwel], endo[:ce_t]]  =  1
    Γ0[eq[:eq_entwel], endo[:re_t]]  = -1
    Γ1[eq[:eq_entwel], endo[:ce_t]]  =  1

    Γ0[eq[:eq_entwel], endo[:c_t]]  =  -m[:wopen]*m[:τ]
    Γ0[eq[:eq_entwel], endo[:r_t]]  =   m[:wopen]
    Γ1[eq[:eq_entwel], endo[:c_t]]  =  -m[:wopen]*m[:τ]

    ### 8. Factor Prices: Capital

    Γ0[eq[:eq_fperp], endo[:re_t]]   = -1
    Γ0[eq[:eq_fperp], endo[:r_t]]    =  1
    Γ0[eq[:eq_fperp], endo[:lev_t]]  =  m[:erp]
    Γ0[eq[:eq_fperp], endo[:tau_t]]  =  m[:L]

    ### 9. Factor Prices: Labor

    Γ0[eq[:eq_fplab], endo[:w_t]]   = -1
    Γ0[eq[:eq_fplab], endo[:y_t]]   =  1
    Γ0[eq[:eq_fplab], endo[:n_t]]   = -1
    Γ0[eq[:eq_fplab], endo[:tau_t]] = -1

    ### 10. Fisher relation

    Γ0[eq[:eq_fisher], endo[:r_t]]   = -1
    Γ1[eq[:eq_fisher], endo[:R_t]]   =  1
    Γ0[eq[:eq_fisher], endo[:π_t]]   = -1

    ### 11. NK Phillips Curve

    Γ0[eq[:eq_phillips], endo[:w_t]]   = -m[:κ]
    Γ0[eq[:eq_phillips], endo[:y_t]]   = m[:κ]
    Γ0[eq[:eq_phillips], endo[:n_t]]   = -m[:κ]
    Γ0[eq[:eq_phillips], endo[:tau_t]]   = -m[:κ]
    Γ0[eq[:eq_phillips], endo[:π_t]]   =  1
    Γ0[eq[:eq_phillips], endo[:Eπ_t1]] = -1/(1+m[:rA]/400)

    ### 12. Monetary Policy Rule

    Γ0[eq[:eq_mp], endo[:y_t]] = -(1-m[:ρ_R])*m[:ψ_2]
    Γ0[eq[:eq_mp], endo[:π_t]] = -(1-m[:ρ_R])*m[:ψ_1]
    Γ0[eq[:eq_mp], endo[:R_t]] =  1
    Γ0[eq[:eq_mp], endo[:g_t]] =  (1-m[:ρ_R])*m[:ψ_2]
    Γ1[eq[:eq_mp], endo[:R_t]] =  m[:ρ_R]
    Ψ[eq[:eq_mp],  exo[:rm_sh]] =  1

    ### 13. Wages and salaries lag

    Γ0[eq[:eq_nw], endo[:nw_t]] = -1
    Γ0[eq[:eq_nw], endo[:w_t]]  =  1
    Γ0[eq[:eq_nw], endo[:n_t]]  =  1

    ### 13. Wages and salaries lag

    Γ0[eq[:eq_nw_t1], endo[:nw_t1]] = 1
    Γ1[eq[:eq_nw_t1], endo[:nw_t]]  = 1

    ### 13. Output lag

    Γ0[eq[:eq_y_t1], endo[:y_t1]] = 1
    Γ1[eq[:eq_y_t1], endo[:y_t]]  = 1

    ### 14. Government spending

    Γ0[eq[:eq_g], endo[:g_t]] = 1
    Γ1[eq[:eq_g], endo[:g_t]] = m[:ρ_g]
    Ψ[eq[:eq_g],  exo[:g_sh]] = 1

    ### 15. Technology

    Γ0[eq[:eq_z], endo[:z_t]] = 1
    Γ1[eq[:eq_z], endo[:z_t]] = m[:ρ_z]
    Ψ[eq[:eq_z],  exo[:z_sh]] = 1

    ### 20. Uncertainty

    Γ0[eq[:eq_xi], endo[:xi_t]] = 1
    Γ1[eq[:eq_xi], endo[:xi_t]] = m[:ρ_xi]
    Ψ[eq[:eq_xi],  exo[:xi_sh]] = 1

    ### 16. Expected output

    Γ0[eq[:eq_Ec], endo[:c_t]]   = 1
    Γ1[eq[:eq_Ec], endo[:Ec_t1]] = 1
    Π[eq[:eq_Ec],  ex[:Ec_sh]]   = 1

    ### 17. Expected inflation

    Γ0[eq[:eq_Eπ], endo[:π_t]]   = 1
    Γ1[eq[:eq_Eπ], endo[:Eπ_t1]] = 1
    Π[eq[:eq_Eπ],  ex[:Eπ_sh]]   = 1

    return Γ0, Γ1, C, Ψ, Π
end
