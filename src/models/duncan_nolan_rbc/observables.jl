function init_observable_mappings!(m::DuncanNolanRBC)

    observables = OrderedDict{Symbol,Observable}()
    population_mnemonic = get(get_setting(m, :population_mnemonic))

    ############################################################################
    ## 1. Real GDP Growth
    ############################################################################
    gdp_fwd_transform = function (levels)
        # FROM: Level of GDP (from FRED)
        # TO: Quarter-to-quarter percent change of real GDP per capita

        levels[:temp] = percapita(m, :GDP, levels)
        gdp = 1000 * nominal_to_real(:temp, levels)
        oneqtrpctchange(gdp)
    end

    gdp_rev_transform = loggrowthtopct_annualized_percapita

    observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, population_mnemonic, :GDPDEF__FRED],
                                       gdp_fwd_transform, gdp_rev_transform,
                                       "Real GDP Growth", "Real GDP Growth Per Capita")

   ############################################################################
   ## 1. Consumption growth
   ############################################################################
   con_fwd_transform = function (levels)
       # FROM: Level of GDP (from FRED)
       # TO: Quarter-to-quarter percent change of real GDP per capita

       levels[:temp] = percapita(m, :PCEC, levels)
       con = 1000 * nominal_to_real(:temp, levels)
       oneqtrpctchange(con)
   end

   gdp_rev_transform = loggrowthtopct_annualized_percapita

   observables[:obs_con] = Observable(:obs_con, [:PCEC__FRED, population_mnemonic, :GDPDEF__FRED],
                                      gdp_fwd_transform, gdp_rev_transform,
                                      "Real Consumption Growth",
                                      "Real Consumption Growth Per Capita")

  ############################################################################
  ## 4. Wages and salaries
  ############################################################################

  wagsal_fwd_transform = function (levels)
      # FROM: Level of Wages and Salaries (from FRED)
      # TO: Quarter-to-quarter percent change of real wages and salaries

      levels[:temp] = percapita(m, :A4102C1Q027SBEA, levels)
      wagsal = 1000 * nominal_to_real(:temp, levels)
      oneqtrpctchange(wagsal)
  end

  wagsal_rev_transform = loggrowthtopct_annualized_percapita

  observables[:obs_wagsal] = Observable(:obs_wagsal, [:A4102C1Q027SBEA__FRED, population_mnemonic, :GDPDEF__FRED],
                                     wagsal_fwd_transform, wagsal_rev_transform,
                                     "Real Labour Income Growth", "Real Labour Income Growth Per Capita")




    # ############################################################################
    # ## 2. CPI Inflation
    # ############################################################################
    #
    # cpi_fwd_transform = function (levels)
    #     # FROM: CPI urban consumers index (from FRED)
    #     # TO: Annualized quarter-to-quarter percent change of CPI index
    #
    #     quartertoannual(oneqtrpctchange(levels[:CPIAUCSL]))
    # end
    #
    # cpi_rev_transform = loggrowthtopct_annualized
    #
    # observables[:obs_cpi] = Observable(:obs_cpi, [:CPIAUCSL__FRED],
    #                                     cpi_fwd_transform, cpi_rev_transform,
    #                                     "CPI Inflation",
    #                                     "CPI Inflation")
    #
    # ############################################################################
    # ## 3. Nominal short-term interest rate (3 months)
    # ############################################################################
    #
    # nominalrate_fwd_transform = function (levels)
    #     # FROM: Nominal effective federal funds rate (aggregate daily data at a
    #     #       quarterly frequency at an annual rate)
    #     # TO:   Nominal effective fed funds rate, at a quarterly rate
    #
    #     levels[:DFF]
    # end
    #
    # nominalrate_rev_transform = identity
    #
    # observables[:obs_nominalrate] = Observable(:obs_nominalrate, [:DFF__FRED],
    #                                            nominalrate_fwd_transform, nominalrate_rev_transform,
    #                                            "Nominal FFR",
    #                                            "Nominal Effective Fed Funds Rate")


    m.observable_mappings = observables
end
