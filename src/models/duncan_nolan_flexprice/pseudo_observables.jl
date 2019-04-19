function init_pseudo_observable_mappings!(m::DuncanNolanFlexPrice)

    pseudo_names = [:y_t, :π_t, :z_t, :NominalFFR, :RealFFR, :nw_t]

    # Create PseudoObservable objects
    pseudo = OrderedDict{Symbol,PseudoObservable}()
    for k in pseudo_names
        pseudo[k] = PseudoObservable(k)
    end

    # Fill in names and reverse transforms
    pseudo[:y_t].name = "Output Growth"
    pseudo[:y_t].longname = "Output Growth Per Capita"

    pseudo[:π_t].name = "Inflation"
    pseudo[:π_t].longname = "Inflation"
    pseudo[:π_t].rev_transform = quartertoannual

    pseudo[:z_t].name     = "z_t"
    pseudo[:z_t].longname = "z_t"

    pseudo[:NominalFFR].name     = "Nominal FFR"
    pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"

    pseudo[:RealFFR].name     = "Real FFR"
    pseudo[:RealFFR].longname = "Real FFR at an annual rate"

    pseudo[:nw_t].name = "Labor income"
    pseudo[:nw_t].longname = "Real growth in wages and salaries"

    # Add to model object
    m.pseudo_observable_mappings = pseudo
end
