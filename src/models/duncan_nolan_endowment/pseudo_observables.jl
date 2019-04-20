function init_pseudo_observable_mappings!(m::DuncanNolanEndowment)

    pseudo_names = [:y_t, :ctot_t, :nw_t]

    # Create PseudoObservable objects
    pseudo = OrderedDict{Symbol,PseudoObservable}()
    for k in pseudo_names
        pseudo[k] = PseudoObservable(k)
    end

    # Fill in names and reverse transforms
    pseudo[:y_t].name = "Output Growth"
    pseudo[:y_t].longname = "Output Growth Per Capita"

    pseudo[:ctot_t].name     = "Consumption Growth"
    pseudo[:ctot_t].longname = "Consumption Growth per capita"

    pseudo[:nw_t].name = "Labor income"
    pseudo[:nw_t].longname = "Real growth in wages and salaries"

    # Add to model object
    m.pseudo_observable_mappings = pseudo
end
