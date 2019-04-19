"""
`init_subspec!(m::DuncanNolanFlexPrice)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::DuncanNolanFlexPrice)
    if subspec(m) == "ss0"
        return
    else
        error("This subspec should be a 0")
    end 
end

