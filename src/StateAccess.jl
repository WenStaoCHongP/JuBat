"""
    representative_temperature(case, state; supplied=nothing)

Return the temperature consumed by an electrochemical model.  The thermal
configuration, rather than the incidental presence of a key in `case.index`,
defines where that value comes from.
"""
function representative_temperature(case::Case, state; supplied=nothing)
    supplied === nothing || return supplied

    thermalmodel = case.opt.thermalmodel
    thermalmodel == "none" && return case.param.cell.T0
    if thermalmodel == "lumped" || thermalmodel == "distributed2D"
        return only(state[case.index["temperature"]])
    end

    throw(ArgumentError("thermal model $(repr(thermalmodel)) has no electrochemical temperature-state rule"))
end
