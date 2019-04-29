
abstract type Output end

struct lammpstrj_output<:Output
    filename::String
    period::Int64
end

function output(universe::Universe,output::Output)
    c=1
end

