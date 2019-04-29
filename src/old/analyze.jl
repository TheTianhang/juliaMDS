
abstract type Analyze end

struct analyze_log<:Analyze
    filename::String
    quantities::Vector{String}
    period::Int64
end

function analyze(universe::Universe,analyze::Analyze)
    c=1
end

