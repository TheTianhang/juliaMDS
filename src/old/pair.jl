

abstract type PairPotential<:ForceField
end
abstract type PairParameter end

abstract type NeighborList end

struct LennardJonesParam
    atomtype_a::Int64
    atomtype_b::Int64
    epsilon::Float64
    sigma::Float64
end

struct CellList <:NeighborList
end

function celllist(universe::Universe,cutpff_distance::Float64)
    
end

struct LennardJonesPair <: PairPotential
    neighborlist::NeighborList
    parameters::Vector{LennardJonesParam}
end

function forcecompute(universe::Universe,lennardjonespair::LennardJonesPair)

end