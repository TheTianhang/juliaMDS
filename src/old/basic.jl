

const kB = 1.3806505e-23 #J/K Boltzmann Constant

const NA = 6.02214179e23 # AVOGADRO'S CONSTANT

const R =  8.3145 #J/mol*K

struct Atom
    atomtype::Int64
    molindex::Int64
    charge::Float64
    mass::Float64
end

struct Bond
    bondtype::Int64
    tag_a::Int64
    tag_b::Int64
end

struct Angle
    angletype::Int64
    tag_a::Int64
    tag_b::Int64
    tag_c::Int64
end

struct Dihedral
    dihedraltype::Int64
    tag_a::Int64
    tag_b::Int64
    tag_c::Int64
    tag_d::Int64
end

struct Molecule
    molindex::Int64
    numberatomspermole::Int64
    atoms::Vector{Atom}
    numberbondtype::Int64
    bonds::Vector{Bond}
    numberangletype::Int64
    angles::Vector{Angle}
    numberdihedraltype::Int64
    dihedrals::Vector{Dihedral}
end

"3D coordinates, e.g. for an atom, in nm."
mutable struct Coordinates
    x::Float64
    y::Float64
    z::Float64
end

"3D velocity values, e.g. for an atom, in nm/ps."
mutable struct Velocity
    x::Float64
    y::Float64
    z::Float64
end

mutable struct Force
    x::Float64
    y::Float64
    z::Float64
end

struct Universe
    numbermoles::Int64
    numberatoms::Int64
    molecules::Vector{Molecule}
    coordinates::Vector{Coordinates}
    velocities::Vector{Velocity}
    forces::Vector{Force}
    temperature::Float64
    box_size::Vector{Float64}
    dt::Float64
    #neighbour_list::Vector{Tuple{Int, Int, Bool}} # i, j and whether they are 1-4 pairs (halved force)
end