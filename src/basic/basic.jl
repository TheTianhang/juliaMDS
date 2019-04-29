const kB = 1.3806505e-23 #J/K Boltzmann Constant

const NA = 6.02214179e23 # AVOGADRO'S CONSTANT

const R =  8.3145 #J/mol*K

mutable struct coordinate
    x::Float64
    y::Float64
    z::Float64
    mass::Int64
end

mutable struct velocity
    vx::Float64
    vy::Float64
    vz::Float64
end

mutable struct image
    ix::Int64
    iy::Int64
    iz::Int64
end

struct tag_atom
    tag::Int64
end

mutable struct force
    fx::Float64
    fy::Float64
    fz::Float64
end

mutable struct virial
    vrx::Float64
    vry::Float64
    vrz::Float64
end

struct axis
    low::Float64
    high::Float64
end

struct box
    x::axis
    y::axis
    z::axis
end

function get_boxsize(boxsize::Float64)
    axissz=axis(-boxsize/2,boxsize/2)
    return box(axissz,axissz,axissz)
end

function get_boxsize(boxsize::Vector{Float64})
    if length(boxsize)==3
        axis_x=axis(-boxsize[1]/2,boxsize[1]/2)
        axis_y=axis(-boxsize[2]/2,boxsize[2]/2)
        axis_z=axis(-boxsize[3]/2,boxsize[3]/2)
        return box(axis_x,axis_y,axis_z)
    elseif length(boxsize)==6
        axis_x=axis(-boxsize[1],boxsize[2])
        axis_y=axis(-boxsize[3],boxsize[4])
        axis_z=axis(-boxsize[5],boxsize[6])
        return box(axis_x,axis_y,axis_z)
    else
        println("ERROR:get_boxsize::length of box size")
    end
end

function get_Natoms(coordinates::Vector{coordinate})
    return length(coordinates)
end
   
function get_NatomTypes(tags_atom::Vector{tag_atom})
    atomtype=[]
    for atomii in 1:length(tags_atom)
        if tags_atom[atomii] in atomtype
            continue
        else
            push!(atomtype,tags_atom[atomii])
        end
    end
    return length(atomtype) 
end        
        
struct basicinfo
    Natoms::Int64
    NatomTypes::Int64
    coordinates::Vector{coordinate}
    velocities::Vector{velocity}
    images::Vector{image}
    tags_atom::Vector{tag_atom}
    forces::Vector{force}
    virials::Vector{virial}
    box::box
end

function get_basicinfo(coordinates::Vector{coordinate},velocities::Vector{velocity},images::Vector{image},tags_atom::Vector{tag_atom},forces::Vector{force},virials::Vector{virial},boxsize)
    Natoms=get_Natoms(coordinates)
    NatomTypes=get_NatomTypes(tags_atom)
    box=get_boxsize(boxsize)
    return basicinfo(Natoms,NatomTypes,coordinates,velocities,images,tags_atom,forces,virials,box)
end