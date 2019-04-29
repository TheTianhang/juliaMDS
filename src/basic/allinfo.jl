function get_InitTimeStep(input::input)
    return input.inittimestep
end

struct allinfo
    basicinfo::basicinfo
    bondinfo::bondinfo
    #angleinfo::angleinfo
    moleculeinfo::moleculeinfo
    inittimestep::Int64
    dt::Float64
end

function get_allinfo(input::input,dt::Float64)
    Natoms=get_Natoms(input.coordinates)
    forces=force[]
    virials=virial[]
    for atomii in 1:Natoms
        push!(forces,force(0.0,0.0,0.0))
        push!(virials,virial(0.0,0.0,0.0))
    end
    velocities=velocity[]
    images=image[]
    for atomii in 1:Natoms
        push!(velocities,velocity(0.0,0.0,0.0))
        push!(images,image(0,0,0))
    end
    coordinates=input.coordinates
    tags_atom=input.tags_atom
    basicinfo=get_basicinfo(coordinates,velocities,images,tags_atom,forces,virials,input.boxsize)
    bondinfo=get_bondinfo(basicinfo,input.bonds,input.tag_bond)
    #angleinfo=get_angleinfo(basicinfo,input.angles,input.tag_angle)
    moleculeinfo=get_moleculeinfo(basicinfo,bondinfo)
    inittimestep=get_InitTimeStep(input)
    return allinfo(basicinfo,bondinfo,moleculeinfo,inittimestep,dt)
end