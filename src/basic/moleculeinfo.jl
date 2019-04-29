struct moleculeinfo
    particle_mole_id::Vector{Int64}
    Nmoles::Int64
end

function get_moleculeinfo(basicinfo::basicinfo,bondinfo::bondinfo)
    Natoms=get_Natoms(basicinfo.coordinates)::Int64
    particle_mole_id=zeros(Natoms)
    m_NumMol=0::Int64
    n_tag_bond=get_n_tag_bond(basicinfo,bondinfo.bonds)
    #for atomii in 1:Natoms
    #    N_b = n_tag_bond[atomii]
    #    
    #end
    particle_mole_id=zeros(Int64,Natoms)
    Nmoles=1
    return moleculeinfo(particle_mole_id,Nmoles)
end