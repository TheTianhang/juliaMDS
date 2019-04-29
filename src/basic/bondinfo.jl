struct bond
    tag_a::Int64
    tag_b::Int64
end

function get_Nbonds(bonds::Vector{bond})
    return length(bonds)
end

function get_NbondTypes(tag_bond::Vector{Int64})
    bondtype=[]
    for bondii in 1:length(tag_bond)
        if tag_bond[bondii] in bondtype
            continue
        else
            push!(bondtype,tag_bond[bondii])
        end
    end
    return length(bondtype) 
end

function get_n_tag_bond(basicinfo::basicinfo,bonds::Vector{bond})
    Natoms=get_Natoms(basicinfo.coordinates)::Int64
    n_tag_bond=zeros(Int64,Natoms)
    Nbonds=get_Nbonds(bonds)::Int64
    for bondii in 1:Nbonds
        tag_a=bonds[bondii].tag_a
        tag_b=bonds[bondii].tag_b
        n_tag_bond[tag_a]=n_tag_bond[tag_a]+1
        n_tag_bond[tag_b]=n_tag_bond[tag_b]+1
    end
    return n_tag_bond
end
    
function get_max_n_tag_bond(basicinfo::basicinfo,bonds::Vector{bond})
    Natoms=get_Natoms(basicinfo.coordinates)::Int64
    n_tag_bond=get_n_tag_bond(basicinfo,bonds)
        
    MaxNumofBond=0::Int64
    for atomii in 1:Natoms
        if n_tag_bond[atomii]>MaxNumofBond
            MaxNumofBond=n_tag_bond[atomii]
        end
    end
    return MaxNumofBond
end
    
struct bondinfo
    bonds::Vector{bond}
    tag_bond::Vector{Int64}
    n_tag_bond::Vector{Int64}
    max_n_tag_bond::Int64
    Nbonds::Int64
    NbondTypes::Int64
end
    
function get_bondinfo(basicinfo::basicinfo,bonds::Vector{bond},tag_bond::Vector{Int64})
    Nbonds=get_Nbonds(bonds)
    NbondTypes=get_NbondTypes(tag_bond)
    n_tag_bond=get_n_tag_bond(basicinfo,bonds)
    max_n_tag_bond=get_max_n_tag_bond(basicinfo,bonds)
    return bondinfo(bonds,tag_bond,n_tag_bond,max_n_tag_bond,Nbonds,NbondTypes)
end
