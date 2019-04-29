
struct HarmonicBondParam <: BondParameter
    K::Float64
    r_0::Float64
end

struct harmonicbond
    dict::Dict
end

function init_harmonicbond()
    return harmonicbond(Dict())
end

function set_HarmonicBond(bondinit::harmonicbond;tag::Int64=1,K::Float64=1.0,r_0::Float64=1.0)
    get!(bondinit.dict,tag,HarmonicBondParam(K,r_0))
end

function forcecompute(allinfo::allinfo,bond::harmonicbond)
    forces=allinfo.basicinfo.forces
    coordinates=allinfo.basicinfo.coordinates
    bonds=allinfo.bondinfo.bonds
    for bondii in 1:allinfo.bondinfo.Nbonds
        bondtype=allinfo.bondinfo.tag_bond[bondii]
        parameter=bond.dict[bondtype]

        vd=vector3d(coordinates[bonds[bondii].tag_a],coordinates[bonds[bondii].tag_b])
        d_2=sqr(vd)
        d=sqrt(d_2)
        potential=0.5*parameter.K*(d-parameter.r_0)*(d-parameter.r_0)
        forcemag_divr=parameter.K * (parameter.r_0/d - 1.0)
        forces[bonds[bondii].tag_a].fx += dx * forcemag_divr;
        forces[bonds[bondii].tag_a].fy += dy * forcemag_divr;
        forces[bonds[bondii].tag_a].fz += dz * forcemag_divr;

        forces[bonds[bondii].tag_b].fx -= dx * forcemag_divr;
        forces[bonds[bondii].tag_b].fy -= dy * forcemag_divr;
        forces[bonds[bondii].tag_b].fz -= dz * forcemag_divr;
    end
end