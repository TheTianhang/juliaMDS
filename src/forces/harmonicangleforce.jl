
struct HarmonicAngleParam <: AngleParameter
    K::Float64
    r_0::Float64
end

struct harmonicangle
    dict::Dict
end

function init_harmonicangle()
    return harmonicangle(Dict())
end

function set_HarmonicAngle(bondinit::harmonicangle;tag::Int64=1,K::Float64=1.0,r_0::Float64=1.0)
    get!(bondinit.dict,tag,HarmonicAngleParam(K,r_0))
end

function forcecompute(allinfo::allinfo,angle::harmonicangle)
    forces=allinfo.basicinfo.forces
    coordinates=allinfo.basicinfo.coordinates
    angles=allinfo.angleinfo.angles
    for angleii in 1:allinfo.angleinfo.Nangles
        angletype=allinfo.angleinfo.tag_angle[angleii]
        parameter=angle.dict[angletype]
        d_ba=vector3d(coordinates[angles[angleii].tag_b],coordinates[angles[angleii].tag_a])
        d_bc=vector3d(coordinates[angles[angleii].tag_b],coordinates[angles[angleii].tag_c])

        pa=

        d_2=dx*dx+dy*dy+dz*dz
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