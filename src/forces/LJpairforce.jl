
struct LJpairParam <: PairParameter
    epsilon::Float64
    sigma::Float64
    cut_off::FLoat64
end

struct LJpair
    dict::Dict
end

function init_LJpair()
    return LJpair(Dict())
end

function set_LJpair(LJpair::LJpair;tag_a::Int64=1,tag_b::Int64=2,epsilon::Float64=1.0,sigma::Float64=1.0,cut_off::Float64=1.0)
    if tag_a<=tag_b
        tag=String("$tag_a")*"-"*String("$tag_b")
    else
        tag=String("$tag_b")*"-"*String("$tag_a")
    end
    get!(LJpair.dict,tag,LJpairParam(epsilon,sigma,cut_off))
end

function forcecompute(allinfo::allinfo,pair::LJpair)
    forces=allinfo.basicinfo.forces
    coordinates=allinfo.basicinfo.coordinates
    bonds=allinfo.bondinfo.bonds
    tags_atom=allinfo.basicinfo.tags_atom
    for atomii in 1:allinfo.basicinfo.Natoms-1
        for atomjj in atomii:allinfo.basicinfo.Natoms
            tag_ii=tags_atom[atomii]
            tag_jj=tags_atom[atomjj]
            if tag_ii<=tag_jj
                tag=String("$tag_ii")*"-"*String("$tag_jj")
            else
                tag=String("$tag_jj")*"-"*String("$tag_ii")
            end
            parameter=pair.dict[tag]
        #bondtype=allinfo.bondinfo.tag_bond[bondii]
            vd=vector3d(coordinates[atomii],coordinates[atomjj])
            d_2=sqr(vd)
            d=sqrt(d_2)
            if d_2>(parameter.cut_off)^2
                d2inv = 0.0
			    dcut2inv = 0.0
            else
                d2inv = 1.0 / d_2
            end
            
            d6inv = d2inv * d2inv * d2inv;
		    force_divr= d2inv * d6inv * (12.0 * 4*parameter.epsilon*(parameter.sigma)^12  * d6inv - 6.0 * 4*parameter.epsilon*(parameter.sigma)^6)
            
            forces[atomii].fx += dx * force_divr
            forces[atomii].fy += dy * force_divr
            forces[atomii].fz += dz * force_divr

            forces[atomjj].fx -= dx * force_divr
            forces[atomjj].fy -= dy * force_divr
            forces[atomjj].fz -= dz * force_divr
    
        end
    end
end