
abstract type ForceField end

abstract type BondPotential<:ForceField
end
abstract type BondParameter end

struct HarmonicBondParam <: BondParameter
    bondtype::Int64
    K::Float64
    r_0::Float64
end

struct HarmonicBond <: BondPotential
    parameters::Vector{HarmonicBondParam}
end

function forcecompute(universe::Universe,bond::HarmonicBond)
    for i in 1:length(bond.parameters)
        parameter=bond.parameters[i]
        numbermoles=universe.numbermoles
        force=universe.forces
        coordinates=universe.coordinates
        K=parameter.K
        r_0=parameter.r_0
        bondtype=parameter.bondtype
        for moleii in 1:numbermoles
            numberbonds=length(universe.molecules[moleii].bonds)
            bonds=universe.molecules[moleii].bonds
            for bondii in 1:numberbonds
                if bonds[bondii].bondtype==bondtype
                    dx=coordinates[bonds[bondii].tag_a].x-coordinates[bonds[bondii].tag_b].x
                    dy=coordinates[bonds[bondii].tag_a].y-coordinates[bonds[bondii].tag_b].y
                    dz=coordinates[bonds[bondii].tag_a].z-coordinates[bonds[bondii].tag_b].z
                    d_2=dx*dx+dy*dy+dz*dz
                    d=sqrt(d_2)
                    potential=0.5*K*(d-r_0)*(d-r_0)
                    forcemag_divr=K * (r_0/d - 1.0)
                    force[bonds[bondii].tag_a].x += dx * forcemag_divr;
                    force[bonds[bondii].tag_a].y += dy * forcemag_divr;
                    force[bonds[bondii].tag_a].z += dz * forcemag_divr;

                    force[bonds[bondii].tag_b].x -= dx * forcemag_divr;
                    force[bonds[bondii].tag_b].y -= dy * forcemag_divr;
                    force[bonds[bondii].tag_b].z -= dz * forcemag_divr;
                end
            end
        end
    end
end