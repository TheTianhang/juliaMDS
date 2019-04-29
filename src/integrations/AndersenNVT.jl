using Random

struct AndersenNVT <: ThermoStat
    kT::Float64
    gamma::Float64
    seed::Int64
end

function set_AndersenNVT(;kT::Float64=1.0,gamma::Float64=1.0,seed::Int64=1111)
    return AndersenNVT(kT,gamma,seed)
end

#update coordinates of atoms in system
function AndersenNVT_firststep(allinfo::allinfo)
    forces=allinfo.basic.forces
    coordinates=allinfo.basic.coordinates
    velocities=allinfo.basic.velocities
    dt=allinfo.dt
    Natoms=allinfo.basic.Natoms
    for atomii in 1:Natoms
        coordinate=coordinates[atomii]
        pos_x=coordinate.x
        pos_y=coordinate.y
        pos_z=coordinate.z
        velocity=velocities[atomii]
        vel_x=velocity.vx
        vel_y=velocity.vy
        vel_z=velocity.vz
        force=forces[atomii]
        f_x=force.fx
        f_y=force.fy
        f_z=force.fz
        accel_x=f_x/coordinate.mass
        accel_y=f_y/coordinate.mass
        accel_z=f_z/coordinate.mass
        vel_x += 0.5 * accel_x * dt
        vel_y += 0.5 * accel_y * dt
        vel_z += 0.5 * accel_z * dt
        pos_x += vel_x * dt
        pos_y += vel_y * dt
        pos_z += vel_z * dt
        
        coordinates[atomii]=Coordinate(pos_x,pos_y,pos_z)
        velocities[atomii]=Velocity(vel_x,vel_y,vel_z)
    end
end

#update velocity under Andersen thermostat
function AndersenNVT_secondstep(allinfo::allinfo,kT::Float64,gamma::Float64,seed::Int64)
    forces=allinfo.basic.forces
    coordinates=allinfo.basic.coordinates
    velocities=allinfo.basic.velocities
    dt=allinfo.dt
    Natoms=allinfo.basic.Natoms
    pi2 = 2.0* acos(-1.0)
    Random.seed!(seed)
    for atomii in 1:Natoms
        coordinate=coordinates[atomii]
        velocity=velocities[atomii]
        vel_x=velocity.vx
        vel_y=velocity.vy
        vel_z=velocity.vz
        force=forces[atomii]
        f_x=force.fx
        f_y=force.fy
        f_z=force.fz
        accel_x=f_x/coordinate.mass
        accel_y=f_y/coordinate.mass
        accel_z=f_z/coordinate.mass
        rn=rand(1)[1]
        if rn<gamma
            sqmass=sqrt(kT//coordinate.mass)
            x1,y1,x2,y2,x3,y3=rand(6)
            dgauss2x = sqrt(-2.0*log(x1)) * cos(pi2 * y1)
            dgauss2y = sqrt(-2.0*log(x2)) * cos(pi2 * y2)
            dgauss2z = sqrt(-2.0*log(x3)) * cos(pi2 * y3)
            vel_x=sqmass*dgauss2x
            vel_y=sqmass*dgauss2y
            vel_z=sqmass*dgauss2z
        else
            vel_x += 0.5 * accel_x * dt
            vel_y += 0.5 * accel_y * dt
            vel_z += 0.5 * accel_z * dt
        end
        velocity=Velocity(vel_x,vel_y,vel_z)
        velocities[atomii]=velocity
    end
end

function integrate(allinfo::allinfo,andersenNVT::AndersenNVT)
    AndersenNVT_firststep(allinfo)
    AndersenNVT_secondstep(allinfo,andersenNVT.kT,andersenNVT.gamma,andersenNVT.seed)
end