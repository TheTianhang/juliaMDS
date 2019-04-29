
abstract type ThermoStat end

struct AndersenNVT <: ThermoStat
    kT::Float64
    gamma::Float64
    seed::Int64
end

#update coordinates of atoms in system
function AndersenNVT_firststep(universe::Universe)
    forces=universe.forces
    coordinates=universe.coordinates
    velocities=universe.velocities
    dt=universe.dt
    for atomii in 1:numberatoms
        coordinate=coordinates[atomii]
        pos_x=coordinate.x
        pos_y=coordinate.y
        pos_z=coordinate.z
        velocity=velocities[atomii]
        vel_x=velocity.x
        vel_y=velocity.y
        vel_z=velocity.z
        force=forces[atomii]
        f_x=force.x
        f_y=force.y
        f_z=force.z
        accel_x=f_x
        accel_y=f_y
        accel_z=f_z
        vel_x += 0.5 * accel_x * dt
        vel_y += 0.5 * accel_y * dt
        vel_z += 0.5 * accel_z * dt
        pos_x += vel_x * dt
        pos_y += vel_y * dt
        pos_z += vel_z * dt
        
        coordinates[atomii]=Coordinates(pos_x,pos_y,pos_z)
        velocities[atomii]=Velocity(vel_x,vel_y,vel_z)
    end
end

#update velocity under Andersen thermostat
function AndersenNVT_secondstep(universe::Universe,kT::Float64,prob::Float64,pi2::Int64)
    forces=universe.forces
    coordinates=universe.coordinates
    velocities=universe.velocities
    dt=universe.dt
    for atomii in 1:numberatoms
        velocity=velocities[atomii]
        vel_x=velocity.x
        vel_y=velocity.y
        vel_z=velocity.z
        force=forces[atomii]
        f_x=force.x
        f_y=force.y
        f_z=force.z
        accel_x=f_x
        accel_y=f_y
        accel_z=f_z
        rn=rand(1)[1]
        if rn<prob
            sqmass=sqrt(kT/1)
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

function integrate(universe::Universe,andersenNVT::AndersenNVT)
    AndersenNVT_firststep(universe)
    AndersenNVT_secondstep(universe,andersenNVT.kT,andersenNVT.gamma,andersenNVT.seed)
end