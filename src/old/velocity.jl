using Distributions

function maxwellboltzmann(mass::Real, T::Real)
    return rand(Normal(0.0, sqrt(8.3144598*T/mass)))
end

function Velocity(mass::Real, T::Real)
    return Velocity(maxwellboltzmann(mass, T),maxwellboltzmann(mass, T),maxwellboltzmann(mass, T))
end