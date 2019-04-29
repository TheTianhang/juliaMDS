struct Simulation
    allinfo::allinfo
    forcefields::Vector{ForceField}
    thermostat::ThermoStat
    #analyze_log::Analyze
    #output::Output
    steps::Int64
end

function run(simulation::Simulation)
    for stepii in 1:simulation.steps
        for forcefieldii in 1:length(simulation.forcefields)
            forcecompute(simulation.allinfo,simulation.forcefields[forcefieldii])
        end
        integrate(simulation.allinfo,simulation.thermostat)
        #if stepii % simulation.analyze_log.period == 0
        #    analyze(simulation.universe,simulation.analyze_log)
        #end
        #if stepii % simulation.output.period==0
        #    output(simulation.universe,simulation.output)
        #end
    end
end