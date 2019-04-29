
universe=JuliaMD.universe()

harmonicbond1=JuliaMD.harmonicparameter(1,2,3)
harmonicbond2=JuliaMD.harmonicparameter(2,2,3)
harmonicbond=JuliaMD.HarmonicBond([harmonicbond1,harmonicbond2])

thermo=JuliaMD.AndersenNVT(kT=1.0,gamma=0.7,seed=1234)

analyze_log=JuliaMD.analyze_log(filename="test.log",quantities=["potential_energy"],period=100)

output=Julia.lammpstrj("lammpstrj.custom",period=100)

simulation=JuliaMD.Simulation(universe,
                            [harmonicbond,harmonicangle,oplsdihedral,LJ,hPF],
                            thermo,analyze_log,output,
                            steps=10000000)

@run simulation