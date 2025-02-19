import mdtraj as md
from openmm.app import PDBFile, ForceField
from openmm import Platform, VerletIntegrator
from openmm.unit import kilojoules_per_mole

pdb_file = '1uvi.pdb'
traj = md.load(pdb_file)

indices = []
print("Starting to compute bond lengths...")
ml = []

for i in range(0,100):
    for j in range(i+1,100):
        temp = []
        temp.append(i)
        temp.append(j)
        indices.append(temp)
        ml.append(temp)

dist = md.compute_distances(traj, ml)

print("Writing bond lengths...")

for i in range(0, len(dist)):
    
    msg = None
    try:
        msg = f"Atom Indices: ( {indices[i][0], indices[i][1]}) -----> {dist[i]}\n"
    except IndexError:
        try: 
            msg = f"Atom Indices: (?, ?) -----> {dist[i]}\n"
        except IndexError:
            msg = f"Atom Indices: (?, ?) -----> ERROR (or) unable to detect\n"

    with open("bond_length.txt", 'w') as out:
        out.write(msg)

print("Bond lengths computed!")
print("Starting to compute bond angles...")

indices = []
ml = []

for i in range(0,100):
    for j in range(i+1,100):
        for k in range(j+1, 100):
            temp = []
            temp.append(i)
            temp.append(j)
            temp.append(k)
            ml.append(temp)
            indices.append(temp)

bond_angle = md.compute_angles(traj, ml)
bond_angle_deg = []
for x in bond_angle:
    bond_angle_deg.append(x[0] * (180/3.14159265358979323846264)) #* To degrees

print("Writing bond angles...")

for i in range(0, len(bond_angle_deg)):
    
    msg = None
    try:
        msg = f"Atom Indices: ( {indices[i][0], indices[i][1], indices[i][2]}) -----> {bond_angle_deg[i]}\n"
    except IndexError:
        try: 
            msg = f"Atom Indices: (?, ?) -----> {bond_angle_deg[i]}\n"
        except IndexError:
            msg = f"Atom Indices: (?, ?) -----> ERROR (or) unable to detect\n"
    
    with open("bond_angles.txt", 'w') as out:
        out.write(msg)

print(f"Bond angle in degrees computed")

pdb = PDBFile(pdb_file)
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology)
platform = Platform.getPlatformByName('CPU')
integrator = VerletIntegrator(0.001)
simulation = openmm.app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
state = simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy()
print(f"Potential Energy: {potential_energy.value_in_unit(kilojoules_per_mole)} kJ/mol")