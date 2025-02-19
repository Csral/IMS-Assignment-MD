import numpy as np
import mdtraj as md
import MDAnalysis as mda

traj = md.load("polyethylene.pdb");
u = mda.Universe("polyethylene.pdb");

dcharges = {
    "C": 0.1,
    "O": -0.3,
    "N": -0.2,
    "H": 0.1,
    "Cl": -0.3
};

charges = np.array([dcharges.get(atom[0], 0.0) for atom in u.atoms.names]);
bonds = [[bond.atom1.index, bond.atom2.index] for bond in traj.topology.bonds];

bondLength = md.compute_distances(traj, bonds)[0];
bondLength[bondLength < 0.1] = 0.1;

angles = [];

for i, bond1 in enumerate(bonds):
    for bond2 in bonds[i + 1:]:
        shared_atom = set(bond1) & set(bond2);
    
        if shared_atom:
            a, b = bond1;
            c, d = bond2;
    
            common = shared_atom.pop();
            angle = [a, common, d] if common == b else [c, common, b];
            if angle not in angles:
                angles.append(angle);

bondAngles = np.degrees(md.compute_angles(traj, angles)[0]);
bondAngles = np.nan_to_num(bondAngles, nan=109.5);

dihedralIndices = [];

for i, angle1 in enumerate(angles):
    for angle2 in angles[i + 1:]:
        shared_atoms = set(angle1) & set(angle2);

        if len(shared_atoms) == 2:
            a, b, c = angle1;
            d = angle2[2] if angle2[0] == b else angle2[0];
            dihedral = [a, b, c, d];
            
            if dihedral not in dihedralIndices:
                dihedralIndices.append(dihedral);

torsions = np.degrees(md.compute_dihedrals(traj, dihedralIndices)[0]);
forceConstants = {f"{u.atoms[b[0]].name}-{u.atoms[b[1]].name}": np.random.uniform(100, 500) for b in bonds};

bond_energy = lambda l, l0, k_b: 0.5 * k_b * (l - l0) ** 2
angle_energy = lambda theta, theta0, k_theta: 0.5 * k_theta * (theta - theta0) ** 2
torsion_energy = lambda omega, V_n, n, gamma: (V_n / 2) * (1 + np.cos(np.radians(n * omega) - gamma))
vdw_energy = lambda r, epsilon, Rij: epsilon * ((Rij / r) ** 12 - 2 * (Rij / r) ** 6)
electrostatic_energy = lambda q_i, q_j, r, epsilon_0: (q_i * q_j) / (4 * np.pi * epsilon_0 * max(r, 0.1))

E_bond = sum(bond_energy(l, 1.54, forceConstants.get(f"{u.atoms[b[0]].name}-{u.atoms[b[1]].name}", 300)) for b, l in zip(bonds, bondLength))
E_angle = sum(angle_energy(theta, 109.5, 120) for theta in bondAngles)
E_torsion = sum(torsion_energy(omega, 12, 3, 0) for omega in torsions)
E_vdw = sum(vdw_energy(r, 0.3, 0.35) for r in bondLength)
E_electrostatic = sum(electrostatic_energy(charges[b[0]], charges[b[1]], l, 8.854e-12) for b, l in zip(bonds, bondLength))
E_total = E_bond + E_angle + E_torsion + E_vdw + E_electrostatic

print(f"\n Main force field params:\n Bonds: {forceConstants}\n Angles: {len(angles)} detected\n Torsions: {len(dihedralIndices)} detected"
f"\n\nEnergy printed in components:\n Bond Energy: {E_bond:.4f} kJ/mol\n Angle Energy: {E_angle:.4f} kJ/mol\n Torsion Energy: {E_torsion:.4f} kJ/mol"
f"\n Van der Waals Energy: {E_vdw:.4f} kJ/mol\n Electrostatic Energy: {E_electrostatic:.4f} kJ/mol"
f"\n\nTotal Potential Energy: {E_total:.4f} kJ/mol")