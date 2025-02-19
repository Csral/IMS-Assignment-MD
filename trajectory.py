import os
import subprocess
import MDAnalysis as mda

def analyze_structure(pdb_file):
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    u = mda.Universe(pdb_file)
    print("Structure Analysis:")
    print(f"Number of atoms: {len(u.atoms)}")
    print(f"Number of residues: {len(u.residues)}")
    print("Residue names:", list(set(u.residues.resnames)))
    return u

def generate_gromacs_topology(pdb_file, force_field="amber99sb-ildn"):
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    base_name = os.path.splitext(pdb_file)[0]
    gro_file = f"{base_name}.gro"
    #topol.top = f"{base_name}.top"
    boxed_gro = "boxed.gro"
    solvated_gro = "solvated.gro"
    ionized_gro = "ionized.gro"
    em_tpr = f"{base_name}.tpr"
    edr_file = f"{base_name}.edr"
    cpt_file = "em_tpr.cpt"

    #MDP parameters for energy minimization
    mdp_content = """
    ; Energy Minimization Parameters
    integrator = steep       ; Steepest descent minimization
    emtol = 1000.0           ; Convergence criterion (kJ/mol/nm)
    emstep = 0.01            ; Minimization step size (nm)
    nsteps = 50000           ; Maximum number of steps
    nstenergy = 1            ; Write energies every step
    """
    # Save to a file
    with open("em.mdp", "w") as f:
        f.write(mdp_content)
    print("MDP file created: em.mdp")

    # Minimal MDP file for ion addition
    ions_mdp_content = """
    ; Dummy MDP file for ion addition
    integrator = steep
    nsteps = 1
    """
    with open("ions.mdp", "w") as f:
        f.write(ions_mdp_content)
    print("Dummy MDP file for ion addition created: ions.mdp")

    nvt_mdp = """
    ; NVT Equilibration Parameters
    integrator = md          ; Leap-frog integrator
    dt = 0.002               ; Time step (ps)
    nsteps = 50000           ; Total simulation time: 100 ps
    nstxout = 500            ; Save coordinates every 1 ps
    nstvout = 500            ; Save velocities every 1 ps
    nstenergy = 500          ; Save energies every 1 ps

    ; Temperature Coupling
    tcoupl = v-rescale       ; Temperature coupling method
    tc-grps = System         ; Apply to the entire system
    tau_t = 0.1              ; Coupling time constant (ps)
    ref_t = 300              ; Reference temperature (K)

    ; Constraints
    constraints = h-bonds    ; Constrain all bonds involving hydrogen
    constraint_algorithm = LINCS
    """
    with open("nvt.mdp", "w") as f:
        f.write(nvt_mdp)
    print("MDP file created: nvt.mdp")

    npt_mdp = """
    ; NPT Equilibration Parameters
    integrator = md          ; Leap-frog integrator
    dt = 0.002               ; Time step (ps)
    nsteps = 50000           ; Total simulation time: 100 ps
    nstxout = 500            ; Save coordinates every 1 ps
    nstvout = 500            ; Save velocities every 1 ps
    nstenergy = 500          ; Save energies every 1 ps

    ; Temperature Coupling
    tcoupl = v-rescale       ; Temperature coupling method
    tc-grps = System         ; Apply to the entire system
    tau_t = 0.1              ; Coupling time constant (ps)
    ref_t = 300              ; Reference temperature (K)

    ; Pressure Coupling
    pcoupl = Parrinello-Rahman
    pcoupltype = isotropic   ; Isotropic pressure coupling
    tau_p = 2.0              ; Coupling time constant (ps)
    ref_p = 1.0              ; Reference pressure (bar)
    compressibility = 4.5e-5 ; Compressibility (bar^-1), typical for water

    ; Constraints
    constraints = h-bonds    ; Constrain all bonds involving hydrogen
    constraint_algorithm = LINCS
    """
    with open("npt.mdp", "w") as f:
        f.write(npt_mdp)
    print("MDP file created: npt.mdp")

    try:
        # Check if GROMACS is installed
        subprocess.run(["gmx", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Step 1: Run pdb2gmx to generate topology and coordinate files
        print("Running pdb2gmx...")
        subprocess.run(
            ["gmx", "pdb2gmx", "-f", pdb_file, "-o", gro_file, "-p", "topol.top", "-ff", force_field],
            check=True
        )
        print(f"Topology and coordinate files generated: {gro_file}, topol.top")

        # Step 2: Create a cubic box
        print("Creating a cubic box...")
        subprocess.run(
            ["gmx", "editconf", "-f", gro_file, "-o", boxed_gro, "-c", "-d", "1.0", "-bt", "cubic"],
            check=True
        )
        print(f"Cubic box created: {boxed_gro}")

        # Step 3: Solvate the system
        print("Solvating the system...")
        subprocess.run(
            ["gmx", "solvate", "-cp", boxed_gro, "-cs", "spc216.gro", "-o", solvated_gro, "-p", "topol.top"],
            check=True
        )
        print(f"Solvated system created: {solvated_gro}")

        # Step 4: Add ions to neutralize the system
        print("Adding ions...")
        subprocess.run(
            ["gmx", "grompp", "-f", "ions.mdp", "-c", solvated_gro, "-p", "topol.top", "-o", "ions.tpr"],
            check=True
        )
        subprocess.run(
            ["gmx", "genion", "-s", "ions.tpr", "-o", ionized_gro, "-p", "topol.top", "-pname", "NA", "-nname", "CL", "-neutral"],
            input=b"SOL\n",  # Automatically select water group for ion replacement
            check=True
        )
        print(f"Ions added: {ionized_gro}")

        # Step 5: Preprocess for energy minimization
        print("Preprocessing for energy minimization...")
        subprocess.run(
            ["gmx", "grompp", "-f", "em.mdp", "-c", ionized_gro, "-p", "topol.top", "-o", "em.tpr"],
            check=True
        )
        print(f"Energy minimization input file created: em.tpr")

        # Step 6: Run the simulation
        print("Running energy minimization...")
        subprocess.run(
            ["gmx", "mdrun", "-deffnm", "em"],
            check=True
        )
        print(f"Simulation completed. EDR file generated: em_edr")

        # Step 7: Energy Analysis
        print("Energy Analysis:")
        subprocess.run(
            ["gmx", "energy", "-f", "em.edr", "-o", "potential.xvg"],
            check=True
        )

        # Step 8: Equilibration
        print("Equilibration:")
        subprocess.run(
            ["gmx", "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro", "-p", "topol.top", "-o", "nvt.tpr"],
            check=True
        )
        print(f"Equilibration input file created: nvt.tpr")

        # Step 9: Run the simulation
        print("Running equilibration...")
        subprocess.run(
            ["gmx", "mdrun", "-v", "-deffnm", "nvt"],
            check=True
        )
        print(f"Equilibration completed. EDR file generated: nvt.edr")
        
        # Step 10: Analysis of equilibration (Temperature)
        print("Equilibration Analysis:")
        subprocess.run(
            ["gmx", "energy", "-f", "nvt.edr", "-o", "temperature.xvg"],
            check=True
        )

        # Step 11: Complete Equilibrium
        print("Equilibration:")
        subprocess.run(
            ["gmx", "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro", "-t", "nvt.cpt", "-p", "topol.top", "-o", "npt.tpr"],
            check=True
        )
        subprocess.run(
            ["gmx", "mdrun", "-deffnm", "npt"],
            check=True
        )

        # Step 12: Analysis of equilibration (Pressure)
        print("Equilibration Analysis:")
        subprocess.run(
            ["gmx", "energy", "-f", "npt.edr", "-o", "pressure.xvg"],
            check=True
        )

        # Step 13: Analysis of equilibration (Density)
        print("Equilibration Analysis:")
        subprocess.run(
            ["gmx", "energy", "-f", "npt.edr", "-o", "density.xvg"],
            check=True
        )

        # Step 14: Production MD
        print("Production MD:")
        subprocess.run(
            ["gmx", "grompp", "-f", "md.mdp", "-c", "npt.gro", "-t", "npt.cpt", "-p", "topol.top", "-o", "md_0_1.tpr"],
            check=True
        )

        # Step 15: Execute the MD production
        print("Running MD production...")
        subprocess.run(
            ["gmx", "mdrun", "-deffnm", "md_0_1"],
            check=True
        )

        # Step 16: Trajectory Simulation
        print("Trajectory Simulation:")
        subprocess.run(
            ["gmx", "trjconv", "-s", "md_0_1.tpr", "-f", "md_0_1.xtc", "-o", "md_0_1_noPBC.xtc", "-pbc", "mol", "-center"], # Enter 1 and 0
            check=True
        )

        # Step 17: RMSD Analysis
        print("RMSD Analysis:")
        subprocess.run(
            ["gmx", "rms", "-s", "md_0_1.tpr", "-f", "md_0_1_noPBC.xtc", "-o", "rmsd.xvg", "-tu", "ns"], # Enter 4 and 4
            check=True
        )

        # Step 18: Radius of Gyration:
        print("Radius of Gyration:")
        subprocess.run(
            ["gmx", "gyrate", "-s", "md_0_1.tpr", "-f", "md_0_1_noPBC.xtc", "-o", "gyrate.xvg"], # Enter 1
            check=True
        )

    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        print("Check the output above for details.")
    except FileNotFoundError:
        print("GROMACS 'gmx' command not found. Please install GROMACS and add it to your PATH.")

def main():
    # Specify the path to the local PDB file
    pdb_file = "4f51.pdb"  # Replace with the path to your local PDB file

    # Step 1: Analyze the structure
    universe = analyze_structure(pdb_file)

    # Step 2: Generate GROMACS topology and run energy minimization
    force_field = "amber99sb-ildn"  # Example force field
    generate_gromacs_topology(pdb_file, force_field)

if __name__ == "__main__":
    main()