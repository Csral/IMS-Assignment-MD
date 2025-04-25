# IMS-Assignment-MD

**IMS-Assignment-MD** is a Python-based Molecular Dynamics (MD) simulation toolkit designed to automate and streamline the setup, execution, and analysis of molecular dynamics simulations. It integrates with GROMACS and RDKit to facilitate tasks such as system preparation, simulation parameterization, and trajectory analysis.

---

## Features

- **Automated Simulation Setup**: Scripts to prepare molecular systems, including solvation, ion addition, and energy minimization.
- **Parameter File Management**: Includes `.mdp` files for various simulation stages: energy minimization (`em.mdp`), equilibration (`nvt.mdp`, `npt.mdp`), and production runs (`md.mdp`).
- **Trajectory Analysis**: Tools to process and analyze simulation trajectories, extracting meaningful insights.
- **RDKit Integration**: Utilizes RDKit for molecular manipulations and property calculations.
- **Logging and Output Management**: Organized directories (`out_logs`, `graphs`) to store simulation outputs and logs.

---

## Directory Structure

```

IMS-Assignment-MD/
├── build/               # Build scripts and related files
├── excess/              # Additional resources or data
├── graphs/              # Generated plots and graphs
├── out_logs/            # Simulation output logs
├── rdkit/               # RDKit-related scripts and data
├── 2zwj.pdb             # Sample PDB structure file
├── 2zwjc.pdb            # Sample PDB structure file
├── em.mdp               # Energy minimization parameters
├── ions.mdp             # Ion addition parameters
├── md.mdp               # Production MD parameters
├── mdout.mdp            # Output MD parameters
├── npt.mdp              # NPT equilibration parameters
├── nvt.mdp              # NVT equilibration parameters
├── prepare_build.sh     # Shell script to prepare the build
├── main.py              # Main execution script
├── non_gromacs.py       # Scripts not related to GROMACS
├── prop.py              # Property calculation scripts
├── trajectory.py        # Trajectory analysis scripts
└── xvgpy.py             # XVG file processing scripts
```


---

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/Csral/IMS-Assignment-MD.git
   cd IMS-Assignment-MD
   ```


2. **Set Up a Virtual Environment** (Optional but recommended):
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```


3. **Install Dependencies**


   *Note: Ensure that GROMACS is installed and accessible via the command line.*

---

## Usage

1. **Prepare the System**:
   Use the provided PDB files (`2zwj.pdb`, `2zwjc.pdb`) or your own to set up the molecular system.

2. **Run the Main Script**:
   ```bash
   python main.py
   ```


   This script will guide you through the simulation setup, execution, and analysis processes.

3. **Analyze Results**:
   Post-simulation, use the scripts in `trajectory.py` and `xvgpy.py` to analyze trajectories and generate plots.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- [GROMACS](http://www.gromacs.org/) for molecular dynamics simulations.
- [RDKit](https://www.rdkit.org/) for cheminformatics functionalities.

---

For more details and updates, visit the [IMS-Assignment-MD GitHub Repository](https://github.com/Csral/IMS-Assignment-MD).

--- 
