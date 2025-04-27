#!/bin/bash

# --- 0. ファイルパスと名前の確認 (必要なら編集) ---

INITIAL_GRO="decane_bulk_prep/decane_bulk_100.gro"
TOPOLOGY_FILE="decane_bulk_prep/decane.acpype/decane_GMX.top"
MDP_MINIM="min.mdp"
MDP_NVT="nvt.mdp"
MDP_NPT="npt.mdp"
MDP_PROD="md_prod.mdp"

# --- 1. エネルギー最小化 (Energy Minimization) ---
echo "Starting Energy Minimization (EM)..."
gmx grompp -f ${MDP_MINIM} -c ${INITIAL_GRO} -p ${TOPOLOGY_FILE} -o em.tpr -maxwarn 1
if [ $? -ne 0 ]; then echo "ERROR during grompp for EM. Exiting."; exit 1; fi

gmx mdrun -v -deffnm em
# --- GPUを使う場合の例 (EMでは効果は限定的) ---
# gmx mdrun -v -deffnm em -nb gpu
if [ $? -ne 0 ]; then echo "ERROR during mdrun for EM. Exiting."; exit 1; fi
echo "EM finished successfully."
echo "----------------------------------------"

# --- 2. NVT 平衡化 ---
echo "Starting NVT Equilibration..."
gmx grompp -f ${MDP_NVT} -c em.gro -r em.gro -p ${TOPOLOGY_FILE} -o nvt.tpr -maxwarn 1
#  ^^^^^^ EM後の座標(em.gro)を入力に使う
if [ $? -ne 0 ]; then echo "ERROR during grompp for NVT. Exiting."; exit 1; fi

gmx mdrun -v -deffnm nvt
# --- GPUを使う場合の例 (推奨) ---
# gmx mdrun -v -deffnm nvt -nb gpu -pme gpu
if [ $? -ne 0 ]; then echo "ERROR during mdrun for NVT. Exiting."; exit 1; fi
echo "NVT Equilibration finished successfully."
echo "----------------------------------------"

# --- 3. NPT 平衡化 ---
echo "Starting NPT Equilibration..."
# NVTの状態(座標、速度、ボックスサイズ等)を引き継ぐために -t nvt.cpt を使う
gmx grompp -f ${MDP_NPT} -c nvt.gro -r nvt.gro -t nvt.cpt -p ${TOPOLOGY_FILE} -o npt.tpr -maxwarn 1
#  ^^^^^^ NVT後の座標(nvt.gro)とチェックポイント(nvt.cpt)を入力に使う
if [ $? -ne 0 ]; then echo "ERROR during grompp for NPT. Exiting."; exit 1; fi

gmx mdrun -v -deffnm npt
# --- GPUを使う場合の例 (推奨) ---
# gmx mdrun -v -deffnm npt -nb gpu -pme gpu
if [ $? -ne 0 ]; then echo "ERROR during mdrun for NPT. Exiting."; exit 1; fi
echo "NPT Equilibration finished successfully."
echo "----------------------------------------"

# --- 4. 本計算 (Production MD) ---
echo "Starting Production MD run..."
echo "*** WARNING: This step can take a very long time! ***"
# NPTの状態を引き継ぐために -t npt.cpt を使う
gmx grompp -f ${MDP_PROD} -c npt.gro -t npt.cpt -p ${TOPOLOGY_FILE} -o md_prod.tpr -maxwarn 1
#  ^^^^^^ NPT後の座標(npt.gro)とチェックポイント(npt.cpt)を入力に使う
if [ $? -ne 0 ]; then echo "ERROR during grompp for Production MD. Exiting."; exit 1; fi

gmx mdrun -v -deffnm md_prod
# --- GPUを使う場合の例 (強く推奨) ---
# gmx mdrun -v -deffnm md_prod -nb gpu -pme gpu
if [ $? -ne 0 ]; then echo "ERROR during mdrun for Production MD. Exiting."; exit 1; fi
echo "Production MD finished successfully."
echo "----------------------------------------"
echo "You can now calculate viscosity using the command:"
echo "gmx energy -f md_prod.edr -vis viscosity_gk.xvg -evisco viscosity_einstein.xvg -o energy.xvg -b 2000"
#!/bin/bash

# ==============================================================================
# GROMACS Simulation Workflow Script
# ==============================================================================
#
# Description:
#   Runs a standard GROMACS simulation workflow:
#   1. Energy Minimization (EM)
#   2. NVT Equilibration
#   3. NPT Equilibration
#   4. Production MD
#   Takes the initial structure (.gro) and topology (.top) files as input arguments.
#
# Dependencies:
#   - GROMACS (gmx command)
#   - Corresponding .mdp files (min.mdp, nvt.mdp, npt.mdp, md_prod.mdp) in the
#     directory where the script is run.
#
# Usage:
#   bash run_gmx_simulation.sh <initial_structure.gro> <topology.top>
#
# Arguments:
#   $1: Initial Structure File (.gro) - Path to the starting coordinate file.
#   $2: Topology File (.top)         - Path to the GROMACS topology file.
#
# Example:
#   bash run_gmx_simulation.sh decane_bulk_prep/decane_bulk_100.gro decane_bulk_prep/decane.acpype/decane_GMX.top
#
# ==============================================================================

# --- Strict Mode ---
set -e # Exit immediately if a command exits with a non-zero status.
set -o pipefail # Causes a pipeline to return the exit status of the last command that exited with a non-zero status.

# --- 0. Input Arguments and File Definitions ---

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <initial_structure.gro> <topology.top>"
    echo "Error: Incorrect number of arguments provided."
    echo "Example: $0 system.gro system.top"
    exit 1
fi

# Assign arguments to variables
INITIAL_GRO="$1"
TOPOLOGY_FILE="$2"

# Check if the input files exist
if [ ! -f "${INITIAL_GRO}" ]; then
    echo "Error: Initial GRO file not found: ${INITIAL_GRO}" >&2
    exit 1
fi
if [ ! -f "${TOPOLOGY_FILE}" ]; then
    echo "Error: Topology TOP file not found: ${TOPOLOGY_FILE}" >&2
    exit 1
fi

# Define MDP file names (assuming they exist in the current directory)
MDP_MINIM="min.mdp"
MDP_NVT="nvt.mdp"
MDP_NPT="npt.mdp"
MDP_PROD="md_prod.mdp"

# Check if MDP files exist
for mdp in ${MDP_MINIM} ${MDP_NVT} ${MDP_NPT} ${MDP_PROD}; do
    if [ ! -f "$mdp" ]; then
        echo "Error: MDP file not found: $mdp" >&2
        echo "Please ensure min.mdp, nvt.mdp, npt.mdp, and md_prod.mdp exist in the current directory." >&2
        exit 1
    fi
done

echo "--- Starting GROMACS Simulation Workflow ---"
echo "Initial GRO file: ${INITIAL_GRO}"
echo "Topology file:    ${TOPOLOGY_FILE}"
echo "MDP files:        ${MDP_MINIM}, ${MDP_NVT}, ${MDP_NPT}, ${MDP_PROD}"
echo "----------------------------------------"

# --- 1. エネルギー最小化 (Energy Minimization) ---
echo "Starting Energy Minimization (EM)..."
# Use input arguments for -c and -p
gmx grompp -f ${MDP_MINIM} -c ${INITIAL_GRO} -p ${TOPOLOGY_FILE} -o em.tpr -maxwarn 2
if [ $? -ne 0 ]; then echo "ERROR during grompp for EM. Exiting." >&2; exit 1; fi

echo "Running mdrun for EM..."
gmx mdrun -v -deffnm em
# --- GPUを使う場合の例 (EMでは効果は限定的) ---
# gmx mdrun -v -deffnm em -nb gpu
if [ $? -ne 0 ]; then echo "ERROR during mdrun for EM. Exiting." >&2; exit 1; fi
echo "EM finished successfully."
echo "Output files: em.gro, em.log, em.edr, em.trr (optional)"
echo "----------------------------------------"

# --- 2. NVT 平衡化 ---
echo "Starting NVT Equilibration..."
# Use EM output (em.gro) for coordinates and topology file argument for -p
gmx grompp -f ${MDP_NVT} -c em.gro -r em.gro -p ${TOPOLOGY_FILE} -o nvt.tpr -maxwarn 2
if [ $? -ne 0 ]; then echo "ERROR during grompp for NVT. Exiting." >&2; exit 1; fi

echo "Running mdrun for NVT..."
gmx mdrun -v -deffnm nvt
# --- GPUを使う場合の例 (推奨) ---
# gmx mdrun -v -deffnm nvt -nb gpu -pme gpu -bonded gpu
if [ $? -ne 0 ]; then echo "ERROR during mdrun for NVT. Exiting." >&2; exit 1; fi
echo "NVT Equilibration finished successfully."
echo "Output files: nvt.gro, nvt.log, nvt.edr, nvt.trr, nvt.cpt"
echo "----------------------------------------"

# --- 3. NPT 平衡化 ---
echo "Starting NPT Equilibration..."
# Use NVT output (nvt.gro, nvt.cpt) for coordinates/state and topology file argument for -p
gmx grompp -f ${MDP_NPT} -c nvt.gro -r nvt.gro -t nvt.cpt -p ${TOPOLOGY_FILE} -o npt.tpr -maxwarn 2
if [ $? -ne 0 ]; then echo "ERROR during grompp for NPT. Exiting." >&2; exit 1; fi

echo "Running mdrun for NPT..."
gmx mdrun -v -deffnm npt
# --- GPUを使う場合の例 (推奨) ---
# gmx mdrun -v -deffnm npt -nb gpu -pme gpu -bonded gpu
if [ $? -ne 0 ]; then echo "ERROR during mdrun for NPT. Exiting." >&2; exit 1; fi
echo "NPT Equilibration finished successfully."
echo "Output files: npt.gro, npt.log, npt.edr, npt.trr, npt.cpt"
echo "----------------------------------------"

# --- 4. 本計算 (Production MD) ---
echo "Starting Production MD run..."
echo "*** WARNING: This step can take a very long time! ***"
# Use NPT output (npt.gro, npt.cpt) for coordinates/state and topology file argument for -p
gmx grompp -f ${MDP_PROD} -c npt.gro -t npt.cpt -p ${TOPOLOGY_FILE} -o md_prod.tpr -maxwarn 2
if [ $? -ne 0 ]; then echo "ERROR during grompp for Production MD. Exiting." >&2; exit 1; fi

echo "Running mdrun for Production MD..."
gmx mdrun -v -deffnm md_prod
# --- GPUを使う場合の例 (強く推奨) ---
# gmx mdrun -v -deffnm md_prod -nb gpu -pme gpu -bonded gpu
if [ $? -ne 0 ]; then echo "ERROR during mdrun for Production MD. Exiting." >&2; exit 1; fi
echo "Production MD finished successfully."
echo "Output files: md_prod.gro, md_prod.log, md_prod.edr, md_prod.xtc, md_prod.cpt"
echo "----------------------------------------"

# --- Final Message ---
echo "Simulation workflow completed."
echo "You can now analyze the results, for example, calculate viscosity:"
echo "gmx energy -f md_prod.edr -evisco viscosity_einstein.xvg -o energy.xvg -b 2000" # Example command
echo "python plot_viscosity.py viscosity_einstein.xvg" # Example command to plot viscosity
echo "----------------------------------------"
exit 0