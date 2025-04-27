#!/bin/bash

# ==============================================================================
# GROMACS Bulk Structure Preparation Script (v3.1 - with top file update)
# ==============================================================================
#
# Description:
#   Generates a GROMACS bulk system configuration (.gro) and topology (.top, .itp)
#   for a given molecule specified by its SMILES string.
#   Automatically estimates an appropriate initial box size based on the
#   number of molecules and density estimated from the SMILES string using
#   a companion Python script (calc_cell_length.py).
#   Assigns GAFF2 parameters using AmberTools (antechamber, parmchk2) via acpype.
#   Packs molecules using Packmol.
#   Updates the molecule count in the acpype-generated GMX top file. <--- Added Feature
#   Logs all output (stdout & stderr) to a file while also displaying it.
#
# Dependencies:
#   - AmberTools (antechamber, parmchk2)
#   - acpype
#   - OpenBabel (obabel)
#   - Packmol
#   - GROMACS (gmx)
#   - Python 3.x with RDKit library
#   - calc_cell_length.py (Python script in the same directory as this script)
#
# Usage:
#   bash create_gmx_bulk_v3.1.sh "SMILES_STRING" MOLECULE_NAME OUTPUT_PREFIX [NUM_MOLECULES] [NET_CHARGE]
#
# Arguments:
#   $1: SMILES_STRING   - SMILES string of the molecule (enclose in quotes).
#   $2: MOLECULE_NAME   - Short, unique name for the molecule (e.g., residue name, max 3-4 chars recommended).
#   $3: OUTPUT_PREFIX   - Prefix for all generated output files and the working directory.
#   $4: NUM_MOLECULES   - Number of molecules to pack in the box (optional, default: 100).
#   $5: NET_CHARGE      - Net integer charge of the molecule (optional, default: 0).
#
# Example:
#   bash create_gmx_bulk_v3.1.sh "CCCCCCCCCC" DEC decane 100 0
#
# ==============================================================================

# --- Strict Mode ---
set -e # Exit immediately if a command exits with a non-zero status.
set -o pipefail # Causes a pipeline to return the exit status of the last command that exited with a non-zero status.

# --- Input Arguments ---
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 \"SMILES_STRING\" MOLECULE_NAME OUTPUT_PREFIX [NUM_MOLECULES, default=100] [NET_CHARGE, default=0]"
    exit 1
fi

SMILES="$1"
MOLECULE_NAME="$2"
OUTPUT_PREFIX="$3"
NUM_MOLECULES=${4:-100}
NET_CHARGE=${5:-0}

# --- Script Location and Python Script Path ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PYTHON_CALC_SCRIPT="${SCRIPT_DIR}/calc_cell_length.py"

if [ ! -f "$PYTHON_CALC_SCRIPT" ]; then
    echo "ERROR: Python script 'calc_cell_length.py' not found in the script directory: $SCRIPT_DIR" >&2
    exit 1
fi
if ! command -v python &> /dev/null; then
    echo "ERROR: Python command not found. Please ensure Python 3 is installed and in your PATH." >&2
    exit 1
fi

# --- Working Directory ---
WORKDIR="${OUTPUT_PREFIX}_bulk_prep"
echo "[Info] Creating working directory: $WORKDIR"
mkdir -p "$WORKDIR"
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create working directory: $WORKDIR" >&2
    exit 1
fi
cd "$WORKDIR" || { echo "ERROR: Failed to change directory to $WORKDIR" >&2; exit 1; }
echo "[Info] Changed working directory to: $(pwd)"


# --- Setup Logging ---
LOG_FILE="$(pwd)/${OUTPUT_PREFIX}_preparation.log"
[ -e "$LOG_FILE" ] && rm -f "$LOG_FILE" # Start with a fresh log file each time
exec > >(tee "$LOG_FILE") 2> >(tee -a "$LOG_FILE" >&2) # Log stdout/stderr to file and terminal
echo "[Info] Logging stdout and stderr to: $LOG_FILE"
echo "[Info] Script started at: $(date)"


# --- Calculate Box Size using Python script ---
echo "----------------------------------------------------"
echo "[Info] Calculating estimated box size for N=$NUM_MOLECULES molecules..."
PYTHON_OUTPUT=$(python "$PYTHON_CALC_SCRIPT" "$SMILES" "$NUM_MOLECULES" 2>&1)
PYTHON_EXIT_CODE=$?

if [ $PYTHON_EXIT_CODE -ne 0 ] || ! [[ "$PYTHON_OUTPUT" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "####################################################" >&2
    echo "ERROR: Failed to calculate box size using calc_cell_length.py." >&2
    echo "       SMILES: '$SMILES', N=$NUM_MOLECULES" >&2
    echo "       Python script path: $PYTHON_CALC_SCRIPT" >&2
    echo "       Python script output/error:" >&2
    echo "$PYTHON_OUTPUT" >&2
    echo "####################################################" >&2
    exit 1
fi

CALCULATED_BOX_NM="$PYTHON_OUTPUT"
echo "[Info] Calculated initial estimated box size: ${CALCULATED_BOX_NM} nm"
echo "----------------------------------------------------"


# --- Print Initial Setup ---
echo "--- Starting GROMACS Bulk Structure Preparation ---"
echo "Molecule Name:   $MOLECULE_NAME"
echo "SMILES:          $SMILES"
echo "Molecule Count:  $NUM_MOLECULES"
echo "Net Charge:      $NET_CHARGE"
echo "Output Prefix:   $OUTPUT_PREFIX"
echo "Working Dir:     $(pwd)"
echo "Est. Box Size:   $CALCULATED_BOX_NM nm x $CALCULATED_BOX_NM nm x $CALCULATED_BOX_NM nm"
echo "----------------------------------------------------"

# Define Step Count (updated)
TOTAL_STEPS=7

# === Step 1/7: Generate 3D structure from SMILES ===
echo "[Step 1/${TOTAL_STEPS}] Generating initial 3D structure from SMILES (using obabel)..."
obabel -:"$SMILES" --gen3d -omol2 -O "${OUTPUT_PREFIX}_initial.mol2" --title "$MOLECULE_NAME"
if [ $? -ne 0 ]; then echo "ERROR: obabel failed to generate 3D structure!" >&2; exit 1; fi
echo "[Step 1/${TOTAL_STEPS}] Completed: ${OUTPUT_PREFIX}_initial.mol2"

# === Step 2/7: Assign Atomic Charges and Atom Types ===
echo "[Step 2/${TOTAL_STEPS}] Assigning atomic charges (AM1-BCC) and GAFF2 atom types (using antechamber)..."
antechamber -i "${OUTPUT_PREFIX}_initial.mol2" -fi mol2 \
            -o "${OUTPUT_PREFIX}_charged.mol2" -fo mol2 \
            -c bcc -s 2 \
            -at gaff2 \
            -nc "$NET_CHARGE"
if [ $? -ne 0 ]; then echo "ERROR: antechamber failed! Check sqm.out and antechamber output for details." >&2; exit 1; fi
echo "[Step 2/${TOTAL_STEPS}] Completed: ${OUTPUT_PREFIX}_charged.mol2"

# === Step 3/7: Check/Generate Missing Force Field Parameters ===
echo "[Step 3/${TOTAL_STEPS}] Checking for missing GAFF2 parameters (using parmchk2)..."
parmchk2 -i "${OUTPUT_PREFIX}_charged.mol2" -f mol2 -o "${OUTPUT_PREFIX}.frcmod" -a Y -s gaff2
if [ $? -ne 0 ]; then echo "ERROR: parmchk2 failed!" >&2; exit 1; fi
echo "[Step 3/${TOTAL_STEPS}] Completed: ${OUTPUT_PREFIX}.frcmod"

# === Step 4/7: Convert Amber files to GROMACS Format ===
echo "[Step 4/${TOTAL_STEPS}] Converting topology to GROMACS format (using acpype)..."
# --frcmod オプションは acpype では認識されないため削除。
# acpype (内部の tleap) は通常、カレントディレクトリにある対応する .frcmod ファイルを自動で読み込む。
acpype -i "${OUTPUT_PREFIX}_charged.mol2" -o gmx -b "$OUTPUT_PREFIX" -n "$NET_CHARGE" -a gaff2 --verbose # <-- --frcmod を削除
if [ $? -ne 0 ]; then echo "ERROR: acpype execution reported an error! Check acpype logs/output." >&2; exit 1; fi

# Define expected output paths
ACPYPE_OUTPUT_DIR="${OUTPUT_PREFIX}.acpype"
SINGLE_MOL_GRO="${ACPYPE_OUTPUT_DIR}/${OUTPUT_PREFIX}_GMX.gro"
SINGLE_MOL_ITP="${ACPYPE_OUTPUT_DIR}/${OUTPUT_PREFIX}_GMX.itp"
# *** Define the top file path ***
TOP_FILE="${ACPYPE_OUTPUT_DIR}/${OUTPUT_PREFIX}_GMX.top" # The topology file generated by acpype

# Check if crucial acpype output files exist, including the .top file
if [ ! -d "$ACPYPE_OUTPUT_DIR" ]; then
    echo "ERROR: Expected output directory '$ACPYPE_OUTPUT_DIR' not found after running acpype." >&2
    exit 1
fi
# Check for necessary files including the .top file now
if [ ! -f "$SINGLE_MOL_GRO" ] || [ ! -f "$SINGLE_MOL_ITP" ] || [ ! -f "$TOP_FILE" ]; then
    echo "ERROR: Did not find expected .gro, .itp, and/or .top files in '$ACPYPE_OUTPUT_DIR'." >&2
    echo "       Please check the directory content and acpype output." >&2
    echo "       Looking for GRO: $SINGLE_MOL_GRO" >&2
    echo "       Looking for ITP: $SINGLE_MOL_ITP" >&2
    echo "       Looking for TOP: $TOP_FILE" >&2
    ls -l "$ACPYPE_OUTPUT_DIR" # List directory content for debugging
    exit 1
fi
echo "[Step 4/${TOTAL_STEPS}] Completed. GROMACS files generated in: $ACPYPE_OUTPUT_DIR"
echo "  Single Mol GRO: $SINGLE_MOL_GRO"
echo "  Single Mol ITP: $SINGLE_MOL_ITP"
echo "  GROMACS TOP:    $TOP_FILE"

# === Step 5/7: Create Bulk Structure using Packmol ===
echo "[Step 5/${TOTAL_STEPS}] Creating bulk structure with Packmol..."

# 5.1 Convert single molecule GRO to PDB for Packmol input
SINGLE_MOL_PDB="${OUTPUT_PREFIX}_single.pdb"
echo "[Info] Converting single molecule GRO to PDB for Packmol..."
gmx editconf -f "$SINGLE_MOL_GRO" -o "$SINGLE_MOL_PDB"
if [ $? -ne 0 ]; then echo "ERROR: gmx editconf (gro->pdb conversion) failed!" >&2; exit 1; fi
echo "[Info] Created $SINGLE_MOL_PDB"

# 5.2 Create Packmol input file (pack.inp)
PACKMOL_INPUT="pack.inp"
PACKMOL_OUTPUT_PDB="packed_temp.pdb"
BOX_SIZE_A=$(echo "$CALCULATED_BOX_NM * 10" | bc) # Convert nm to Angstrom for Packmol
echo "[Info] Creating Packmol input file ($PACKMOL_INPUT) with box size ${BOX_SIZE_A} Å..."

# Write Packmol input commands
echo "tolerance 2.0" > $PACKMOL_INPUT
echo "filetype pdb" >> $PACKMOL_INPUT
echo "output $PACKMOL_OUTPUT_PDB" >> $PACKMOL_INPUT
echo "" >> $PACKMOL_INPUT
echo "# Define the molecule to be packed" >> $PACKMOL_INPUT
echo "structure $SINGLE_MOL_PDB" >> $PACKMOL_INPUT
echo "  number $NUM_MOLECULES" >> $PACKMOL_INPUT
echo "  inside cube 0. 0. 0. $BOX_SIZE_A" >> $PACKMOL_INPUT
echo "end structure" >> $PACKMOL_INPUT

echo "[Info] Generated Packmol input file ($PACKMOL_INPUT):"
cat $PACKMOL_INPUT
echo "----------------------------------------"

# 5.3 Run Packmol
echo "[Info] Running Packmol..."
packmol < $PACKMOL_INPUT
if [ $? -ne 0 ]; then echo "ERROR: packmol failed! Check packmol output above." >&2; exit 1; fi
if [ ! -f "$PACKMOL_OUTPUT_PDB" ]; then
    echo "ERROR: Packmol finished but did not generate the expected output file: $PACKMOL_OUTPUT_PDB" >&2
    exit 1
fi
echo "[Step 5/${TOTAL_STEPS}] Completed: Packmol created $PACKMOL_OUTPUT_PDB"

# === Step 6/7: Generate Final GROMACS .gro File ===
echo "[Step 6/${TOTAL_STEPS}] Converting Packmol PDB output to final GROMACS .gro file..."
FINAL_BULK_GRO="${OUTPUT_PREFIX}_bulk_${NUM_MOLECULES}.gro"

# Use gmx editconf to convert PDB to GRO and set the box vectors
echo "[Info] Running gmx editconf to create $FINAL_BULK_GRO with box size $CALCULATED_BOX_NM nm..."
gmx editconf -f "$PACKMOL_OUTPUT_PDB" -o "$FINAL_BULK_GRO" -box $CALCULATED_BOX_NM $CALCULATED_BOX_NM $CALCULATED_BOX_NM -c
if [ $? -ne 0 ]; then echo "ERROR: gmx editconf (pdb->gro conversion) failed!" >&2; exit 1; fi
echo "[Step 6/${TOTAL_STEPS}] Completed: Final bulk structure file $FINAL_BULK_GRO created."

# === Step 7/7: Modify Molecule Count in GROMACS .top File ===
echo "[Step 7/${TOTAL_STEPS}] Modifying molecule count in GROMACS top file ($TOP_FILE)..."
if [ -f "$TOP_FILE" ]; then
    # sed コマンドで [ molecules ] セクションの行を探す。
    # 行頭の任意個の空白(^[\t ]*)、OUTPUT_PREFIX、1つ以上の空白([[:space:]]+)に続く
    # 数字の '1' で終わる行を見つけ、'1' を NUM_MOLECULES に置換する。
    # -i.bak でバックアップを作成しつつファイルを直接編集する。
    # パターン: 行頭任意空白(\(^[[:space:]]*\)) OUTPUT_PREFIX 1つ以上の空白(\(${OUTPUT_PREFIX}[[:space:]]+\)) 1 行末(\(1[[:space:]]*$\))
    # 置換: キャプチャしたグループ1(\1) キャプチャしたグループ2(\2) NUM_MOLECULES
    # よりシンプルな方法：OUTPUT_PREFIX を含む行で末尾の 1 を置換
    sed -i.bak "s/^\([[:space:]]*${OUTPUT_PREFIX}[[:space:]]\+\)1[[:space:]]*$/\1${NUM_MOLECULES}/" "$TOP_FILE" # <- この行を修正

    # Check if sed command was successful
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to modify molecule count in $TOP_FILE using sed." >&2
        echo "       Pattern may not have matched. Check the content of [ molecules ] section." >&2
        echo "       Attempting to restore from backup..." >&2
        if [ -f "${TOP_FILE}.bak" ]; then
            mv "${TOP_FILE}.bak" "$TOP_FILE"
            echo "       Restored from ${TOP_FILE}.bak" >&2
        else
            echo "       ERROR: Backup file ${TOP_FILE}.bak not found!" >&2
        fi
        exit 1 # Exit indicating failure
    else
        echo "[Info] Molecule count in $TOP_FILE successfully updated to $NUM_MOLECULES."
        # Remove the backup file if sed was successful
        rm -f "${TOP_FILE}.bak"
        echo "[Step 7/${TOTAL_STEPS}] Completed."
    fi
else
    # This case indicates an earlier problem if acpype was expected to create the file.
    echo "WARNING: Topology file $TOP_FILE not found! Cannot update molecule count." >&2
    echo "[Step 7/${TOTAL_STEPS}] Skipped (File not found)."
fi

# === Final Summary ===
echo "----------------------------------------------------"
echo "--- Preparation Complete! ---"
echo "Timestamp: $(date)"
echo "Outputs generated in directory: $(pwd)"
echo "----------------------------------------------------"

# Optional: Clean up intermediate files
# echo "Cleaning up intermediate files..."
# rm -f "${OUTPUT_PREFIX}_initial.mol2" "${OUTPUT_PREFIX}_charged.mol2" "${OUTPUT_PREFIX}.frcmod" "$SINGLE_MOL_PDB" "$PACKMOL_INPUT" "$PACKMOL_OUTPUT_PDB" ANTECHAMBER* ATOMTYPE.INF leap.log sqm.* NEWPDB.PDB PREPIN.* *.log *.bak

# Return to the original directory from where the script was called
cd ..
echo "[Info] Returned to directory: $(pwd)"

exit 0 # Indicate successful completion of the script