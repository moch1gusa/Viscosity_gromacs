1. create environments:

conda create -n env_gro -c conda-forge ambertools openbabel
conda activate env_gro
pip install acpype rdkit


2. create molecule bulk:


Usage:
bash create_gmx_bulk_v4.sh "SMILES_STRING" MOLECULE_NAME OUTPUT_PREFIX [NUM_MOLECULES] [NET_CHARGE]


ex. 100 molecules decane

bash 01_create_gmx_bulk_v4.sh "CCCCCCCCCC" DEC decane 100 0


3. EMD and calculate viscosity :

minimize -> nvt -> npt -> nvt


bash 02_vis_calc_gro.sh
