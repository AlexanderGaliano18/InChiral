# ========================
# 1. Instalar RDKit en Colab
# ========================
#!pip install rdkit

# ========================
# 2. Importar librerías
# ========================
from rdkit import Chem
from rdkit.Chem import AllChem
from google.colab import files
import os, zipfile

# ========================
# 3. Subir archivo.smi
# ========================
uploaded = files.upload()  # selecciona tu archivo.smi
input_file = "archivo.smi"

output_folder = "xyz_files_rdkit"
zip_name = "xyz_results_rdkit.zip"
os.makedirs(output_folder, exist_ok=True)

# ========================
# 4. Función para convertir SMILES → XYZ
# ========================
def smiles_to_xyz(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Error: SMILES inválido {smiles}")
        return

    # Agregar hidrógenos
    mol = Chem.AddHs(mol)

    # Generar geometría inicial con ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 42  # reproducible
    if AllChem.EmbedMolecule(mol, params) != 0:
        print(f"⚠️ No se pudo generar conformación 3D para {smiles}")
        return

    # Optimizar con MMFF94 (si está disponible)
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.UFFOptimizeMolecule(mol)

    # Guardar en formato XYZ
    conf = mol.GetConformer()
    with open(filename, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n{smiles}\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")

# ========================
# 5. Procesar cada línea del archivo.smi
# ========================
with open(input_file, "r") as f:
    smiles_list = [line.strip() for line in f if line.strip()]

print(f"🔎 Se encontraron {len(smiles_list)} moléculas")

for i, smi in enumerate(smiles_list, start=1):
    out_file = os.path.join(output_folder, f"mol_{i}.xyz")
    smiles_to_xyz(smi, out_file)

# ========================
# 6. Comprimir en ZIP y descargar
# ========================
with zipfile.ZipFile(zip_name, "w") as zipf:
    for file in os.listdir(output_folder):
        zipf.write(os.path.join(output_folder, file), file)

files.download(zip_name)

print(f"\n📦 Archivo '{zip_name}' generado con éxito 🎉")
