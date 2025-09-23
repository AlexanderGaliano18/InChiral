import streamlit as st
import itertools
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import zipfile
import os

# -------------------
# 1. Generar estereoisómeros
# -------------------
def generar_estereoisomeros(smiles: str):
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                posiciones.append((i, True))  # @@
                i += 2
            else:
                posiciones.append((i, False)) # @
                i += 1
        else:
            i += 1

    n = len(posiciones)
    if n == 0:
        return [], "⚠️ El SMILES no tiene centros quirales."
    elif n > 3:
        return [], "❌ El SMILES tiene más de 3 centros quirales. No se generarán isómeros."

    combinaciones = list(itertools.product(["@", "@@"], repeat=n))

    resultados = []
    for comb in combinaciones:
        chars = list(smiles)
        offset = 0
        for (pos, era_doble), val in zip(posiciones, comb):
            real_pos = pos + offset
            if era_doble:
                chars[real_pos:real_pos+2] = list(val)
                offset += len(val) - 2
            else:
                chars[real_pos:real_pos+1] = list(val)
                offset += len(val) - 1
        resultados.append("".join(chars))

    return resultados, f"✅ Se generaron {len(resultados)} estereoisómeros."


# -------------------
# 2. SMILES → XYZ
# -------------------
def smiles_to_xyz(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        return False

    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    with open(filename, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n{smiles}\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")
    return True


# -------------------
# 3. Visualización 3D con py3Dmol
# -------------------
def visualizar_3d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)
    block = Chem.MolToMolBlock(mol)

    viewer = py3Dmol.view(width=400, height=300)
    viewer.addModel(block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer


# -------------------
# 4. Interfaz Streamlit
# -------------------
st.title("🧪 Generador de Estereoisómeros y Visualización 3D")
st.write("Ingresa un **código SMILES** para generar estereoisómeros y descargarlos en formato XYZ.")

smiles_input = st.text_input("👉 Ingresa el SMILES:", "CC(C)Br")

if st.button("Generar"):
    isomeros, mensaje = generar_estereoisomeros(smiles_input)
    st.info(mensaje)

    if isomeros:
        st.write("### Ejemplos de estereoisómeros generados:")
        for s in isomeros[:5]:
            st.code(s)

        # Visualización del primero
        viewer = visualizar_3d(isomeros[0])
        if viewer:
            viewer.show()
            st.components.v1.html(viewer._make_html(), height=350)

        # Guardar todos en carpeta y comprimir
        output_folder = "xyz_files"
        zip_name = "isomeros_xyz.zip"
        os.makedirs(output_folder, exist_ok=True)

        for i, smi in enumerate(isomeros, start=1):
            out_file = os.path.join(output_folder, f"mol_{i}.xyz")
            smiles_to_xyz(smi, out_file)

        with zipfile.ZipFile(zip_name, "w") as zipf:
            for file in os.listdir(output_folder):
                zipf.write(os.path.join(output_folder, file), file)

        with open(zip_name, "rb") as f:
            st.download_button("📥 Descargar ZIP con todos los .xyz", f, file_name=zip_name)
