import itertools
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import io
import zipfile
import os

# ========================
# Función 1: Generar estereoisómeros
# ========================
def generar_estereoisomeros(smiles: str):
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i + 1] == "@":
                posiciones.append((i, True))  # ya es @@
                i += 2
            else:
                posiciones.append((i, False))  # es @ simple
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


# ========================
# Función 2: SMILES → XYZ
# ========================
def smiles_to_xyz(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None

    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    xyz_str = f"{mol.GetNumAtoms()}\n{smiles}\n"
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz_str += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"

    return xyz_str


# ========================
# Interfaz Streamlit
# ========================
st.title("🌀 Inchiral - Estereoisómeros y Conversión XYZ")
st.write("Ingresa un **SMILES** para detectar centros quirales, generar estereoisómeros y exportarlos en formatos `.smi` y `.xyz`.")

smiles = st.text_input("👉 Ingresa el código SMILES:", "")

if smiles:
    isomeros, mensaje = generar_estereoisomeros(smiles)
    st.info(mensaje)

    if isomeros:
        st.subheader("Ejemplos de estereoisómeros:")
        for s in isomeros[:5]:
            st.code(s, language="text")

        # Descargar archivo .smi
        smi_buffer = io.StringIO("\n".join(isomeros))
        st.download_button(
            "📥 Descargar archivo .smi",
            smi_buffer.getvalue(),
            file_name="estereoisomeros.smi",
            mime="text/plain"
        )

        # Convertir a XYZ y empaquetar en ZIP
        xyz_zip = io.BytesIO()
        with zipfile.ZipFile(xyz_zip, "w") as zipf:
            for i, s in enumerate(isomeros, start=1):
                xyz_data = smiles_to_xyz(s)
                if xyz_data:
                    zipf.writestr(f"mol_{i}.xyz", xyz_data)
        xyz_zip.seek(0)

        st.download_button(
            "📦 Descargar ZIP con XYZ",
            xyz_zip,
            file_name="xyz_results.zip",
            mime="application/zip"
        )
