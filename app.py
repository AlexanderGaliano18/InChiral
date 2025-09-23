import itertools
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import os, zipfile
from io import BytesIO
import py3Dmol

# ========================
# 1. Generador de estereoisómeros
# ========================
def generar_estereoisomeros(smiles: str):
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                posiciones.append((i, True))  # ya es @@
                i += 2
            else:
                posiciones.append((i, False)) # es @ simple
                i += 1
        else:
            i += 1

    n = len(posiciones)
    if n == 0:
        return [], "⚠️ El SMILES no tiene centros quirales. No se generarán isómeros."
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

    return resultados, f"✅ Total estereoisómeros generados: {len(resultados)}"

# ========================
# 2. SMILES → XYZ con RDKit
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
    content = f"{mol.GetNumAtoms()}\n{smiles}\n"
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        content += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
    return content, mol

# ========================
# 3. Visualización 3D con py3Dmol
# ========================
def mostrar_molecula(mol):
    block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer

# ========================
# 4. Streamlit App
# ========================
st.title("🧪 Generador de estereoisómeros y visualizador 3D")

smiles = st.text_input("👉 Ingresa un código SMILES:")

if st.button("Generar"):
    if not smiles.strip():
        st.warning("Por favor, ingresa un SMILES válido.")
    else:
        isomeros, msg = generar_estereoisomeros(smiles)
        st.info(msg)

        if isomeros:
            st.write("Ejemplos de estereoisómeros:", isomeros[:5])

            # Crear ZIP en memoria
            memory_file = BytesIO()
            with zipfile.ZipFile(memory_file, "w") as zipf:
                for i, smi in enumerate(isomeros, start=1):
                    xyz_content, mol = smiles_to_xyz(smi)
                    if xyz_content:
                        zipf.writestr(f"mol_{i}.xyz", xyz_content)

                        # Mostrar solo la primera molécula como preview
                        if i == 1 and mol:
                            st.subheader("👀 Vista previa 3D del primer isómero")
                            viewer = mostrar_molecula(mol)
                            viewer.show()
                            st.components.v1.html(viewer._make_html(), height=420)

            memory_file.seek(0)

            st.download_button(
                label="📦 Descargar ZIP con XYZ",
                data=memory_file,
                file_name="xyz_results.zip",
                mime="application/zip"
            )
