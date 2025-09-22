import streamlit as st
import itertools
import os
import zipfile
import subprocess

# ----------------------------
# Función: Generar estereoisómeros
# ----------------------------
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

    combinaciones = list(itertools.product(["@", "@@"], repeat=len(posiciones)))
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
    return resultados


# ----------------------------
# Interfaz Streamlit
# ----------------------------
st.set_page_config(page_title="InChiral", page_icon="🧬", layout="centered")
st.title("🧬 Generador de Estereoisómeros - InChiral")

smiles = st.text_input("👉 Ingresa el código SMILES:")

if st.button("Generar"):
    if smiles.strip() == "":
        st.error("⚠️ Debes ingresar un SMILES válido.")
    else:
        # Generar isómeros
        isomeros = generar_estereoisomeros(smiles)
        st.success(f"✅ Total estereoisómeros generados: {len(isomeros)}")

        # Mostrar algunos ejemplos
        st.write("Ejemplos generados:")
        st.code("\n".join(isomeros[:5]))

        # Guardar archivo .smi
        filename = "archivo.smi"
        with open(filename, "w") as f:
            for s in isomeros:
                f.write(s + "\n")

        # Carpeta de salida
        output_folder = "xyz_files"
        os.makedirs(output_folder, exist_ok=True)

        # Ejecutar OpenBabel con subprocess
        try:
            subprocess.run(
                ["obabel", filename, "-O", f"{output_folder}/isomero.xyz", "--gen3D", "-m"],
                check=True
            )
        except Exception as e:
            st.error("❌ Error ejecutando OpenBabel. Asegúrate de que esté instalado.")
            st.exception(e)
        else:
            # Comprimir resultados
            zip_name = "xyz_results.zip"
            with zipfile.ZipFile(zip_name, "w") as zipf:
                for file in os.listdir(output_folder):
                    zipf.write(os.path.join(output_folder, file), file)

            # Botón de descarga
            with open(zip_name, "rb") as f:
                st.download_button("📦 Descargar ZIP con isómeros", f, file_name=zip_name)
