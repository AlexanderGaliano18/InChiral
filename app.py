import os
import zipfile
import itertools
from google.colab import files

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
    print(f"🔎 Se encontraron {n} centros quirales (@).")

    # Verificación: solo aceptar exactamente 3
    if n > 3:
        print("❌ El SMILES tiene más de 3 centros quirales. No se generarán isómeros.")
        return []

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

    return resultados


# -------------------
# INTERACTIVO
# -------------------
smiles = input("👉 Ingresa el código SMILES: ")

isomeros = generar_estereoisomeros(smiles)

if isomeros:
    print(f"\n✅ Total estereoisómeros generados: {len(isomeros)}")
    print("Ejemplos:")
    for s in isomeros[:5]:
        print(s)


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
    print(f"🔎 Se encontraron {n} centros quirales (@).")

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

    return resultados


# -------------------
# INTERACTIVO
# -------------------
smiles = input("👉 Ingresa el código SMILES: ")

isomeros = generar_estereoisomeros(smiles)

print(f"\n✅ Total estereoisómeros generados: {len(isomeros)}")
print("Ejemplos:")
for s in isomeros[:5]:
    print(s)

# Guardar archivo
filename = "archivo.smi"
with open(filename, "w") as f:
    for s in isomeros:
        f.write(s + "\n")

print(f"\n📁 Archivo '{filename}' guardado con éxito.")

# Descargar archivo en Colab
files.download(filename)




# Subir archivo.smi desde tu PC
uploaded = files.upload()  # Selecciona "archivo.smi"

input_file = "archivo.smi"
output_folder = "xyz_files"
zip_name = "xyz_results.zip"

# Crear carpeta de salida
os.makedirs(output_folder, exist_ok=True)

# Generar un archivo .xyz por cada SMILES usando OpenBabel
# -m = múltiple salida (un archivo por molécula)
!obabel {input_file} -O {output_folder}/isomero.xyz --gen3D -m

# Comprimir los XYZ en un ZIP
with zipfile.ZipFile(zip_name, "w") as zipf:
    for file in os.listdir(output_folder):
        zipf.write(os.path.join(output_folder, file), file)

# Descargar el ZIP
files.download(zip_name)

print(f"\n📦 Archivo '{zip_name}' generado con éxito y listo para descargar.")
