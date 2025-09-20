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
