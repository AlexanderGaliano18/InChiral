# Instalar OpenBabel (solo una vez por sesión en Colab)
# !apt-get install -y openbabel

import os
import zipfile
from google.colab import files

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
