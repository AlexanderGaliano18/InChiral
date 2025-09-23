import streamlit as st
import itertools
import os
import zipfile
import tempfile
import io
import sys

# Configuración para evitar warnings de RDKit
import warnings
warnings.filterwarnings('ignore')

# Manejo de importación de RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    st.error("❌ RDKit no está instalado. Por favor instala RDKit para usar la funcionalidad de conversión a XYZ.")
    st.info("Instala con: pip install rdkit")
    RDKIT_AVAILABLE = False

def generar_estereoisomeros(smiles: str):
    """
    Genera todos los estereoisómeros posibles de un SMILES dado
    """
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                posiciones.append((i, True))  # ya es @@
                i += 2
            else:
                posiciones.append((i, False))  # es @ simple
                i += 1
        else:
            i += 1
    
    n = len(posiciones)
    
    # Verificación: aceptar 1, 2 o 3; rechazar > 3
    if n == 0:
        st.warning("⚠️ El SMILES no tiene centros quirales. No se generarán isómeros.")
        return [], n
    elif n > 3:
        st.error("❌ El SMILES tiene más de 3 centros quirales. No se generarán isómeros.")
        return [], n
    
    # Generar todas las combinaciones posibles
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
    
    return resultados, n

def smiles_to_xyz(smiles, mol_id):
    """
    Convierte un SMILES a formato XYZ y retorna el contenido como string
    """
    if not RDKIT_AVAILABLE:
        return None, "❌ RDKit no está disponible"
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, f"❌ Error: SMILES inválido {smiles}"
        
        # Agregar hidrógenos
        mol = Chem.AddHs(mol)
        
        # Generar geometría inicial con ETKDG
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  # reproducible
        
        embed_result = AllChem.EmbedMolecule(mol, params)
        if embed_result != 0:
            # Intentar con parámetros más flexibles
            params.useRandomCoords = True
            embed_result = AllChem.EmbedMolecule(mol, params)
            if embed_result != 0:
                return None, f"⚠️ No se pudo generar conformación 3D para {smiles}"
        
        # Optimizar con MMFF94 (si está disponible)
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception as opt_error:
            # Si falla la optimización, continuar con la estructura no optimizada
            pass
        
        # Crear contenido XYZ
        conf = mol.GetConformer()
        xyz_content = f"{mol.GetNumAtoms()}\n{smiles}\n"
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            xyz_content += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
        
        return xyz_content, f"✅ Molécula {mol_id} procesada correctamente"
        
    except Exception as e:
        return None, f"❌ Error procesando {smiles}: {str(e)}"

def crear_archivo_zip(archivos_xyz):
    """
    Crea un archivo ZIP con los archivos XYZ
    """
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for filename, content in archivos_xyz.items():
            zip_file.writestr(filename, content)
    
    return zip_buffer.getvalue()

def main():
    st.set_page_config(
        page_title="Generador de Estereoisómeros",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 Generador de Estereoisómeros")
    st.markdown("**Genera todos los estereoisómeros posibles y convierte a formato XYZ**")
    
    # Sidebar con información
    st.sidebar.title("ℹ️ Información")
    st.sidebar.markdown("""
    **Instrucciones:**
    1. Ingresa un código SMILES con centros quirales (@)
    2. El sistema acepta máximo 3 centros quirales
    3. Genera automáticamente todos los estereoisómeros
    4. Opcionalmente convierte a formato XYZ para visualización 3D
    
    **Ejemplo de SMILES:**
    - `C[C@H](O)[C@@H](N)C`
    - `N[C@@H](C)C(=O)O`
    """)
    
    # Input del usuario
    st.subheader("📝 Entrada de Datos")
    smiles_input = st.text_input(
        "👉 Ingresa el código SMILES:",
        placeholder="Ejemplo: C[C@H](O)[C@@H](N)C",
        help="Ingresa un código SMILES que contenga centros quirales marcados con @ o @@"
    )
    
    if smiles_input:
        # Generar estereoisómeros
        with st.spinner("🔄 Generando estereoisómeros..."):
            isomeros, n_centros = generar_estereoisomeros(smiles_input)
        
        if isomeros:
            st.success(f"🔎 Se encontraron {n_centros} centros quirales (@)")
            st.success(f"✅ Total estereoisómeros generados: {len(isomeros)}")
            
            # Mostrar isómeros en columnas
            st.subheader("🔬 Estereoisómeros Generados")
            
            # Crear tabs para mejor organización
            tab1, tab2, tab3 = st.tabs(["📋 Lista Completa", "💾 Descargar SMI", "🧪 Convertir a XYZ"])
            
            with tab1:
                # Mostrar todos los isómeros
                col1, col2 = st.columns(2)
                for i, isomero in enumerate(isomeros):
                    if i % 2 == 0:
                        col1.code(f"{i+1}. {isomero}")
                    else:
                        col2.code(f"{i+1}. {isomero}")
            
            with tab2:
                st.markdown("**📁 Descargar archivo SMI**")
                # Preparar contenido del archivo SMI
                smi_content = "\n".join(isomeros)
                
                st.download_button(
                    label="📥 Descargar archivo.smi",
                    data=smi_content,
                    file_name="estereoisomeros.smi",
                    mime="text/plain",
                    help="Descarga todos los estereoisómeros en formato SMI"
                )
                
                # Preview del contenido
                with st.expander("👀 Vista previa del archivo SMI"):
                    st.text(smi_content)
            
            with tab3:
                st.markdown("**🧪 Conversión a formato XYZ**")
                
                if not RDKIT_AVAILABLE:
                    st.error("❌ RDKit no está instalado. No se puede convertir a XYZ.")
                    st.info("Para instalar RDKit, usa: `pip install rdkit` o `conda install -c conda-forge rdkit`")
                else:
                    st.info("⚠️ La conversión a XYZ puede tardar unos segundos por molécula")
                    
                    if st.button("🚀 Convertir todos a XYZ", type="primary"):
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                        
                        archivos_xyz = {}
                        mensajes_log = []
                        
                        for i, smiles in enumerate(isomeros):
                            progress = (i + 1) / len(isomeros)
                            progress_bar.progress(progress)
                            status_text.text(f"Procesando molécula {i+1}/{len(isomeros)}: {smiles}")
                            
                            xyz_content, mensaje = smiles_to_xyz(smiles, i+1)
                            mensajes_log.append(mensaje)
                            
                            if xyz_content:
                                archivos_xyz[f"mol_{i+1}.xyz"] = xyz_content
                        
                        progress_bar.progress(1.0)
                        status_text.text("✅ Proceso completado!")
                        
                        # Mostrar log de procesamiento
                        with st.expander("📋 Log de procesamiento"):
                            for mensaje in mensajes_log:
                                if "❌" in mensaje or "⚠️" in mensaje:
                                    st.error(mensaje)
                                else:
                                    st.success(mensaje)
                        
                        if archivos_xyz:
                            # Crear archivo ZIP
                            zip_data = crear_archivo_zip(archivos_xyz)
                            
                            st.success(f"✅ {len(archivos_xyz)} archivos XYZ generados correctamente")
                            
                            st.download_button(
                                label="📦 Descargar archivos XYZ (ZIP)",
                                data=zip_data,
                                file_name="estereoisomeros_xyz.zip",
                                mime="application/zip",
                                help="Descarga todos los archivos XYZ comprimidos en un ZIP"
                            )
                            
                            # Mostrar preview de un archivo XYZ
                            if len(archivos_xyz) > 0:
                                with st.expander("👀 Vista previa del primer archivo XYZ"):
                                    primer_archivo = list(archivos_xyz.values())[0]
                                    st.code(primer_archivo, language=None)
                        else:
                            st.error("❌ No se pudieron generar archivos XYZ")
    
    # Footer
    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center'>
            <small>🧬 Generador de Estereoisómeros | Desarrollado con Streamlit y RDKit</small>
        </div>
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()
