# ==============================================================================
# BLOQUE 9: (main)
# ==============================================================================
# Descripción: La función `main` coordina todos los pasos del análisis.
# Versión con inicialización robusta de variables y control de flujo.
# ------------------------------------------------------------------------------

import os
from Bio.Seq import Seq # Necesario para Seq()

# Importaciones de tus módulos personalizados
from . import constants
from .feature_extraction import cargar_secuencia_fasta, cargar_y_extraer_features_rcrs
from .alignment_and_variant_calling import realizar_alineamiento, extraer_variantes_crudas
from .annotation_and_hotspots import anotar_locus_variante, obtener_hvs_region
from .hgvs_and_nomenclature import normalizar_y_nombrar_hgvs, formatear_estilo_mitomaster, formatear_variantes_empop
from .report_data_preparation import generar_datos_para_informe_y_consola
from .report_generation import generar_informe_html, convertir_html_a_pdf
from .track_viewer import crear_track_viewer_interactivo 


def main(query_fasta_file_path: str):
    print("Iniciando análisis de secuencia mitocondrial...")
    os.makedirs("secuencias", exist_ok=True)
    
    # --- 1. Gestión de Archivos Dummy de Referencia (si no existen) ---
    if not os.path.exists(constants.RCRS_FASTA_PATH):
        print(f"ADVERTENCIA: Archivo rCRS FASTA no encontrado. Creando dummy: {constants.RCRS_FASTA_PATH}")
        dummy_rcrs_seq_content = "GATTACA" * (16569 // 7) + "GATTACA"[:16569 % 7]
        with open(constants.RCRS_FASTA_PATH, "w") as f:
            f.write(f">NC_012920.1 rCRS_dummy_sequence\\n{dummy_rcrs_seq_content}\\n")
    
    if not os.path.exists(constants.RCRS_GENBANK_PATH):
        print(f"ADVERTENCIA: Archivo rCRS GenBank no encontrado. Creando dummy: {constants.RCRS_GENBANK_PATH}")
        dummy_gb_content = f"""LOCUS       NC_012920.1            16569 bp    DNA     circular CON 28-MAY-2025
DEFINITION  Dummy Homo sapiens mitochondrion, rCRS, for testing.
FEATURES             Location/Qualifiers
     source          1..16569
     D-loop          complement(join(1..576,16024..16569))
     gene            648..1601
                     /gene="MT-RNR1"
                     /product="12S ribosomal RNA"
     gene            1671..3229
                     /gene="MT-RNR2"
                     /product="16S ribosomal RNA"
     tRNA            577..647
                     /gene="MT-TF"
                     /product="tRNA-Phe"
     tRNA            1602..1670
                     /gene="MT-TV"
                     /product="tRNA-Val"
     CDS             3307..4262
                     /gene="MT-ND1"
                     /product="NADH dehydrogenase subunit 1"
     CDS             5904..7445
                     /gene="MT-CO1"
                     /product="cytochrome c oxidase subunit I"
     tRNA            3230..3304
                     /gene="MT-TL1"
                     /product="tRNA-Leu(UUR)"
ORIGIN      
        1 {'a'*1000}
        {'g'*1000}
        {'t'*1000}
        {'c'*1000}
        {'a'*1000}
        {'g'*1000}
        {'t'*1000}
        {'c'*1000}
        {'a'*1000}
        {'g'*1000}
        {'t'*1000}
        {'c'*1000}
        {'a'*1000}
        {'g'*1000}
        {'t'*1000}
        {'c'*569}
//"""
        with open(constants.RCRS_GENBANK_PATH, "w") as f:
            f.write(dummy_gb_content)
    
    # --- 2. Carga de Secuencias y Features de Referencia ---
    rcrs_fasta_record = cargar_secuencia_fasta(constants.RCRS_FASTA_PATH, id_esperado="NC_012920.1")
    query_record = cargar_secuencia_fasta(query_fasta_file_path)
    features_rcrs_gb = cargar_y_extraer_features_rcrs(constants.RCRS_GENBANK_PATH)

    if not rcrs_fasta_record or not query_record:
        print("Error crítico: No se pudieron cargar las secuencias FASTA. Abortando.")
        return
    
    query_id_original = query_record.id
    query_id_sanitized = "".join(c if c.isalnum() or c in ('_','-','.') else '_' for c in query_id_original)
    
    # --- 3. Definición de Rutas de Salida ---
    output_html_temp_path = f"temp_informe_variantes_{query_id_sanitized}.html"
    output_pdf_final_path = f"Informe_Variantes_{query_id_sanitized}_vs_rCRS.pdf"
    output_html_track_viewer_path = f"Track_Viewer_{query_id_sanitized}_vs_rCRS.html"
    
    print(f"El archivo PDF de salida se guardará como: {output_pdf_final_path}")
    print(f"El informe HTML detallado se guardará como: {output_html_temp_path}")
    print(f"El Track Viewer interactivo se guardará como: {output_html_track_viewer_path}")

    # --- 4. Preparación y Alineamiento de Secuencias ---
    rcrs_sequence_str = str(rcrs_fasta_record.seq).upper()
    query_sequence_str = str(query_record.seq).upper()
    
    modo_de_alineamiento = 'global' if len(query_sequence_str) >= len(rcrs_sequence_str) * 0.8 else 'local'
    print(f"Usando modo de alineamiento: {modo_de_alineamiento}")
    
    # Inicialización de variables de alineamiento y extracción de variantes crudas
    # con valores por defecto seguros.
    mejor_alineamiento = None
    variantes_crudas_lista = []
    als_ref_str = ""
    als_query_str = ""
    offset_ref_0b = 0
    offset_query_0b = 0
    
# --- 5. Inicialización de TODAS las variables de resultado de variantes y formateo ---
    variantes_crudas_con_locus_lista = []
    variantes_hgvs_norm_lista = []
    variantes_mitomaster_formato_lista = [] # <-- ¡Esta es la que faltaba!
    empop_query_final_str = ""
    num_empop_variantes_final = 0
    empop_variantes_list_for_tv = [] # Lista separada de strings EMPOP para el Track Viewer

    mejor_alineamiento = realizar_alineamiento(Seq(rcrs_sequence_str), Seq(query_sequence_str), modo_alineamiento=modo_de_alineamiento)
    
    if mejor_alineamiento: # Si el alineamiento fue exitoso, intentar extraer variantes
        variantes_crudas_lista, als_ref_str, als_query_str, offset_ref_0b, offset_query_0b = extraer_variantes_crudas(
            mejor_alineamiento, rcrs_fasta_record.id, query_id_original, len(rcrs_sequence_str)
        )
    else: # Si el alineamiento falló, se emite un aviso y se procede con listas vacías.
        print("El alineamiento falló o no se encontraron alineamientos. Las listas de variantes estarán vacías.")
        # Las variables ya están inicializadas a vacío/cero arriba, así que no se necesita re-inicializar aquí.
    
    # --- 5. Inicialización de TODAS las variables de resultado de variantes y formateo ---
    # ¡ESTA ES LA SECCIÓN MÁS IMPORTANTE PARA EVITAR NameError!
    # Asegura que todas estas listas y variables existan, incluso si no hay variantes procesadas.
    variantes_crudas_con_locus_lista = []
    variantes_hgvs_norm_lista = []
    variantes_mitomaster_formato_lista = []
    empop_query_final_str = ""
    num_empop_variantes_final = 0
    empop_variantes_list_for_tv = [] # Lista separada de strings EMPOP para el Track Viewer
    
    if variantes_crudas_lista: # Solo procesar variantes si la lista de variantes crudas NO está vacía
        # --- 5.1. Anotación de Locus y HVS para Variantes Crudas ---
        for var_cruda_dict in variantes_crudas_lista:
            pos_var_cruda = var_cruda_dict.get('pos')
            tipo_var_cruda = var_cruda_dict.get('type', '')
            longitud_evento = 0
            pos_1based_para_locus = -1

            if tipo_var_cruda == 'insertion':
                if isinstance(pos_var_cruda, int): pos_1based_para_locus = pos_var_cruda + 1
                longitud_evento = len(var_cruda_dict.get('alt', ''))
            elif tipo_var_cruda == 'deletion':
                if isinstance(pos_var_cruda, int): pos_1based_para_locus = pos_var_cruda
                longitud_evento = len(var_cruda_dict.get('ref', ''))
            else: # Sustitución
                if isinstance(pos_var_cruda, int): pos_1based_para_locus = pos_var_cruda
                longitud_evento = 1

            locus_anotado = "N/A"
            hvs_region_anotada = None
            # Asegúrate de que features_rcrs_gb no esté vacío antes de intentar anotar
            if features_rcrs_gb and pos_1based_para_locus > 0:
                locus_anotado = anotar_locus_variante(pos_1based_para_locus, tipo_var_cruda, longitud_evento, features_rcrs_gb)
                hvs_region_anotada = obtener_hvs_region(pos_1based_para_locus)

            var_cruda_copia = var_cruda_dict.copy()
            var_cruda_copia['locus'] = locus_anotado
            var_cruda_copia['hvs_region'] = hvs_region_anotada
            variantes_crudas_con_locus_lista.append(var_cruda_copia)
            
        # --- 5.2. Normalización y Formateo de Variantes (HGVS, Mitomaster, EMPOP) ---
        ref_acc_hgvs = rcrs_fasta_record.id
        variantes_hgvs_norm_lista = normalizar_y_nombrar_hgvs(variantes_crudas_con_locus_lista, ref_acc_hgvs, rcrs_sequence_str)
        variantes_mitomaster_formato_lista = formatear_estilo_mitomaster(variantes_hgvs_norm_lista, ref_acc_hgvs, rcrs_sequence_str, variantes_crudas_lista)
        
        # formatear_variantes_empop devuelve 3 valores (cadena, número, lista de strings)
        empop_query_final_str, num_empop_variantes_final, empop_variantes_list_for_tv = formatear_variantes_empop(variantes_crudas_con_locus_lista, rcrs_sequence_str)
    
    # --- 6. Resumen de Variantes para Consola (Mitomaster y EMPOP) ---
    # Estas listas ya han sido inicializadas y rellenadas (o quedaron vacías) en el paso 5.
    vm_validas_html = [vm for vm in variantes_mitomaster_formato_lista if not any(err in vm for err in ["ERROR", "FALLIDO", "FORMATO_", "UNPROCESSED", "N/A"])]
    num_vm_validas_html = len(vm_validas_html)
    str_vm_validas_html = ", ".join(vm_validas_html) if num_vm_validas_html > 0 else "N/A"
    
    print("\nLista final de variantes (estilo Mitomaster):")
    if vm_validas_html: print(str_vm_validas_html)
    else: print("No hay variantes formateadas Mitomaster válidas para mostrar.")
    
    print(f"\nCadena de Consulta EMPOP ({num_empop_variantes_final} variantes):\n")
    print(empop_query_final_str if empop_query_final_str else "N/A")
    
    # --- 7. Generación de Informe Tabular (HTML y PDF) ---
    # Todos los argumentos pasados aquí ya han sido inicializados y/o rellenados.
    df_detallado_final = generar_datos_para_informe_y_consola(
        variantes_crudas_con_locus_lista,
        variantes_hgvs_norm_lista,
        variantes_mitomaster_formato_lista,
        rcrs_sequence_str,
        als_ref_str,
        als_query_str,
        offset_ref_0b,
        offset_query_0b
    )
    
    if df_detallado_final is not None:
        ruta_html_generado = generar_informe_html(
            df_detallado=df_detallado_final,
            query_id=query_id_original,
            num_variantes_mitomaster_validas=num_vm_validas_html,
            str_variantes_mitomaster_validas=str_vm_validas_html,
            num_variantes_empop=num_empop_variantes_final,
            str_empop_query=empop_query_final_str,
            output_html_path=output_html_temp_path
        )
        
        if ruta_html_generado:
            convertir_html_a_pdf(ruta_html_generado, output_pdf_final_path)
            # Opcional: eliminar el archivo HTML temporal después de la conversión
            # try:
            #     os.remove(ruta_html_generado)
            #     print(f"Archivo HTML temporal '{ruta_html_generado}' eliminado.")
            # except OSError as e_del:
            #     print(f"Error al eliminar el archivo HTML temporal '{ruta_html_generado}': {e_del}")
    # --- 8. Generación del Track Viewer Interactivo ---
    print("\n--- Iniciando generación de Track Viewer interactivo ---")
    crear_track_viewer_interactivo( # Llama directamente a la función del módulo track_viewer
        caracteristicas_rcrs=features_rcrs_gb,
        lista_variantes_crudas_con_locus=variantes_crudas_con_locus_lista,
        variantes_hgvs_normalizadas=variantes_hgvs_norm_lista,
        variantes_mitomaster_formateadas=variantes_mitomaster_formato_lista,
        variantes_empop_formateadas=empop_variantes_list_for_tv,
        rcrs_id=rcrs_fasta_record.id,
        query_id=query_record.id,
        rcrs_sequence_str=str(rcrs_fasta_record.seq).upper(),
        aligned_ref_full=als_ref_str,
        aligned_query_full=als_query_str,
        alignment_offset_ref_0based=offset_ref_0b,
        alignment_offset_query_0based=offset_query_0b,
        output_html_path_tv=output_html_track_viewer_path
    )

    print("\nAnálisis completado.\n")

if __name__ == "__main__":
    # Definir la ruta al archivo FASTA de la secuencia query para el análisis
    ruta_query_analisis = constants.RCRS_FASTA_PATH.replace("_rCRS.fasta", "_1.fasta") # Usando la constante para la ruta de la dummy query

    # Crear un archivo de consulta dummy si no existe
    if not os.path.exists(ruta_query_analisis):
        print(f"ADVERTENCIA: Archivo de consulta '{ruta_query_analisis}' no encontrado. Creando dummy.")
        os.makedirs(os.path.dirname(ruta_query_analisis), exist_ok=True)
        # Secuencia dummy con algunas diferencias respecto a una rCRS hipotética (e.g., primeros 100pb)
        dummy_query_seq_content = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTG" # ~160pb
        dummy_query_seq_content += "A" * (16569 - len(dummy_query_seq_content) -100) # Relleno
        dummy_query_seq_content += "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA" # Final diferente
        with open(ruta_query_analisis, "w") as f:
            f.write(f">Query_Dummy_Sequence\\n{dummy_query_seq_content[:16569]}\\n")

    main(query_fasta_file_path=ruta_query_analisis)

    print("\nAnálisis completado.\n")