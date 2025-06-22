# ==============================================================================
# BLOQUE 3: UTILIDADES DE ALINEAMIENTO Y EXTRACCIÓN DE VARIANTES (V4)
# ==============================================================================
# Descripción: Funciones para realizar el alineamiento entre secuencias
# (referencia y consulta) y para extraer las diferencias (variantes crudas)
# a partir del alineamiento resultante.
# SE HA ELIMINADO CÓDIGO REDUNDANTE/MUERTO AL FINAL DE extraer_variantes_crudas.
# ------------------------------------------------------------------------------
from Bio import Align
from Bio.Seq import Seq # Necesario para Seq(rcrs_seq) en main
import traceback
from . import constants # Para POS_3107_BLACKLISTED y otras

def realizar_alineamiento(rcrs_seq: Seq, query_seq: Seq, modo_alineamiento: str = 'local') -> Align.Alignment | None:
    print(f"\\n--- Iniciando Paso 3: Alineamiento en modo '{modo_alineamiento}' ---")
    aligner = Align.PairwiseAligner()
    aligner.mode = modo_alineamiento

    # Matriz de sustitución personalizada (más adecuada para secuencias cercanas o ADN)
    # Puntuaciones: Match=3, Mismatch=-3, N vs Base=1, N vs N=1
    matrix_data = {
        ('A','A'):3, ('A','T'):-3, ('A','G'):-3, ('A','C'):-3, ('A','N'):1,
        ('T','T'):3, ('T','A'):-3, ('T','G'):-3, ('T','C'):-3, ('T','N'):1,
        ('G','G'):3, ('G','A'):-3, ('G','T'):-3, ('G','C'):-3, ('G','N'):1,
        ('C','C'):3, ('C','A'):-3, ('C','T'):-3, ('C','G'):-3, ('C','N'):1,
        ('N','N'):1
    }
    # Asegurar simetría en la matriz
    for k_tuple in list(matrix_data.keys()):
        if (k_tuple[1], k_tuple[0]) not in matrix_data:
            matrix_data[(k_tuple[1], k_tuple[0])] = matrix_data[k_tuple]
    
    custom_matrix = Align.substitution_matrices.Array(data=matrix_data)
    aligner.substitution_matrix = custom_matrix
    
    # Penalizaciones por gap (ajustar según necesidad, estos son valores comunes)
    aligner.open_gap_score = -7.0
    aligner.extend_gap_score = -2.0
    
    print(f"Configuración del alineador: Modo={aligner.mode}, OpenGap={aligner.open_gap_score}, ExtendGap={aligner.extend_gap_score}")

    try:
        alignments = list(aligner.align(rcrs_seq, query_seq)) # Convertir generador a lista
        if alignments:
            best_alignment = alignments[0] # Biopython devuelve el mejor (o uno de los mejores) primero
            print(f"Alineamiento completado. Puntuación: {best_alignment.score}")
            return best_alignment
        else:
            print("No se encontraron alineamientos.")
            return None
    except Exception as e:
        print(f"Ocurrió un error durante el alineamiento: {e}")
        import traceback
        traceback.print_exc()
        return None

def get_substitution_type(ref_base: str, alt_base: str) -> str:

    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    ref_b, alt_b = ref_base.upper(), alt_base.upper()

    if not (ref_b in "ACGTN" and alt_b in "ACGTN"): return "desconocido"
    if ref_b == 'N' or alt_b == 'N': return "substitution (con N)"
    if ref_b == alt_b: return "identico"

    if (ref_b in purines and alt_b in purines) or \
       (ref_b in pyrimidines and alt_b in pyrimidines):
        return "transition"
    elif (ref_b in purines and alt_b in pyrimidines) or \
         (ref_b in pyrimidines and alt_b in purines):
        return "transversion"
    return "substitution"

def extraer_variantes_crudas(alineamiento: Align.Alignment, rcrs_id: str, query_id: str, rcrs_len: int) -> tuple[list, str, str, int, int]:
    print("\\n--- Iniciando Paso 4: Extracción de Variantes Crudas ---")
    
    raw_variants = []
    als_ref = ""
    als_query = ""
    alignment_offset_ref = 0 
    alignment_offset_query = 0 

    if not alineamiento:
        print("Error: Objeto de alineamiento no proporcionado.")
        return raw_variants, als_ref, als_query, alignment_offset_ref, alignment_offset_query

    offset_determinad_con_exito = False
    if hasattr(alineamiento, 'coordinates') and hasattr(alineamiento.coordinates, 'shape'):
        if alineamiento.coordinates.shape[0] >= 1 and alineamiento.coordinates.shape[1] > 0:
            try:
                alignment_offset_ref = int(alineamiento.coordinates[0,0])
                offset_determinad_con_exito = True 
            except (IndexError, TypeError, ValueError) as e:
                print(f"Advertencia: No se pudo extraer offset_ref de coordinates (shape={alineamiento.coordinates.shape}): {e}")
        if alineamiento.coordinates.shape[0] >= 2 and alineamiento.coordinates.shape[1] > 0:
            try:
                alignment_offset_query = int(alineamiento.coordinates[1,0])
            except (IndexError, TypeError, ValueError) as e:
                print(f"Advertencia: No se pudo extraer offset_query de coordinates (shape={alineamiento.coordinates.shape}): {e}")
    
    if not offset_determinad_con_exito and hasattr(alineamiento, 'aligned') and isinstance(alineamiento.aligned, tuple) and len(alineamiento.aligned) == 2:
        print("Información: Intentando obtener offsets desde 'alineamiento.aligned'.")
        if isinstance(alineamiento.aligned[0], tuple) and len(alineamiento.aligned[0]) > 0 and \
           isinstance(alineamiento.aligned[0][0], tuple) and len(alineamiento.aligned[0][0]) == 2:
            try: 
                alignment_offset_ref = alineamiento.aligned[0][0][0]
                offset_determinad_con_exito = True
            except: pass 
        if isinstance(alineamiento.aligned[1], tuple) and len(alineamiento.aligned[1]) > 0 and \
           isinstance(alineamiento.aligned[1][0], tuple) and len(alineamiento.aligned[1][0]) == 2:
            try: alignment_offset_query = alineamiento.aligned[1][0][0]
            except: pass
    
    if not offset_determinad_con_exito:
         print(f"ADVERTENCIA FINAL: No se pudieron determinar los offsets de inicio del alineamiento de forma precisa. Asumiendo offset_ref=0 y offset_query=0.")
         alignment_offset_ref = 0 
         alignment_offset_query = 0

    ref_pos_in_align_segment = 0 
    try:
        als_ref = str(alineamiento[0])
        als_query = str(alineamiento[1])
        if not als_ref or not als_query or len(als_ref) != len(als_query):
            print(f"Error en obtención de secuencias alineadas: Longitudes no coinciden. Ref:{len(als_ref)}, Query:{len(als_query)}")
            return [], "", "", 0, 0 
        
        i = 0
        while i < len(als_ref):
            r_base, q_base = als_ref[i], als_query[i]
            align_col_idx = i 
            consumes_ref_current_char = (r_base != '-')
            current_abs_pos_1_based_ref = alignment_offset_ref + ref_pos_in_align_segment + 1
            
            if r_base != q_base:
                if r_base != '-' and q_base != '-': 
                    var_type = get_substitution_type(r_base, q_base)
                    raw_variants.append({'pos': current_abs_pos_1_based_ref, 'ref': r_base.upper(), 'alt': q_base.upper(), 'type': var_type, 'align_idx': align_col_idx})
                    i += 1
                elif r_base == '-': 
                    ins_seq = ""
                    start_align_idx = align_col_idx
                    insertion_point_ref_pos_0_based = alignment_offset_ref + ref_pos_in_align_segment -1
                    if ref_pos_in_align_segment == 0 and alignment_offset_ref == 0: 
                        insertion_point_ref_pos_0_based = -1 
                    
                    while i < len(als_ref) and als_ref[i] == '-':
                        ins_seq += als_query[i]
                        i += 1
                    raw_variants.append({'pos': insertion_point_ref_pos_0_based, 'ref': '-', 'alt': ins_seq.upper(), 'type': 'insertion', 'align_idx': start_align_idx})
                elif q_base == '-': 
                    del_seq = ""
                    start_align_idx = align_col_idx
                    del_start_abs_pos_1_based = current_abs_pos_1_based_ref
                    temp_ref_pos_tracker_del_segment = ref_pos_in_align_segment 
                    while i < len(als_ref) and als_query[i] == '-':
                        del_seq += als_ref[i]
                        i += 1
                        if als_ref[i-1] != '-':
                             temp_ref_pos_tracker_del_segment +=1 
                    raw_variants.append({'pos': del_start_abs_pos_1_based, 'ref': del_seq.upper(), 'alt': '-', 'type': 'deletion', 'align_idx': start_align_idx})
                    ref_pos_in_align_segment = temp_ref_pos_tracker_del_segment 
                    continue 
            else: 
                i += 1
            
            if consumes_ref_current_char: 
                ref_pos_in_align_segment += 1
        
        print(f"Se encontraron {len(raw_variants)} diferencias crudas (antes de filtrar).")

        # --- INICIO DE LA LÓGICA DE FILTRADO MODIFICADA PARA ARTEFACTOS Y 3107 ---
    
        variants_final_filtradas = [] 
        num_artefactos_filtrados_total = 0
        
        for var_cruda in raw_variants:
            filtrar_esta_variante = False
            
            tipo_var = var_cruda.get('type', '')
            pos_var = var_cruda.get('pos')
            ref_var = var_cruda.get('ref', '').upper()
            alt_var = var_cruda.get('alt', '').upper()

            # Condición 0: Filtrar la posición 3107 según SWGDAM/EMPOP
            if pos_var == constants.POS_3107_BLACKLISTED:
                filtrar_esta_variante = True
                # print(f"DEBUG: Filtrando por Condición 0 (Posición {POS_3107_BLACKLISTED}): {var_cruda}")

            # Condición 1: Artefacto m.3107delN (aunque ya cubierto por Condición 0, lo mantenemos por si acaso)
            elif tipo_var == 'deletion' and pos_var == 3107 and ref_var == 'N':
                filtrar_esta_variante = True
                # print(f"DEBUG: Filtrando por Condición 1 (3107delN): {var_cruda}")

            # Condición 2: Artefacto C3106N (Sustitución C>N en pos 3106)
            # (La rCRS tiene C en 3106. Si la query tiene N allí)
            elif tipo_var.startswith('substitution') and pos_var == 3106 and ref_var == 'C' and alt_var == 'N':
                filtrar_esta_variante = True
                # print(f"DEBUG: Filtrando por Condición 2 (C3106N): {var_cruda}")

            # Condición 3: Artefacto N3107X (Sustitución N>[ACGT] en pos 3107)
            # (La rCRS tiene N en 3107. Si la query tiene una base canónica allí)
            elif tipo_var.startswith('substitution') and pos_var == 3107 and ref_var == 'N' and alt_var in ['A', 'C', 'G', 'T']:
                filtrar_esta_variante = True
                # print(f"DEBUG: Filtrando por Condición 3 (N3107X): {var_cruda}")

            if filtrar_esta_variante:
                num_artefactos_filtrados_total += 1
            else:
                variants_final_filtradas.append(var_cruda)
        
        if num_artefactos_filtrados_total > 0:
            print(f"Se filtraron {num_artefactos_filtrados_total} variantes problemáticas/artefactos (ej. 3107, C3106N, N3107X). Variantes restantes: {len(variants_final_filtradas)}")
        # --- FIN DE LA LÓGICA DE FILTRADO MODIFICADA ---
        
        return variants_final_filtradas, als_ref, als_query, alignment_offset_ref, alignment_offset_query
    
    except Exception as e: 
        print(f"Error extrayendo variantes: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return [], "", "", alignment_offset_ref, alignment_offset_query