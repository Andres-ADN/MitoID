# ==============================================================================
# BLOQUE 6: PREPARACIÓN DE DATOS Y SALIDA A CONSOLA (SIN PDF DIRECTO)
# ==============================================================================
# Descripción: Esta función ahora se enfoca en preparar los datos para el informe
# y generar la salida a consola. La generación de PDF se manejará externamente
# a partir de un HTML.
# ------------------------------------------------------------------------------

import pandas as pd
from . import constants
from .annotation_and_hotspots import es_variante_en_hotspot

def generar_datos_para_informe_y_consola( 
    variantes_crudas_con_locus: list,
    variantes_hgvs_normalizadas: list,
    variantes_mitomaster: list,
    rcrs_seq_str: str, 
    aligned_ref_full: str, 
    aligned_query_full: str,
    alignment_offset_ref_0based: int,
    alignment_offset_query_0based: int
) -> pd.DataFrame | None:
    print("\n--- Iniciando Paso 6: Preparación de Datos para Informe y Salida a Consola ---")
    
    max_len_variantes = 0
    if variantes_crudas_con_locus: max_len_variantes = len(variantes_crudas_con_locus)
    
    hgvs_norm_list_actualizada = list(variantes_hgvs_normalizadas)
    mitomaster_list_actualizada = list(variantes_mitomaster)

    while len(hgvs_norm_list_actualizada) < max_len_variantes: hgvs_norm_list_actualizada.append("N/A (Faltante)")
    while len(mitomaster_list_actualizada) < max_len_variantes: mitomaster_list_actualizada.append("N/A (Faltante)")

    data_for_detailed_df = []
    column_names_for_detailed_df = [
        "Posición (rCRS)", "Ref (rCRS)", "Query", "Tipo (Mutación)",
        "Región Mitocondrial", "HGVS Normalizado",
        "Formato Mitomaster", "Alineamiento" # La columna Alineamiento será un texto preformateado
    ]

    for i in range(max_len_variantes):
        var_cruda = variantes_crudas_con_locus[i] if i < len(variantes_crudas_con_locus) else {}
        hgvs_norm_str = hgvs_norm_list_actualizada[i]
        mito_fmt_str = mitomaster_list_actualizada[i] 
        
        locus_display = var_cruda.get('locus', 'N/A')
        pos_cruda_val = var_cruda.get('pos') 
        tipo_cruda_display = var_cruda.get('type', 'N/A')
        alt_cruda = var_cruda.get('alt', 'N/A')
        ref_cruda = var_cruda.get('ref', 'N/A')

        pos_display_str = str(pos_cruda_val) if pos_cruda_val is not None else "N/A"
        ref_display_str = ref_cruda if ref_cruda is not None and ref_cruda != '-' else "-"
        alt_display_str = alt_cruda if alt_cruda is not None and alt_cruda != '-' else "-"
        
        longitud_evento_crudo = 1
        if tipo_cruda_display == 'insertion':
            if isinstance(pos_cruda_val, int):
                 pos_1based_anterior = pos_cruda_val + 1 
                 pos_display_str = f"{pos_1based_anterior}_{pos_1based_anterior+1}(ins)"
                 if pos_cruda_val >= 0 and pos_cruda_val < len(rcrs_seq_str): 
                     ref_display_str = rcrs_seq_str[pos_cruda_val]
                 else: 
                     ref_display_str = "-" 
            else: 
                pos_display_str = "N/A (ins)"
                ref_display_str = "-"
            longitud_evento_crudo = len(alt_cruda) if alt_cruda != '-' else 0
        
        elif tipo_cruda_display == 'deletion':
            longitud_evento_crudo = len(ref_cruda) if ref_cruda != '-' else 0

        if mito_fmt_str == "C309CCC":
            pos_display_str = "309"
            if 0 <= 308 < len(rcrs_seq_str): ref_display_str = rcrs_seq_str[308]
            else: ref_display_str = "?"
            alt_display_str = "CCC"
        elif mito_fmt_str == "C16193CC":
            pos_display_str = "16193"
            if 0 <= 16192 < len(rcrs_seq_str): ref_display_str = rcrs_seq_str[16192]
            else: ref_display_str = "?"
            alt_display_str = "CC"
        
        mito_fmt_display_con_hotspot = mito_fmt_str 
        pos_para_hotspot_check = -1
        if tipo_cruda_display == 'insertion' and isinstance(pos_cruda_val, int): 
            pos_para_hotspot_check = pos_cruda_val 
        elif tipo_cruda_display == 'deletion' and isinstance(pos_cruda_val, int): 
            pos_para_hotspot_check = pos_cruda_val
        
        if pos_para_hotspot_check != -1 and \
           es_variante_en_hotspot(pos_para_hotspot_check, tipo_cruda_display, longitud_evento_crudo):
            mito_fmt_display_con_hotspot += "*"
        
        contexto_alineamiento_str = "N/A" # Placeholder para el texto de alineamiento
        align_idx_crudo_en_als = var_cruda.get('align_idx')

        if align_idx_crudo_en_als is not None and aligned_ref_full and aligned_query_full:
            # (La lógica para generar contexto_alineamiento_str se mantiene igual)
            align_idx_crudo_en_als = int(align_idx_crudo_en_als)
            len_evento_en_alineamiento = 1 
            if tipo_cruda_display == 'insertion': len_evento_en_alineamiento = len(var_cruda.get('alt', '-'))
            elif tipo_cruda_display == 'deletion': len_evento_en_alineamiento = len(var_cruda.get('ref', '-'))
            len_evento_en_alineamiento = max(1, len_evento_en_alineamiento)
            snip_start_in_als = max(0, align_idx_crudo_en_als - constants.ALIGNMENT_CONTEXT_WINDOW_PDF)
            snip_end_in_als = min(len(aligned_ref_full), align_idx_crudo_en_als + len_evento_en_alineamiento + constants.ALIGNMENT_CONTEXT_WINDOW_PDF)
            ref_snip_als = aligned_ref_full[snip_start_in_als:snip_end_in_als]
            query_snip_als = aligned_query_full[snip_start_in_als:snip_end_in_als]
            match_snip_chars_list = []
            for k_match in range(len(ref_snip_als)):
                if ref_snip_als[k_match] == query_snip_als[k_match] and ref_snip_als[k_match] != '-': match_snip_chars_list.append("|")
                elif ref_snip_als[k_match] != '-' and query_snip_als[k_match] != '-': match_snip_chars_list.append(".")
                else: match_snip_chars_list.append(" ")
            match_snip_str = "".join(match_snip_chars_list)
            ref_bases_antes_snip_en_als = len(aligned_ref_full[:snip_start_in_als].replace('-', ''))
            query_bases_antes_snip_en_als = len(aligned_query_full[:snip_start_in_als].replace('-', ''))
            start_pos_ref_snip_1based = alignment_offset_ref_0based + ref_bases_antes_snip_en_als + 1
            end_pos_ref_snip_1based = start_pos_ref_snip_1based + len(ref_snip_als.replace('-', '')) - 1
            start_pos_query_snip_1based = alignment_offset_query_0based + query_bases_antes_snip_en_als + 1
            end_pos_query_snip_1based = start_pos_query_snip_1based + len(query_snip_als.replace('-', '')) - 1
            if end_pos_ref_snip_1based < start_pos_ref_snip_1based: end_pos_ref_snip_1based = start_pos_ref_snip_1based 
            if end_pos_query_snip_1based < start_pos_query_snip_1based: end_pos_query_snip_1based = start_pos_query_snip_1based
            s_ref_pos_start_str = str(start_pos_ref_snip_1based).rjust(constants.ALIGNMENT_POS_NUM_WIDTH_PDF)
            s_ref_pos_end_str = str(end_pos_ref_snip_1based).ljust(constants.ALIGNMENT_POS_NUM_WIDTH_PDF)
            s_query_pos_start_str = str(start_pos_query_snip_1based).rjust(constants.ALIGNMENT_POS_NUM_WIDTH_PDF)
            s_query_pos_end_str = str(end_pos_query_snip_1based).ljust(constants.ALIGNMENT_POS_NUM_WIDTH_PDF)
            padding_entre_pos_y_seq = "  "
            line1_ref_text = f"{'Ref:'.ljust(constants.ALIGNMENT_LABEL_WIDTH_PDF)}{s_ref_pos_start_str}{padding_entre_pos_y_seq}{ref_snip_als}{padding_entre_pos_y_seq}{s_ref_pos_end_str}"
            match_line_prefix_padding = ' ' * (constants.ALIGNMENT_POS_NUM_WIDTH_PDF + len(padding_entre_pos_y_seq))
            line2_match_text = f"{'Match:'.ljust(constants.ALIGNMENT_LABEL_WIDTH_PDF)}{match_line_prefix_padding}{match_snip_str}" 
            line3_query_text = f"{'Query:'.ljust(constants.ALIGNMENT_LABEL_WIDTH_PDF)}{s_query_pos_start_str}{padding_entre_pos_y_seq}{query_snip_als}{padding_entre_pos_y_seq}{s_query_pos_end_str}"
            # Para HTML, podríamos querer las líneas separadas o usar <pre>
            contexto_alineamiento_str = f"{line1_ref_text}\n{line2_match_text}\n{line3_query_text}"


        data_for_detailed_df.append({
            "Posición (rCRS)": pos_display_str, 
            "Ref (rCRS)": ref_display_str,
            "Query": alt_display_str,
            "Tipo (Mutación)": tipo_cruda_display,
            "Región Mitocondrial": locus_display,
            "HGVS Normalizado": hgvs_norm_str,
            "Formato Mitomaster": mito_fmt_display_con_hotspot,
            "Alineamiento": contexto_alineamiento_str # Se guardará como texto multilínea
        })
    
    df_detallado = pd.DataFrame(data_for_detailed_df, columns=column_names_for_detailed_df)
    
    # --- Salida a consola  ---
    if not df_detallado.empty:
        df_consola_resumen = df_detallado[[
            "Posición (rCRS)", "Ref (rCRS)", "Query", "Tipo (Mutación)", 
            "Formato Mitomaster", "Región Mitocondrial"
        ]].copy()
        print("\nTabla de variantes (resumen para consola con locus y Mitomaster):\n")
        try:
            print(df_consola_resumen.to_markdown(index=False))
        except ImportError: 
            print(df_consola_resumen.to_string(index=False))
    else:
        print("\nNo se encontraron variantes para mostrar en la tabla de resumen de consola.")
        
    print("\nPreparación de datos para informe completada.")
    return df_detallado