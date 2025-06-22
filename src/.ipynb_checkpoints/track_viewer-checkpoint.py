# ==============================================================================
# BLOQUE 8: FUNCIONES DEL TRACK VIEWER INTERACTIVO
# ==============================================================================
# Descripción: Funciones para crear el Track Viewer interactivo usando Plotly.
# Recibe datos ya procesados del flujo principal.
# ------------------------------------------------------------------------------


import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from . import constants 

def crear_track_viewer_interactivo(
    caracteristicas_rcrs: list,
    lista_variantes_crudas_con_locus: list,
    variantes_hgvs_normalizadas: list,
    variantes_mitomaster_formateadas: list,
    variantes_empop_formateadas: list,
    rcrs_id: str,
    query_id: str,
    rcrs_sequence_str: str,
    aligned_ref_full: str,
    aligned_query_full: str,
    alignment_offset_ref_0based: int,
    alignment_offset_query_0based: int,
    output_html_path_tv: str = None
):
    print("\n--- Iniciando Generación de Track Viewer Interactivo ---")

    RCRS_LENGTH_EXPECTED = len(rcrs_sequence_str) # Usar la longitud real de la secuencia rCRS cargada
    
    # Configuración de pistas: Ajustar posiciones 'y' para DISTRIBUCIÓN CLARA
    pistas_config_tv = {
        "Query":  {"y": 7.0, "color": "lightcoral", "label_full": f"Secuencia Query: {query_id}", "height": 0.8},
        "gene":   {"y": 5.0, "color": "royalblue", "label_full": "Genes Codificantes de Proteínas", "height": 0.8},
        "rRNA":   {"y": 3.0, "color": "forestgreen", "label_full": "rRNA (Ribosomal RNA)", "height": 0.8},
        "tRNA":   {"y": 1.0, "color": "orchid", "label_full": "tRNA (Transfer RNA)", "height": 0.6},
        "D-loop": {"y": -1.0, "color": "gold", "label_full": "D-loop (Región de Control)", "height": 0.8}
    }
    
    # Orden de las pistas en el eje Y (de abajo a arriba para el gráfico)
    track_order = ["D-loop", "tRNA", "rRNA", "gene", "Query"] 

    fig = go.Figure()

    # --- Dibujar Línea de Referencia Base de rCRS ---
    fig.add_trace(go.Scatter(
        x=[0, RCRS_LENGTH_EXPECTED],
        y=[-2.5, -2.5], # Línea base aún más abajo
        mode='lines',
        line=dict(color='darkgrey', width=1.5),
        showlegend=False,
        hoverinfo='none',
        name=f"rCRS: {rcrs_id}"
    ))
    fig.add_annotation( # Etiqueta para rCRS en la línea base
        x=RCRS_LENGTH_EXPECTED/2,
        y=-2.5,
        text=f"Referencia: {rcrs_id}",
        showarrow=False,
        yshift=-10,
        font=dict(size=10, color="darkgrey")
    )


    # --- Dibujar Features de rCRS como barras (go.Bar) ---
    for feat_type_key in track_order:
        if feat_type_key == "Query":
            continue
        
        config = pistas_config_tv.get(feat_type_key)
        if not config:
            continue

        y_pos = config["y"]
        height = config["height"]
        color = config["color"]
        label_full = config["label_full"]

        features_of_this_type = [f for f in caracteristicas_rcrs if f["tipo"] == feat_type_key]
        
        hover_texts = []
        
        for feat in features_of_this_type:
            # Hovertext para la barra
            hover_info = f"<b>Característica:</b> {feat['nombre']}<br>" + \
                         f"<b>Tipo:</b> {feat['tipo']}<br>" + \
                         f"<b>Posición:</b> {feat['inicio']}-{feat['fin']} (1-based)<br>" + \
                         f"<b>Longitud:</b> {feat['longitud']} pb<br>"
            if feat['tipo'] in ["gene", "rRNA", "tRNA"]:
                hover_info += f"<b>Hebra:</b> {'Directa (+)' if feat['hebra'] == 1 else 'Complementaria (-)'}<extra></extra>"
            else:
                hover_info += "<extra></extra>"

            hover_texts.append(hover_info)
            
            # Dibujar flechas para hebra de genes y rRNA
            if feat["tipo"] in ["gene", "rRNA"]:
                arrow_length_pixels = min(50, feat["longitud"] * 0.1)
                if arrow_length_pixels > feat["longitud"] / 2:
                    arrow_length_pixels = feat["longitud"] / 2

                if feat["hebra"] == 1: # Directa (+)
                    fig.add_annotation(
                        x=feat["fin"], y=y_pos,
                        ax=feat["fin"] - arrow_length_pixels, ay=y_pos,
                        xref="x", yref="y", axref="x", ayref="y",
                        showarrow=True,
                        arrowhead=2,
                        arrowsize=1.5,
                        arrowwidth=1,
                        arrowcolor="black",
                        opacity=0.8
                    )
                elif feat["hebra"] == -1: # Complementaria (-)
                    fig.add_annotation(
                        x=feat["inicio"], y=y_pos,
                        ax=feat["inicio"] + arrow_length_pixels, ay=y_pos,
                        xref="x", yref="y", axref="x", ayref="y",
                        showarrow=True,
                        arrowhead=2,
                        arrowsize=1.5,
                        arrowwidth=1,
                        arrowcolor="black",
                        opacity=0.8
                    )

        # Añadir las barras para este tipo de feature
        if features_of_this_type:
            fig.add_trace(go.Bar(
                x=[f["inicio"] for f in features_of_this_type],
                y=[height] * len(features_of_this_type),
                width=[f["longitud"] for f in features_of_this_type],
                offset=[0] * len(features_of_this_type),
                base=[y_pos - height/2] * len(features_of_this_type),
                marker_color=color,
                marker_line_color='black',
                marker_line_width=0.5,
                opacity=0.7,
                name=label_full,
                hovertemplate=hover_texts,
                showlegend=True,
                orientation='v'
            ))

            # Añadir etiquetas de texto encima de las barras (como anotaciones)
            for feat in features_of_this_type:
                # Umbral de longitud para mostrar etiqueta
                min_len_for_label_gene_rrna = RCRS_LENGTH_EXPECTED * 0.005 # ~80bp para 16kb
                min_len_for_label_tRNA = 40 # 40 bases para tRNAs
                
                show_label = False
                if feat["tipo"] == "D-loop": show_label = True
                elif feat["tipo"] in ["gene", "rRNA"] and feat["longitud"] >= min_len_for_label_gene_rrna: show_label = True
                elif feat["tipo"] == "tRNA" and feat["longitud"] >= min_len_for_label_tRNA: show_label = True

                if show_label:
                    font_size = 8
                    if feat["tipo"] == "tRNA": font_size = 6
                    
                    texto_label = feat["nombre"]
                    if feat["tipo"] in ["gene", "rRNA"] and feat["hebra"] is not None and feat["hebra"] != 0 and feat["longitud"] > (RCRS_LENGTH_EXPECTED * 0.015):
                        texto_label = f"{feat['nombre']} ({'+' if feat['hebra'] == 1 else '-'})"

                    fig.add_annotation(
                        x=feat["inicio"] + feat["longitud"] / 2, y=y_pos,
                        text=texto_label,
                        showarrow=False,
                        font=dict(size=font_size, color="black"),
                        xanchor="center", yanchor="middle",
                        align="center",
                        valign="middle",
                    )


    # --- Dibujar la Pista de la Query (como barra) ---
    query_track_config = pistas_config_tv.get("Query")
    if query_track_config and query_id != 'N/A' and aligned_ref_full:
        query_start_on_ref = alignment_offset_ref_0based
        query_width_on_ref_alignment = len(aligned_ref_full.replace('-', ''))
        
        if query_width_on_ref_alignment > 0:
            fig.add_trace(go.Bar(
                x=[query_start_on_ref],
                y=[query_track_config["height"]],
                width=[query_width_on_ref_alignment],
                offset=[0],
                base=[query_track_config["y"] - query_track_config["height"]/2],
                marker_color=query_track_config["color"],
                marker_line_color='black',
                marker_line_width=0.5,
                opacity=0.7,
                name=query_track_config["label_full"],
                hovertemplate=f"<b>Muestra:</b> {query_id}<br>" +
                              f"<b>Alineada a rCRS (1-based):</b> {query_start_on_ref+1}-{query_start_on_ref + query_width_on_ref_alignment} <br>" +
                              f"<b>Longitud de Query alineada:</b> {query_width_on_ref_alignment} bp<extra></extra>",
                showlegend=True,
                orientation='v'
            ))


    # --- Dibujar Marcadores de Variantes con Tooltips Detallados y Líneas ---
    print(f"Dibujando marcadores para {len(lista_variantes_crudas_con_locus)} variantes en TV...\n")
    

    for i, var_cruda in enumerate(lista_variantes_crudas_con_locus):
        pos_cruda_0based = var_cruda.get('pos')
        type_cruda = var_cruda.get('type')
        ref_cruda = var_cruda.get('ref', '').upper()
        alt_cruda = var_cruda.get('alt', '').upper()
        locus_display = var_cruda.get('locus', 'N/A')
        hvs_region_display = var_cruda.get('hvs_region')
        
        hgvs_norm_str = variantes_hgvs_normalizadas[i] if i < len(variantes_hgvs_normalizadas) else "N/A (error HGVS)"
        mitomaster_format_str = variantes_mitomaster_formateadas[i] if i < len(variantes_mitomaster_formateadas) else "N/A (error Mitomaster)"
        empop_format_str = variantes_empop_formateadas[i] if i < len(variantes_empop_formateadas) else "N/A (error EMPOP)"


        if pos_cruda_0based is None or pos_cruda_0based < 0:
            continue

        pos_visual_1based = pos_cruda_0based + 1

        marker_color = 'red'
        marker_symbol = 'circle'
        marker_size = 10
        line_dash_style = 'dot'
        line_width = 1.5

        if type_cruda == 'substitution' or type_cruda.startswith('substitution'):
            mutation_description = f"Sustitución: {ref_cruda} > {alt_cruda}"
            marker_color = 'red'
            marker_symbol = 'circle'
        elif type_cruda == 'insertion':
            mutation_description = f"Inserción: {alt_cruda}"
            marker_color = 'blue'
            marker_symbol = 'triangle-up'
            pos_visual_1based = pos_cruda_0based + 1
            line_dash_style = 'dash'
        elif type_cruda == 'deletion':
            mutation_description = f"Deleción: {ref_cruda}"
            marker_color = 'darkorange'
            marker_symbol = 'triangle-down'
            pos_visual_1based = var_cruda.get('pos')
            line_dash_style = 'dashdot'
        
        if pos_visual_1based == 0:
            pos_visual_1based = 0.5

        # --- HTML para el tooltip (hovertemplate) - ¡SIN ALINEAMIENTO CONTEXTUAL! ---
        tooltip_content = f"<b>Variante: {ref_cruda}{pos_visual_1based}{alt_cruda if alt_cruda != '-' else 'del'}</b><br>" + \
                       f"<span style='color:#555;'>Tipo:</span> {type_cruda.capitalize()}<br>" + \
                       f"<span style='color:#555;'>Locus:</span> {locus_display}<br>"
        
        if hvs_region_display:
            tooltip_content += f"<span style='color:#555;'>Región Hipervariable:</span> {hvs_region_display}<br>"

        tooltip_content += f"<span style='color:#555;'>HGVS:</span> {hgvs_norm_str}<br>" + \
                        f"<span style='color:#555;'>Mitomaster:</span> {mitomaster_format_str}<br>" + \
                        f"<span style='color:#555;'>EMPOP:</span> {empop_format_str}"

        # Añadir marcador de la variante (Scatter trace)
        fig.add_trace(go.Scatter(
            x=[pos_visual_1based],
            y=[pistas_config_tv["Query"]["y"]],
            mode='markers',
            marker=dict(
                symbol=marker_symbol,
                size=marker_size,
                color=marker_color,
                line=dict(width=1, color='black')
            ),
            name=f"Variante: {hgvs_norm_str}",
            hovertemplate=tooltip_content, # ¡El contenido HTML va aquí!
            showlegend=False,
            # El estilo de la caja del tooltip se define en hoverlabel
            hoverlabel=dict(
                bgcolor="white",
                bordercolor="#ccc",
                font=dict(color="black", size=11),
                namelength=0 # Ocultar el nombre de la traza en el tooltip
            ),
            customdata=[var_cruda]
        ))
        
        # Añadir etiqueta de posición directamente sobre el marcador
        fig.add_annotation(
            x=pos_visual_1based,
            y=pistas_config_tv["Query"]["y"] + (pistas_config_tv["Query"]["height"] / 2) + 0.1,
            text=str(pos_visual_1based),
            showarrow=False,
            font=dict(size=8, color=marker_color, weight="bold"),
            xanchor="center",
            yanchor="bottom"
        )

        # Dibujar una línea vertical desde la variante
        fig.add_shape(
            type="line", # CORRECCIÓN: asegúrate de que esto sea "line" sin barra invertida
            x0=pos_visual_1based, y0=pistas_config_tv["Query"]["y"] - (pistas_config_tv["Query"]["height"]/2),
            x1=pistas_config_tv["D-loop"]["y"] + (pistas_config_tv["D-loop"]["height"]/2),
            line=dict(color=marker_color, width=line_width, dash=line_dash_style),
            opacity=0.7,
            layer="below",
        )

    # --- Configuración del Layout de Plotly ---
    y_vals_in_use = [cfg["y"] for cfg in pistas_config_tv.values()]
    heights_in_use = [cfg["height"] for cfg in pistas_config_tv.values()]
    
    y_min = -2.5 - 0.5 # Margen inferior para la línea base rCRS
    y_max = max(y_vals_in_use) + max(heights_in_use)/2 + 0.5 # Margen superior para el título de la pista más alta

    fig.update_layout(
        title={
            'text': f"Track Viewer de Variantes Mitocondriales: {query_id} vs {rcrs_id}",
            'yref': 'paper', 'y': 0.95, 'xanchor': 'center', 'yanchor': 'top',
            'font': dict(size=18)
        },
        xaxis=dict(
            title="Posición en rCRS (1-based)",
            range=[0, RCRS_LENGTH_EXPECTED + 1],
            showgrid=True,
            gridcolor='lightgrey',
            gridwidth=0.5,
            zeroline=False
        ),
        yaxis=dict(
            tickmode='array',
            tickvals=[pistas_config_tv[key]["y"] for key in track_order],
            ticktext=[pistas_config_tv[key]["label_full"] for key in track_order],
            range=[y_min, y_max],
            showgrid=False,
            zeroline=False,
            title="Pistas Genómicas",
            automargin=True
        ),
        hovermode="closest", # Mantiene "closest" para evitar el "shampoo"
        plot_bgcolor='white',
        margin=dict(l=150, r=50, t=100, b=50),
        # Configuración de tamaño para hacerlo responsivo
        autosize=True,
        # height=750, # Puedes eliminar estas líneas si usas autosize=True
        # width=1200, # Puedes eliminar estas líneas si usas autosize=True
        shapes=[ # Lineas horizontales para separar las pistas (background lines)
            dict(
                type="line", # CORRECCIÓN: asegúrate de que sea "line" sin barra invertida
                xref="paper", yref="y",
                x0=0, y0=pistas_config_tv[key]["y"] + pistas_config_tv[key]["height"]/2 + 0.1, x1=1, y1=pistas_config_tv[key]["y"] + pistas_config_tv[key]["height"]/2 + 0.1,
                line=dict(color="lightgrey", width=0.5, dash="dot"),
                layer="below"
            ) for key in track_order # Dibuja líneas para todas las pistas, incluyendo la de Query
        ],
        barmode='overlay',
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )

    fig.show()
    if output_html_path_tv:
        try:
            # Para que el HTML sea responsivo, Plotly lo maneja automáticamente cuando se guarda
            fig.write_html(output_html_path_tv, auto_open=False)
            print(f"Track Viewer interactivo guardado como: {output_html_path_tv}")
        except Exception as e:
            print(f"Error al guardar el Track Viewer HTML: {e}")

    print("Track Viewer interactivo generado.")