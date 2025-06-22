# ==============================================================================
# BLOQUE 7: GENERACIÓN DE INFORME HTML Y CONVERSIÓN A PDF 
# ==============================================================================
# Descripción: Funciones para generar un informe HTML a partir de los datos
# y luego convertir ese HTML a PDF usando WeasyPrint.
# Ajustes: Ocultar columna HGVS en PDF mediante CSS @media print.
# ------------------------------------------------------------------------------
import pandas as pd
from weasyprint import HTML as WeasyHTML # Deja FontConfiguration comentado por ahora
import os
import traceback

def generar_informe_html(
    df_detallado: pd.DataFrame, 
    query_id: str,
    num_variantes_mitomaster_validas: int,
    str_variantes_mitomaster_validas: str,
    num_variantes_empop: int,
    str_empop_query: str,
    output_html_path: str = "informe_variantes_temp.html"
) -> str:
    """
    Genera un informe HTML con las tablas de variantes.
    Devuelve la ruta al archivo HTML generado.
    """
    print(f"\n--- Iniciando Generación de Informe HTML: {output_html_path} ---")

    df_html_display = df_detallado.copy()
    if 'Alineamiento' in df_html_display.columns:
        def format_alignment_for_html(text_block):
            if pd.isna(text_block) or not isinstance(text_block, str):
                return "N/A"
            escaped_text = text_block.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
            html_with_line_breaks = escaped_text.replace('\n', '<br />')
            return f"<pre style='margin:0; padding:0; white-space: pre-wrap; font-family: \"DejaVu Sans Mono\", \"Consolas\", \"Courier New\", Courier, monospace; font-size:0.85em; line-height:1.2;'>{html_with_line_breaks}</pre>"
        
        df_html_display['Alineamiento'] = df_html_display['Alineamiento'].apply(format_alignment_for_html)

    # La tabla detallada generada por pandas tendrá las clases "table", "table-striped", "table-bordered"
    tabla_detallada_html = df_html_display.to_html(index=False, escape=False, border=1, classes="table table-striped table-bordered")
    
    body_content_html = f"""
        <h1>Análisis de Variantes Mitocondriales: {query_id}</h1>

        <h2>Tabla 1: Resumen de Variantes (Formato Mitomaster)</h2>
        <table class="summary-table">
            <tr><td>ID de la Muestra:</td><td>{query_id}</td></tr>
            <tr><td>Total Variantes (Mitomaster):</td><td>{num_variantes_mitomaster_validas}</td></tr>
            <tr><td>Lista de Variantes (Mitomaster):</td><td>{str_variantes_mitomaster_validas}</td></tr>
        </table>

        <h2>Tabla 2: Resumen de Variantes (Formato EMPOP)</h2>
        <table class="summary-table">
            <tr><td>ID de la Muestra:</td><td>{query_id}</td></tr>
            <tr><td>Total Variantes (EMPOP):</td><td>{num_variantes_empop}</td></tr>
            <tr><td>Cadena de Consulta EMPOP:</td><td>{str_empop_query}</td></tr>
        </table>

        <h2>Tabla 3: Detalle Completo de Variantes</h2>
        {tabla_detallada_html}

        <p class="footer-note">* Indica variante en un hotspot de longitud conocido (interpretar con cautela en comparaciones forenses y estudios de población que dependen de tasas de mutación estándar).</p>
    """

    html_template = f"""
    <!DOCTYPE html>
    <html lang="es">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Informe de Variantes Mitocondriales: {query_id}</title>
        <style>
            body {{ 
                font-family: Arial, Helvetica, sans-serif; 
            }}
            h1, h2 {{ 
                color: #333333; 
                border-bottom: 1px solid #cccccc; 
                padding-bottom: 5px;
            }}
            h1 {{ font-size: 1.8em; text-align: center; margin-bottom: 1em;}}
            h2 {{ font-size: 1.4em; margin-top: 1.5em; margin-bottom: 0.8em;}}

            table {{ 
                border-collapse: collapse; 
                width: 100%; 
                margin-bottom: 25px; 
                font-size: 0.75em; 
            }}
            th, td {{ 
                border: 1px solid #cccccc; 
                padding: 5px; 
                text-align: left; 
                vertical-align: top; 
            }}
            th {{ 
                background-color: #e9e9e9; 
                font-weight: bold;
            }}
            .summary-table td:first-child {{ 
                font-weight: bold; 
                width: 25%; 
                background-color: #f7f7f7;
            }}
            td pre {{ 
                white-space: pre-wrap;   
                word-break: break-word;  
                font-family: "DejaVu Sans Mono", "Consolas", "Courier New", Courier, monospace; 
                font-size:0.85em;       
                line-height:1.2;        
                margin: 0;               
                padding: 0;              
            }}
            .footer-note {{ 
                font-size: 0.7em; 
                font-style: italic; 
                margin-top: 30px; 
                text-align: center; 
                color: #555555;
            }}

            @page {{
                size: A4 landscape; 
                margin: 1cm;      
            }}
            
            .table-striped tr, .table-bordered tr {{ /* Clases que usa pandas.to_html */
                page-break-inside: avoid !important;
                break-inside: avoid !important;
            }}

            /* --- NUEVA REGLA CSS PARA OCULTAR COLUMNA EN PDF/IMPRESIÓN --- */
            @media print {{
                /* "HGVS Normalizado" es la 6ª columna (índice 5 en la lista de columnas) */
                /* Asegúrate de que las clases .table-striped y .table-bordered se apliquen a tu tabla detallada */
                .table.table-striped.table-bordered th:nth-child(6),
                .table.table-striped.table-bordered td:nth-child(6) {{
                    display: none;
                }}
            }}
            /* --- FIN DE LA NUEVA REGLA CSS --- */

        </style>
    </head>
    <body>
        {body_content_html}
    </body>
    </html>
    """
    
    try:
        with open(output_html_path, "w", encoding="utf-8") as f:
            f.write(html_template)
        print(f"Informe HTML generado exitosamente: {output_html_path}")
        return output_html_path
    except Exception as e:
        print(f"Error al generar el archivo HTML: {e}")
        return None

def convertir_html_a_pdf(ruta_html: str, output_pdf_path: str):
    if not ruta_html:
        print("Error: No se proporcionó ruta al archivo HTML para la conversión a PDF.")
        return

    print(f"\n--- Iniciando Conversión de HTML a PDF: {output_pdf_path} ---")
    try:
        WeasyHTML(filename=ruta_html).write_pdf(output_pdf_path)
        print(f"Informe PDF generado exitosamente desde HTML: {output_pdf_path}")
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo HTML '{ruta_html}' para convertir a PDF.")
        print("Asegúrate de que WeasyPrint y sus dependencias (Pango, Cairo, GDK-PixBuf) estén instaladas correctamente.")
    except Exception as e:
        print(f"Error crítico al convertir HTML a PDF con WeasyPrint: {type(e).__name__} - {e}")
        import traceback
        traceback.print_exc()