# ==============================================================================
# BLOQUE 2: UTILIDADES DE CARGA Y EXTRACCIÓN DE CARACTERÍSTICAS (FEATURES)
# ==============================================================================
# Descripción: Funciones dedicadas a la carga de secuencias desde archivos
# FASTA y a la extracción y procesamiento de características (features)
# desde archivos GenBank, específicamente para la secuencia rCRS.
# Se ha refactorizado la lógica de mapeo de features a tipos de pistas
# para asegurar la correcta clasificación de Genes, rRNAs y tRNAs,
# utilizando las nomenclaturas del GenBank y las de Mitomaster/SWGDAM.
# ------------------------------------------------------------------------------



# ==============================================================================
# Importaciones necesarias

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation # CompoundLocation es importante para D-loop
import os
from . import constants #importaciones de constants en la misma carpeta
# ==============================================================================

def cargar_secuencia_fasta(ruta_archivo: str, id_esperado: str = None) -> SeqRecord | None:
    """Carga una secuencia desde un archivo FASTA."""
    try:
        record = SeqIO.read(ruta_archivo, "fasta")
        if id_esperado and record.id != id_esperado:
            print(f"Advertencia: El ID ('{record.id}') no coincide con el esperado ('{id_esperado}').")
        print(f"Secuencia '{record.id}' cargada desde '{ruta_archivo}' (Longitud: {len(record.seq)} pb).\n")
        return record
    except FileNotFoundError:
        print(f"Error: No se encontró FASTA en: '{ruta_archivo}'")
        return None
    except ValueError:
        print(f"Error: FASTA '{ruta_archivo}' mal formateado o con >1 secuencia.")
        return None
    except Exception as e:
        print(f"Error al leer '{ruta_archivo}': {e}")
        return None

def cargar_y_extraer_features_rcrs(gb_path: str) -> list:
    """
    Carga un archivo GenBank de rCRS y extrae las características genómicas relevantes,
    estandarizando sus nombres utilizando una definición maestra basada en la tesis,
    incluyendo la hebra y la longitud para el Track Viewer.
    """
    try:
        rcrs_record_gb = SeqIO.read(gb_path, "genbank")
        print(f"Archivo GenBank rCRS '{rcrs_record_gb.id}' cargado. Extrayendo features...")
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo GenBank rCRS en: '{gb_path}'")
        return []
    except Exception as e:
        print(f"Error cargando el archivo GenBank rCRS '{gb_path}': {e}")
        return []

    features_list = []
    
    # Lista para almacenar features que ya hemos "usado" del GenBank para evitar duplicados.
    # se usa (type, start, end) para identificar unicidad.
    processed_genbank_features_keys = set()

    # --- 1. Procesar D-loop del GenBank (manejo especial por CompoundLocation) ---
    for feature in rcrs_record_gb.features:
        if feature.type == "D-loop":
            if isinstance(feature.location, CompoundLocation):
                for part_location in feature.location.parts:
                    part_start_1based = int(part_location.start) + 1
                    part_end_1based = int(part_location.end)
                    feature_key = ("D-loop", part_start_1based, part_end_1based)
                    if feature_key not in processed_genbank_features_keys:
                        features_list.append({
                            "nombre": "D-loop",
                            "inicio": part_start_1based,
                            "fin": part_end_1based,
                            "tipo": "D-loop",
                            "hebra": part_location.strand,
                            "longitud": part_end_1based - part_start_1based + 1
                        })
                        processed_genbank_features_keys.add(feature_key)
            else: # Simple D-loop feature
                start_1based = int(feature.location.start) + 1
                end_1based = int(feature.location.end)
                feature_key = ("D-loop", start_1based, end_1based)
                if feature_key not in processed_genbank_features_keys:
                    features_list.append({
                        "nombre": "D-loop",
                        "inicio": start_1based,
                        "fin": end_1based,
                        "tipo": "D-loop",
                        "hebra": feature.location.strand,
                        "longitud": end_1based - start_1based + 1
                    })
                    processed_genbank_features_keys.add(feature_key)
            # No break, ya que puede haber múltiples features de D-loop (aunque rCRS suele tener 1 Compound)
            
    # --- 2. Procesar los Genes, rRNAs y tRNAs de la DEFINICIÓN MAESTRA (TU TESIS) ---
    for master_feat in constants.MASTER_RCRS_FEATURES_DEFINITION:
        display_name = master_feat["nombre_display"]
        final_track_type = master_feat["tipo"]
        start_pos = master_feat["inicio"]
        end_pos = master_feat["fin"]
        hebra_info = master_feat["hebra"]
        
        # Verificar si este feature de la tesis ya fue cubierto por una parte del D-loop GenBank (poco probable)
        feature_key = (final_track_type, start_pos, end_pos)
        if feature_key in processed_genbank_features_keys:
            continue # Ya procesado, saltar

        features_list.append({
            "nombre": display_name,
            "inicio": start_pos,
            "fin": end_pos,
            "tipo": final_track_type,
            "hebra": hebra_info,
            "longitud": end_pos - start_pos + 1
        })
        processed_genbank_features_keys.add(feature_key)


    # --- 3. Ordenar la lista final de features para la visualización ---
    # Re-ordenar por el orden de pista final (final_track_types_order) y luego por posición de inicio.
    features_list.sort(key=lambda x: (constants.final_track_types_order.index(x["tipo"]) if x["tipo"] in constants.final_track_types_order else 99, x["inicio"]))

    print(f"Se extrajeron {len(features_list)} características relevantes de rCRS GenBank.")
    return features_list
    