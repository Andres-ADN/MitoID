# ==============================================================================
# BLOQUE 4: ANOTACIÓN DE LOCUS Y UTILIDADES DE HOTSPOT
# ==============================================================================
# Descripción: Funciones para anotar la región o gen mitocondrial afectado
# por una variante y para identificar si una variante cae en una región hotspot.
# ------------------------------------------------------------------------------

from . import constants

def anotar_locus_variante(pos_variante_1based: int, tipo_variante: str, longitud_evento: int, features_rcrs: list) -> str:
    """
    Anota el locus de una variante basándose en su posición y las features de rCRS.
    'longitud_evento' es len(alt) para inserciones/sustituciones, len(ref) para deleciones.
    """
    if not features_rcrs: return "N/A (Features no disponibles)"
    if not isinstance(pos_variante_1based, int) or pos_variante_1based <= 0: return "N/A (Posición inválida)"

    # Determinar el rango de posiciones afectado por la variante en coordenadas 1-based
    start_affected_pos = pos_variante_1based
    end_affected_pos = pos_variante_1based

    if tipo_variante == 'insertion':
        # Una inserción ocurre *después* de pos_variante_1based (si se ancla a la base anterior)
        # o *entre* pos_variante_1based y pos_variante_1based+1.
        # Para la búsqueda de locus, consideramos que afecta el punto de inserción.
        # Si la pos_variante_1based es la base *anterior* a la inserción, entonces la inserción
        # ocurre entre pos_variante_1based y (pos_variante_1based + 1).
        # Consideramos que la variante "toca" ambas posiciones para el solapamiento.
        start_affected_pos = pos_variante_1based 
        end_affected_pos = pos_variante_1based + 1 # El "espacio" donde se inserta
    elif tipo_variante == 'deletion':
        end_affected_pos = pos_variante_1based + longitud_evento - 1
    elif tipo_variante in ["transition", "transversion", "substitution", "substitution (con N)"]:
        # Para sustituciones, solo afecta a la posición dada.
        pass # start_affected_pos y end_affected_pos ya están correctos.
    
    regiones_solapantes = []
    for feature in features_rcrs:
        # Comprobar solapamiento: !(feature_fin < var_inicio || feature_inicio > var_fin)
        if not (feature['fin'] < start_affected_pos or feature['inicio'] > end_affected_pos):
            regiones_solapantes.append(feature)

    if not regiones_solapantes:
        # Si no hay solapamiento directo, buscar regiones intergénicas.
        feature_anterior = None
        feature_siguiente = None
        for f_idx, f_val in enumerate(features_rcrs):
            if f_val['fin'] < start_affected_pos: # Completamente antes
                if feature_anterior is None or f_val['fin'] > feature_anterior['fin']:
                    feature_anterior = f_val
            elif f_val['inicio'] > end_affected_pos: # Completamente después
                if feature_siguiente is None or f_val['inicio'] < feature_siguiente['inicio']:
                    feature_siguiente = f_val
                break # Las features están ordenadas, no es necesario seguir.
        
        if feature_anterior and feature_siguiente:
            # Caso especial: inserción exactamente entre dos genes contiguos
            if tipo_variante == 'insertion' and \
               start_affected_pos == feature_anterior['fin'] and \
               end_affected_pos == feature_siguiente['inicio']:
                return f"Entre ({feature_anterior['nombre']}, {feature_siguiente['nombre']})"
            
            # Si está cerca de ambas, es intergénica. Distancia umbral (e.g. 10bp)
            distancia_intergenica = feature_siguiente['inicio'] - feature_anterior['fin'] -1 
            if 0 <= distancia_intergenica <= 10 : # Espacio intergénico pequeño
                 return f"Intergénica ({feature_anterior['nombre']}-{feature_siguiente['nombre']})"
            else: # Espacio intergénico grande o no definido claramente
                return "Región no anotada (intergénica amplia)"
        elif feature_anterior and start_affected_pos > feature_anterior['fin']:
            return f"Posterior a {feature_anterior['nombre']}"
        elif feature_siguiente and end_affected_pos < feature_siguiente['inicio']:
            return f"Anterior a {feature_siguiente['nombre']}"
        return "Región desconocida (extremos del genoma o error)"

    # Priorizar tipos funcionales (CDS, rRNA, tRNA) si hay múltiples solapamientos
    tipos_prioritarios = constants.PRIORITY_FEATURE_TYPES
    regiones_prioritarias = [f['nombre'] for f in regiones_solapantes if f['tipo'] in tipos_prioritarios]
    
    nombres_finales = []
    if regiones_prioritarias:
        nombres_finales = list(set(regiones_prioritarias)) # Nombres únicos
    else: # Si no hay funcionales, usar cualquier tipo de región solapante (e.g. D-loop)
        nombres_finales = list(set(f['nombre'] for f in regiones_solapantes))

    nombres_finales.sort() # Ordenar alfabéticamente para consistencia

    if not nombres_finales: return "Región no anotada (solapamiento sin nombre)"
    return ", ".join(nombres_finales)


def es_variante_en_hotspot(posicion_variante: int, tipo_variante_cruda: str, longitud_evento: int = 1) -> bool:
    """
    Verifica si una variante (indel) cae dentro de una región hotspot definida.
    'posicion_variante' es 1-based. Para inserciones, es la base *anterior*.
    'longitud_evento' es relevante para deleciones.
    """
    if tipo_variante_cruda not in ["insertion", "deletion"]:
        return False
    if not isinstance(posicion_variante, int): return False

    # Las regiones HOTSPOT_REGIONS son 1-based.
    # Para inserciones: la posición 'posicion_variante' es la base 0-based *antes* de la inserción,
    # por lo que la inserción ocurre *después* de (posicion_variante + 1) en 1-based.
    # Para deleciones: 'posicion_variante' es el inicio 1-based de la deleción.
    
    start_affected_region = -1
    end_affected_region = -1

    if tipo_variante_cruda == "insertion":
        # Inserción ocurre *después* de 'posicion_variante' (si 0-based),
        # o entre 'posicion_variante' y 'posicion_variante + 1' (si 1-based).
        # Si la 'posicion_variante' que llega es 0-based (como en variantes crudas),
        # el punto de inserción 1-based es (posicion_variante + 1).
        # El evento de inserción se considera en la "unión" o "gap" entre bases.
        # Para hotspots de indel, se suele referenciar a las bases adyacentes.
        # ej. Inserción en 309.1C -> hotspot 303-315. pos_cruda es 301 (0-based) -> 302 (1-based)
        # Si la 'posicion_variante' es la pos 0-based ANTERIOR (lo que usa 'extraer_variantes_crudas'), entonces el punto de inserción es (posicion_variante + 1)
        # Este es el punto donde se inserta, así que es el que se compara.
        # La 'posicion_variante' que llega a 'generar_tabla_y_exportar' para 'es_variante_en_hotspot'
        # viene de `var_cruda.get('pos')`. Para inserciones, esto es 0-based.
        # Entonces, el punto real de inserción en 1-based es `posicion_variante + 1`.
        point_of_insertion_1based = posicion_variante + 1 
        start_affected_region = point_of_insertion_1based
        end_affected_region = point_of_insertion_1based # Un punto
    elif tipo_variante_cruda == "deletion":
        # 'posicion_variante' es el inicio 1-based de la deleción (como en variantes crudas)
        start_affected_region = posicion_variante
        end_affected_region = posicion_variante + longitud_evento - 1
    
    if start_affected_region == -1: return False

    for hs_start, hs_end, _ in constants.HOTSPOT_REGIONS:
        # Comprobar si la región afectada por la variante [start_affected_region, end_affected_region]
        # solapa con la región hotspot [hs_start, hs_end].
        # Solapamiento si: !(hs_end < start_affected_region || hs_start > end_affected_region)
        if not (hs_end < start_affected_region or hs_start > end_affected_region):
            return True
    return False

def obtener_hvs_region(pos_variante_1based: int) -> str | None:
    """
    Determina si una posición cae en una región hipervariable (HVS-I, HVS-II, HVS-III).
    Retorna el nombre de la región (ej. "HVS-I") o None si no cae en ninguna.
    """
    for hvs in constants.HVS_REGIONS: 
        if hvs["inicio"] <= pos_variante_1based <= hvs["fin"]:
            return hvs["nombre"]
    return None