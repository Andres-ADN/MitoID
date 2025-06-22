# -*- coding: utf-8 -*-
# ==============================================================================
# 1 CONSTANTS: CONFIGURACIÓN GLOBAL
# ==============================================================================
# Descripción: Este bloque contiene todas las constantes globales
# y configuraciones iniciales del entorno.
# ------------------------------------------------------------------------------


# --- Rutas a archivos de referencia ---
RCRS_FASTA_PATH = "data/NC_012920.1_rCRS.fasta"
RCRS_GENBANK_PATH = "data/NC_012920.1_rCRS.gb"

# --- Archivo de secuencia de consulta de ejemplo (si no se proporciona uno) ---
DEFAULT_QUERY_FASTA_PATH = "data/LC733703.1_full.fasta"


# --- Definición de regiones Hotspot (para informe en PDF) ---
# Formato: (inicio_1_based, fin_1_based, tipo_hotspot)
HOTSPOT_REGIONS = [
    (16184, 16193, 'indel'), (303, 315, 'indel'), (452, 455, 'indel'),
    (568, 573, 'indel'), (956, 960, 'indel'), (5895, 5899, 'indel'),
    (8272, 8276, 'indel'), (8281, 8285, 'indel')
]

# --- Regiones Hipervariables (HVS) para anotación en Track Viewer ---
# Basado en tu Tabla 1 de la tesis (pág. 6/23)
HVS_REGIONS = [
    {"nombre": "HVS-II", "inicio": 57, "fin": 372},
    {"nombre": "HVS-III", "inicio": 438, "fin": 574},
    {"nombre": "HVS-I", "inicio": 16024, "fin": 16383}
]

# --- Constantes para Nomenclatura SWGDAM/EMPOP ---
POS_3107_BLACKLISTED = 3107 # La posición 3107 no debe reportarse
HVI_C_STRETCH_START = 16182 # Ajustado para incluir inserciones como 16183.1C que EMPOP re-ancla a 16193.1C
HVI_C_STRETCH_END = 16193 # Rango del C-stretch de HVI
INSERTION_199_ANCHOR_RANGE_START = 198 # Rango para inserciones que se anclan a 199
INSERTION_199_ANCHOR_RANGE_END = 199
DELETION_249_ANCHOR_POS = 249 # Posición de anclaje para deleciones en 248/249
HVII_C_STRETCH_START = 302 # Ajustado para incluir inserciones como 302.1C que EMPOP re-ancla a 309.1C
HVII_C_STRETCH_END = 315 # Rango del C-stretch de HVII (considerando 302-310 y 311-315)

RCRS_HVI_C_COUNT_ANCHOR_16189 = 0
RCRS_HVII_C_COUNT_303_309 = 7
RCRS_HVII_C_COUNT_311_315 = 5
INSERTION_291_ANCHOR_POS = 291
INSERTION_291_DETECTION_RANGE_START = 289
INSERTION_291_DETECTION_RANGE_END = 295
AC_REPEAT_MOTIF_START = 513
AC_REPEAT_MOTIF_END = 524
AC_REPEAT_MOTIF_ANCHOR_523 = 523
AC_REPEAT_MOTIF_ANCHOR_524 = 524
MITOMASTER_HVII_DEL_CT_ANCHOR_309 = 309
MITOMASTER_HVI_INS_ANCHOR_16193 = 16193
MITOMASTER_INS_198_ANCHOR = 198
MITOMASTER_INS_291_ANCHOR = 290

# Hotspots de indels que EMPOP ignora por defecto en búsquedas
EMPOP_IGNORED_INDEL_HOTSPOTS_POS = {16193, 309, 455, 463, 573, 960, 5899, 8276, 8285}

# Mapeo de bases IUPAC extendidas (solo las más comunes para variantes)
IUPAC_MIXED_BASES = {
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
}

# ==============================================================================
# 1.1 CONSTANTS: Mapeos y constantes globales de features
# ==============================================================================

# --- Constantes para formateo de alineamiento en PDF (y usado en TV tooltip) ---
ALIGNMENT_CONTEXT_WINDOW_PDF = 13 # Bases a cada lado del evento (para PDF)
ALIGNMENT_LABEL_WIDTH_PDF = 7 # Ancho para etiquetas como "Ref:", "Query:" (para PDF)
ALIGNMENT_POS_NUM_WIDTH_PDF = 5 # Ancho para números de posición en el alineamiento (para PDF)


# --- Mapeos y constantes globales para la extracción y display de Features (rCRS) ---
# Mapeo de códigos de aminoácidos de una letra a tres letras (para referencia interna, no siempre para display)
aa_one_to_three_letter_map = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

# Mapeo de códigos de aminoácidos de una letra para display en el TV (más conciso)
aa_code_map_for_display = {
    'A': 'A', 'R': 'R', 'N': 'N', 'D': 'D', 'C': 'C',
    'Q': 'Q', 'E': 'E', 'G': 'G', 'H': 'H', 'I': 'I',
    'K': 'K', 'M': 'M', 'F': 'F', 'P': 'P', 'T': 'T',
    'W': 'W', 'Y': 'Y', 'V': 'V',
    'L': 'L', # Por defecto 'L', pero luego veremos L(UUR)/L(CUN)
    'S': 'S'  # Por defecto 'S', pero luego veremos S(UCN)/S(AGY)
}

# Estandarización de nombres de genes codificantes de proteínas (Mitomaster style)
protein_genes_map_short = {
    "MT-CO1": "CO1", "MT-COX1": "CO1", "COI": "CO1", "COX1": "CO1",
    "MT-COII": "CO2", "MT-COX2": "CO2", "COII": "CO2", "COX2": "CO2",
    "MT-COIII": "CO3", "MT-COX3": "CO3", "COIII": "CO3", "COX3": "CO3",
    "MT-CYTB": "CYTB", "MT-CYB": "CYTB", "CYTB": "CYTB", "CYB": "CYTB",
    "MT-ATP8": "ATP8", "MT-ATPASE8": "ATP8", "ATP8": "ATP8", "ATPASE8": "ATP8",
    "MT-ATP6": "ATP6", "MT-ATPASE6": "ATP6", "ATP6": "ATP6", "ATPASE6": "ATP6",
    "MT-ND1":"ND1", "ND1":"ND1", "MT-ND2":"ND2", "ND2":"ND2", "MT-ND3":"ND3", "ND3":"ND3",
    "MT-ND4L":"ND4L", "ND4L":"ND4L", "MT-ND4":"ND4", "ND4":"ND4", "MT-ND5":"ND5", "ND5":"ND5",
    "MT-ND6":"ND6", "ND6":"ND6"
}
protein_genes_final_names = set(protein_genes_map_short.values())

# Orden y mapeo de tipos de pistas finales (controla el eje Y en el TV)
final_track_types_order = ["D-loop", "tRNA", "rRNA", "gene"] # Orden de abajo a arriba en el gráfico

genbank_to_track_type_map = {
    "D-loop": "D-loop",
    "tRNA": "tRNA",
    "rRNA": "rRNA",
    "CDS": "gene",
    "gene": "gene"
}

# --- Definición MAESTRA de Regiones Genómicas de rCRS (BASADA EN LA TESIS) ---
# Esto se usará para construir las features_list en cargar_y_extraer_features_rcrs.
MASTER_RCRS_FEATURES_DEFINITION = [
    # Genes de ARNr (según tu Tabla 1, pág. 7/23 en PDF)
    {"nombre_genbank": "MT-RNR1", "nombre_display": "12S", "tipo": "rRNA", "inicio": 648, "fin": 1601, "hebra": 1},
    {"nombre_genbank": "MT-RNR2", "nombre_display": "16S", "tipo": "rRNA", "inicio": 1671, "fin": 3229, "hebra": 1},

    # Genes de ARNt (según tu Tabla 1 y 2, pág. 6-7/23 en PDF)
    {"nombre_genbank": "MT-TF", "nombre_display": "F", "tipo": "tRNA", "inicio": 577, "fin": 647, "hebra": 1},
    {"nombre_genbank": "MT-TV", "nombre_display": "V", "tipo": "tRNA", "inicio": 1602, "fin": 1670, "hebra": 1},
    {"nombre_genbank": "MT-TL1", "nombre_display": "L(UUR)", "tipo": "tRNA", "inicio": 3230, "fin": 3304, "hebra": 1},
    {"nombre_genbank": "MT-TI", "nombre_display": "I", "tipo": "tRNA", "inicio": 4263, "fin": 4331, "hebra": 1},
    {"nombre_genbank": "MT-TQ", "nombre_display": "Q", "tipo": "tRNA", "inicio": 4329, "fin": 4400, "hebra": -1}, # Es (L) complementaria
    {"nombre_genbank": "MT-TM", "nombre_display": "M", "tipo": "tRNA", "inicio": 4402, "fin": 4469, "hebra": 1},
    {"nombre_genbank": "MT-TW", "nombre_display": "W", "tipo": "tRNA", "inicio": 5512, "fin": 5579, "hebra": 1},
    {"nombre_genbank": "MT-TA", "nombre_display": "A", "tipo": "tRNA", "inicio": 5587, "fin": 5655, "hebra": -1}, # Es (L) complementaria
    {"nombre_genbank": "MT-TN", "nombre_display": "N", "tipo": "tRNA", "inicio": 5657, "fin": 5729, "hebra": -1}, # Es (L) complementaria
    {"nombre_genbank": "MT-TC", "nombre_display": "C", "tipo": "tRNA", "inicio": 5761, "fin": 5826, "hebra": -1}, # Es (L) complementaria
    {"nombre_genbank": "MT-TY", "nombre_display": "Y", "tipo": "tRNA", "inicio": 5826, "fin": 5891, "hebra": -1}, # Es (L) complementaria
    {"nombre_genbank": "MT-TS1", "nombre_display": "S(UCN)", "tipo": "tRNA", "inicio": 7446, "fin": 7514, "hebra": -1}, # Es (L) complementaria
    {"nombre_genbank": "MT-TD", "nombre_display": "D", "tipo": "tRNA", "inicio": 7518, "fin": 7585, "hebra": 1},
    {"nombre_genbank": "MT-TK", "nombre_display": "K", "tipo": "tRNA", "inicio": 8295, "fin": 8364, "hebra": 1},
    {"nombre_genbank": "MT-TG", "nombre_display": "G", "tipo": "tRNA", "inicio": 9991, "fin": 10058, "hebra": 1},
    {"nombre_genbank": "MT-TR", "nombre_display": "R", "tipo": "tRNA", "inicio": 10405, "fin": 10469, "hebra": 1},
    {"nombre_genbank": "MT-TH", "nombre_display": "H", "tipo": "tRNA", "inicio": 12138, "fin": 12206, "hebra": 1},
    {"nombre_genbank": "MT-TS2", "nombre_display": "S(AGY)", "tipo": "tRNA", "inicio": 12207, "fin": 12265, "hebra": 1},
    {"nombre_genbank": "MT-TL2", "nombre_display": "L(CUN)", "tipo": "tRNA", "inicio": 12266, "fin": 12336, "hebra": 1},
    {"nombre_genbank": "MT-TE", "nombre_display": "E", "tipo": "tRNA", "inicio": 14674, "fin": 14742, "hebra": -1}, # Es (L) complementaria
    {"nombre_genbank": "MT-TT", "nombre_display": "T", "tipo": "tRNA", "inicio": 15888, "fin": 15953, "hebra": 1},
    {"nombre_genbank": "MT-TP", "nombre_display": "P", "tipo": "tRNA", "inicio": 15956, "fin": 16023, "hebra": -1}, # Es (L) complementaria

    # Genes codificantes de proteínas (según tu Tabla 1, pág. 6/23 en PDF)
    {"nombre_genbank": "MT-ND1", "nombre_display": "ND1", "tipo": "gene", "inicio": 3307, "fin": 4262, "hebra": 1},
    {"nombre_genbank": "MT-ND2", "nombre_display": "ND2", "tipo": "gene", "inicio": 4470, "fin": 5511, "hebra": 1},
    {"nombre_genbank": "MT-CO1", "nombre_display": "CO1", "tipo": "gene", "inicio": 5904, "fin": 7445, "hebra": 1},
    {"nombre_genbank": "MT-CO2", "nombre_display": "CO2", "tipo": "gene", "inicio": 7586, "fin": 8269, "hebra": 1},
    {"nombre_genbank": "MT-ATP8", "nombre_display": "ATP8", "tipo": "gene", "inicio": 8366, "fin": 8572, "hebra": 1},
    {"nombre_genbank": "MT-ATP6", "nombre_display": "ATP6", "tipo": "gene", "inicio": 8527, "fin": 9207, "hebra": 1},
    {"nombre_genbank": "MT-CO3", "nombre_display": "CO3", "tipo": "gene", "inicio": 9207, "fin": 9990, "hebra": 1},
    {"nombre_genbank": "MT-ND3", "nombre_display": "ND3", "tipo": "gene", "inicio": 10059, "fin": 10404, "hebra": 1},
    {"nombre_genbank": "MT-ND4L", "nombre_display": "ND4L", "tipo": "gene", "inicio": 10470, "fin": 10766, "hebra": 1},
    {"nombre_genbank": "MT-ND4", "nombre_display": "ND4", "tipo": "gene", "inicio": 10760, "fin": 12137, "hebra": 1},
    {"nombre_genbank": "MT-ND5", "nombre_display": "ND5", "tipo": "gene", "inicio": 12337, "fin": 14148, "hebra": 1},
    {"nombre_genbank": "MT-ND6", "nombre_display": "ND6", "tipo": "gene", "inicio": 14149, "fin": 14673, "hebra": -1}, # Es (L) complementaria
    {"nombre_genbank": "MT-CYB", "nombre_display": "CYTB", "tipo": "gene", "inicio": 14747, "fin": 15887, "hebra": 1},
]



# Constantes para anotación de locus
PRIORITY_FEATURE_TYPES = ["CDS", "rRNA", "tRNA"]
