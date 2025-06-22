# MitoID: Herramienta de Análisis y Reporte de Variantes de ADN Mitocondrial

## 🧬 Visión General

MitoID es una herramienta bioinformática de código abierto, desarrollada para facilitar el análisis estandarizado y la comprensión visual de variantes en el ADN mitocondrial (ADNmt) humano. Diseñada para ser utilizada en contextos académicos y forenses, MitoID aborda la necesidad de una solución integrada, accesible y didáctica para el estudio de mutaciones del ADNmt a partir de secuencias en formato FASTA.

Este proyecto surge de la investigación realizada para el Trabajo de Fin de Máster en Genética, Física y Química Forense, con el objetivo de ofrecer una alternativa robusta y fácil de usar para la comunidad científica.

## ✨ Características Principales

* **Análisis Estandarizado:** Alineamiento de secuencias FASTA contra la Secuencia de Referencia de Cambridge Revisada (rCRS, NC_012920.1) para la identificación precisa de SNPs e indels.
* **Filtrado Inteligente de Artefactos:** Incluye lógicas avanzadas para el manejo de artefactos comunes, como la ambigua posición 3107N, garantizando la fiabilidad de los perfiles genéticos.
* **Gestión de Regiones de Interés:** Reconocimiento y procesamiento específico de regiones hipervariables (HVS-I, HVS-II, HVS-III), motivos repetitivos AC (np 515-525) y C-stretches.
* **Nomenclatura Consistente:** Estandarización automática de variantes bajo las convenciones HGVS y EMPOP/SWGDAM, facilitando la comparabilidad y comunicación de resultados.
* **Reportes Detallados:** Generación de informes en formato HTML y PDF, incluyendo tablas claras con el detalle de las variantes y su anotación.
* **Visualización Interactiva:** Un "Track Viewer" dinámico que permite explorar la posición de cada mutación en relación con las regiones codificantes y de control del genoma, ideal para fines didácticos y de investigación.

## 🚀 Instalación y Uso (para Usuarios de Linux/WSL)

MitoID está diseñado para funcionar en un entorno local y requiere Python 3.12. Se recomienda encarecidamente utilizar un entorno virtual para gestionar las dependencias de manera limpia.

### Requisitos Previos:

* **Python 3.12:** Asegúrate de tener Python 3.12 instalado en tu sistema.
* **Miniconda o Anaconda:** Herramientas recomendadas para la gestión de entornos y paquetes de Python. Si no las tienes, puedes instalarlas siguiendo las instrucciones en [Miniconda](https://docs.conda.io/en/latest/miniconda.html) o [Anaconda](https://www.anaconda.com/download).
* **Dependencias de WeasyPrint (solo Linux):** Para la generación de PDFs, `WeasyPrint` necesita algunas librerías del sistema operativo. Si usas Linux, instálalas con:
    ```bash
    sudo apt-get update
    sudo apt-get install libpango-1.0-0 libcairo2 libgdk-pixbuf2.0-0
    ```
    (Si usas otra distribución de Linux, busca cómo instalar `pango`, `cairo` y `gdk-pixbuf` para tu sistema).

### Pasos de Instalación:

1.  **Clonar el Repositorio:**
    ```bash
    git clone [https://github.com/Andres-ADN/MitoID.git]
    cd MitoID
    ```
    

2.  **Crear y Activar el Entorno Virtual (Recomendado):**
    Utilizaremos `conda` para crear un entorno limpio y aislado para MitoID.

    ```bash
    # (Opcional) Desactivar cualquier entorno activo
    conda deactivate
    # Eliminar cualquier versión anterior del entorno MitoID_conda (si existe)
    conda env remove -n MitoID_conda -y
    # Crear un nuevo entorno Python 3.12 con los canales recomendados
    conda create -n MitoID_conda python=3.12 -c conda-forge -c bioconda -c defaults -y
    # Activar el entorno
    conda activate MitoID_conda
    ```

3.  **Instalar las Dependencias de Python:**
    Una vez activado el entorno, instala todas las librerías Python necesarias desde el archivo `requirements.txt` generado:

    ```bash
    pip install -r requirements.txt
    ```

### Cómo Ejecutar MitoID:

1.  **Asegúrate de que tu entorno `MitoID_conda` esté activado:**
    ```bash
    conda activate MitoID_conda
    ```

2.  **Coloca tus archivos de referencia y consulta:**
    * Los archivos de referencia **rCRS** (`NC_012920.1_rCRS.fasta` y `NC_012920.1_rCRS.gb`) **deben estar siempre** en la carpeta `data/`.
    * Coloca tus archivos FASTA de las muestras a analizar (las "query") también en la carpeta `data/`. Por ejemplo, `data/mi_muestra.fasta`.
    * MitoID incluye una secuencia de ejemplo (`data/LC733703.1_full.fasta`) que se usará automáticamente si no proporcionas un archivo de consulta o si el archivo proporcionado no se encuentra.

3.  **Ejecuta el Script Principal:**
    Desde la carpeta raíz del proyecto (`MitoID/`), puedes ejecutar la aplicación de las siguientes maneras:

    * **Para analizar una secuencia específica (recomendado):**
        ```bash
        python src/main.py data/nombre_de_tu_muestra.fasta
        ```
        *Ejemplo:*
        ```bash
        python src/main.py data/MW389258_1.fasta
        ```

    * **Para ejecutar con la secuencia de ejemplo por defecto:**
        ```bash
        python src/main.py
        ```
        (MitoID detectará que no se proporcionó un archivo y usará `data/LC733703.1_full.fasta`).

### Resultados:

MitoID generará archivos de salida en la carpeta `output/` de tu proyecto:
* Un informe HTML detallado (ej. `temp_informe_variantes_MW389258.1.html`)
* Un informe PDF profesional (ej. `Informe_Variantes_MW389258.1_vs_rCRS.pdf`)
* Un visualizador interactivo de variantes en formato HTML (ej. `Track_Viewer_MW389258.1_vs_rCRS.html`)

Abre los archivos HTML en tu navegador web para explorar los resultados interactivos.

## 🤝 Contribuciones y Soporte

MitoID es un proyecto de código abierto. Si deseas contribuir o reportar problemas, por favor visita el repositorio de GitHub .


## 🎓 Agradecimientos

Este trabajo ha sido desarrollado como parte del Trabajo de Fin de Máster en Genética, Física y Química Forense, en la Universitat Rovira i Virgili. Agradezco a mi tutor y a la comunidad científica por su inspiración y apoyo.
