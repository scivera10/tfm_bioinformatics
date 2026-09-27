# tfm_snakemake

Flujo de trabajo bioinformático desarrollado en Snakemake para el ensamblaje genómico y la anotación de elementos repetitivos en plantas, a partir de datos de secuenciación PacBio HiFi e Illumina Hi-C.

Este repositorio forma parte del Trabajo de Fin de Máster (TFM) **"Diseño e implementación de un flujo de trabajo bioinformático para el ensamblaje genómico y la anotación de elementos repetitivos en plantas"**.

- **Autor:** Sergio Civera Arroyo
- **Año:** 2026

## Descripción del flujo de trabajo

El pipeline está construido con [Snakemake](https://snakemake.readthedocs.io/) y organizado en cuatro módulos que agrupan las reglas según la etapa del análisis a la que pertenecen:

1. **Preprocesado y control de calidad** (`rules/preprocessing.smk`): control de calidad de las lecturas HiFi crudas, filtrado y recorte de las lecturas Hi-C, y análisis de k-mers para la estimación de las características del genoma.
2. **Ensamblaje de contigs** (`rules/assembly.smk`): ensamblaje primario del genoma en modo *phased* a partir de las lecturas HiFi y Hi-C, evaluación de la calidad del ensamblaje, y eliminación de la redundancia haplotípica.
3. **Mapeo Hi-C y scaffolding cromosómico** (`rules/arima_yahs.smk`): mapeo de las lecturas Hi-C contra el ensamblaje purgado siguiendo el protocolo de Arima Genomics, y scaffolding cromosómico del ensamblaje con YaHS.
4. **Anotación de elementos repetitivos** (`rules/repeat_annotation.smk`): identificación y anotación de retrotransposones LTR, elementos repetitivos de novo, microsatélites y repeticiones en tándem.

Cada regla se ejecuta en su propio entorno conda (definido en `envs/`), lo que garantiza la reproducibilidad y evita conflictos entre las dependencias de las distintas herramientas.

## Estructura del repositorio

```
tfm_snakemake/
├── Snakefile              # Punto de entrada del pipeline; importa los cuatro módulos de reglas
├── config/
│   └── config.yaml        # Rutas a los datos de entrada y parámetros de cada herramienta
├── rules/
│   ├── preprocessing.smk       # Control de calidad y análisis de k-mers
│   ├── assembly.smk            # Ensamblaje de contigs, evaluación y purga de duplicados
│   ├── arima_yahs.smk          # Mapeo Hi-C y scaffolding cromosómico
│   └── repeat_annotation.smk   # Anotación de elementos repetitivos
├── envs/                  # Un archivo .yaml por entorno conda, con las herramientas y versiones fijadas
├── scripts/
│   ├── arima_mapping_pipeline/ # Scripts Perl del protocolo de mapeo Hi-C de Arima Genomics
│   ├── hifiasm/                 # Script auxiliar para convertir la salida de hifiasm (GFA) a FASTA
│   ├── misa/                    # Herramienta MISA para la identificación de microsatélites
│   └── purge_dups/              # Script para generar el histograma de profundidad de cobertura
├── data/                  # Datos de entrada (no incluidos en el repositorio; ver más abajo)
├── logs/                  # Logs de ejecución de las reglas (se genera al ejecutar el pipeline)
└── results/               # Resultados del pipeline (se genera al ejecutar el pipeline)
```

### Descripción de las carpetas

- **`config/`**: contiene `config.yaml`, donde se especifican las rutas a los datos de entrada (lecturas HiFi y Hi-C) y los parámetros específicos de cada herramienta del pipeline.
- **`rules/`**: contiene los cuatro módulos `.smk` en los que se ha organizado el flujo de trabajo, cada uno agrupando las reglas de una etapa del análisis.
- **`envs/`**: contiene un archivo `.yaml` por cada entorno conda necesario para ejecutar las reglas del pipeline, con las versiones exactas de cada herramienta.
- **`scripts/`**: contiene los scripts auxiliares empleados por algunas reglas, organizados en subcarpetas por la herramienta a la que pertenecen.
- **`data/`**: carpeta donde deben colocarse los datos de entrada (lecturas HiFi en `data/hifi/` y lecturas Hi-C en `data/hic/`). No se incluye en el repositorio por el tamaño de los archivos.
- **`logs/`** y **`results/`**: se generan automáticamente al ejecutar el pipeline y contienen, respectivamente, los registros de ejecución y los resultados de cada regla.

## Ejecución

El pipeline se ha desarrollado y probado en un entorno de Google Colab, con Snakemake instalado mediante mamba y con conda (a través de condacolab) como gestor de los entornos por regla.

Antes de ejecutar el pipeline, es necesario:

1. Colocar los datos de entrada en `data/hifi/` y `data/hic/`.
2. Configurar las rutas y parámetros en `config/config.yaml`.

Para comprobar que el flujo de trabajo está correctamente definido, sin ejecutar ningún comando real:

```bash
snakemake -n --use-conda
```

Para ejecutar el pipeline completo, indicando el número de núcleos disponibles (por ejemplo, 4):

```bash
snakemake --use-conda --cores 4
```

Snakemake se encargará de crear automáticamente los entornos conda necesarios para cada regla y de gestionar el orden de ejecución de las tareas en función de sus dependencias.
