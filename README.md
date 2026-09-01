# Dinámica epigenómica durante la progresión del Alzheimer

Proyecto exploratorio que busca analizar datos de snATAC-seq (accesibilidad de cromatina a nivel de núcleo único) para rastrear cambios de accesibilidad asociados a la reactivación de elementos transponibles (TEs) en el contexto de la enfermedad de Alzheimer.

Para conocer el planteamiento, las preguntas de investigación, la metodología y
los resultados, consulta el [reporte del proyecto](docs/reporte-proyecto.md).

## Funcionalidades

- [Funcionalidad disponible 1]
- [Funcionalidad disponible 2]
- [Funcionalidad disponible 3]

<!--
El proyecto está en una etapa muy temprana. Conforme se desarrolle el pipeline, 
describan aquí solo las capacidades que ya funcionen (por ejemplo: descarga y 
preprocesamiento de matrices de snATAC-seq, generación de picos de accesibilidad, 
cruce con anotaciones de elementos transponibles, etc.).
-->

## Estructura del repositorio

```text
proyecto/
├── data/           # Datos de ejemplo o archivos pequeños
├── docs/           # Reporte y documentación
├── notebooks/      # Análisis exploratorios
├── results/        # Tablas y figuras generadas
├── src/            # Código fuente
├── tests/          # Pruebas
├── CITATION.cff    # Información para citar el software
├── LICENSE         # Licencia
└── README.md       # Introducción y guía rápida de uso
```

## Contexto y motivación

Estudios recientes muestran que, en modelos de tauopatía (ratones PS19) y en cerebros de pacientes con Alzheimer, la agregación de tau desestabiliza la heterocromatina constitutiva marcada por H3K9me3, desplazando a HP1α de los nucleosomas y provocando la erosión de la lámina nuclear. Esta pérdida de condensación reactiva elementos transponibles (principalmente LINE y LTR) normalmente silenciados, generando transcritos de doble cadena que adoptan conformación Z (Z-RNA) y activan a ZBP1, desencadenando necroptosis neuronal (Liu et al., *Nature Neuroscience*, 2026). De manera complementaria, atlas de tipo célula-específicos del tejido cerebral en Alzheimer permiten ubicar estos fenómenos dentro de la diversidad celular del cerebro.

Con base en esto, este proyecto busca usar datos de snATAC-seq para explorar si los cambios de accesibilidad de la cromatina a nivel de núcleo único coinciden con las regiones genómicas donde se localizan estos elementos transponibles, como una forma indirecta de rastrear su posible reactivación durante la progresión de la enfermedad.

## Requisitos

- [Lenguaje y versión]
- [Herramienta o biblioteca indispensable]
- [Recurso computacional o condición de acceso]

<!--
Por ejemplo: Python o R, herramientas típicas de snATAC-seq (Signac/ArchR, 
Cell Ranger ATAC, etc.), anotaciones de elementos transponibles (RepeatMasker), 
y el poder de cómputo necesario para manejar matrices de accesibilidad de célula única.
-->

## Datos

Los datos de snATAC-seq utilizados provienen del atlas celular de la iniciativa **SEA-AD** (Seattle Alzheimer's Disease Brain Cell Atlas), referenciado en el artículo del atlas celular de Alzheimer incluido en este proyecto.

[Indiquen aquí el comando o enlace de descarga y la carpeta de destino donde deben colocarse los datos.]

<!--
Falta documentar el procedimiento exacto de descarga (portal, versión del 
dataset, criterios de selección de muestras/donantes) y las anotaciones de 
elementos transponibles que se vayan a cruzar. La procedencia y características 
detalladas se describen en el reporte del proyecto.
-->

## Uso

[Expliquen qué hace el siguiente comando y qué entradas necesita.]

```bash
python src/main.py [argumentos]
```

[Indiquen dónde se guardan los resultados.]

## Reproducción de resultados

Ejecuten los siguientes pasos en el orden indicado:

```bash
python scripts/download_data.py
python src/run_analysis.py
python src/create_figures.py
```

Los resultados esperados se generarán en `results/`.

## Documentación

- [Reporte del proyecto](docs/reporte-proyecto.md)
- [Información para citar el software](CITATION.cff)

## Equipo

- Rafael Díaz Martínez — [Contribución o responsabilidad general]
- Santiago Martínez Enciso — [Contribución o responsabilidad general]
- Sixto Nezahualcoyotl Gonzales Morales — [Contribución o responsabilidad general]

## Citación

Si utilizas este software, consulta [CITATION.cff](CITATION.cff) o la opción
**Cite this repository** de GitHub.

## Licencia

[Indiquen el nombre de la licencia.] Consulta [LICENSE](LICENSE) para conocer
los términos de uso.

## Agradecimientos

[Incluyan reconocimientos institucionales o académicos]