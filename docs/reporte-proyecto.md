# Dinámica epigenómica de elementos transponibles durante la progresión del Alzheimer: análisis de snATAC-seq de SEA-AD

## Información general

| Dato | Información |
|:--|:--|
| Integrante 1 | [Rafael Díaz Martínez — rafadiaz@lcg.unam.mx] |
| Integrante 2 | [Santiago Martínez Enciso — santiago.marenc@gmail.com] |
| Integrante 3 | [Sixto Netzahualcóyotl González Morales — sixtogm@lcg.unam.mx] |
| Fecha de creación | [25/08/2026] |
| Última actualización | [20/10/2026] |
| Repositorio | [https://github.com/rafadmtz/NeuroEpigenomics_AD] |
| Estado | [En desarrollo: pipeline de preprocesamiento implementado (scripts 00–04); faltan el cruce con elementos transponibles y el análisis de accesibilidad diferencial] |

## Resumen del proyecto

<!-- AYUDA: Expliquen en un párrafo el problema, su relevancia, los datos, la
solución propuesta y el resultado principal. Conviene actualizarlo al final.
EJEMPLO: Se analizarán genomas de E. coli disponibles en NCBI para identificar
y comparar genes de resistencia. Se desarrollará un flujo reproducible en Python
que obtenga los datos y genere tablas y visualizaciones comparativas. -->

En modelos de tauopatía y en cerebros de pacientes con enfermedad de Alzheimer (EA), la agregación de tau desestabiliza la heterocromatina y reactiva elementos transponibles (TEs) normalmente silenciados. Estos TEs producen RNA de doble cadena que activa la vía Z-RNA–ZBP1 y la muerte neuronal (Liu et al., 2026). Este proyecto explora si esa pérdida de silenciamiento se refleja en la accesibilidad de la cromatina a nivel de núcleo único. Para ello usa los datos de snATAC-seq del giro temporal medio (MTG) del atlas SEA-AD (Gabitto et al., 2024) y la anotación de repeticiones de RepeatMasker para hg38. Se desarrolla en Python un flujo reproducible que: descarga los datos; integra y valida los metadatos de núcleos y donadores; agrupa los núcleos por estadio neuropatológico (ADNC) y subclase celular; genera perfiles pseudobulk por donador para un análisis de accesibilidad diferencial; y produce visualizaciones de composición celular y embeddings UMAP de las subclases L2/3 IT, Sst y Vip. El paso siguiente es cruzar los picos de accesibilidad con las coordenadas de los TEs.

## 1. Contexto y antecedentes

<!-- AYUDA: Presenten la información necesaria para comprender el proyecto.
Definan conceptos biológicos y computacionales, describan qué se conoce y citen
trabajos, datos o herramientas relacionados.
PREGUNTAS: ¿Cuál es el fenómeno de interés? ¿Qué debe conocer quien lea el
reporte? ¿Qué métodos o herramientas se han utilizado anteriormente?
EJEMPLO: La resistencia antimicrobiana es un problema de salud pública. Aunque
existen bases especializadas, comparar varios genomas requiere integrar datos
procedentes de distintos archivos. -->

Elementos transponibles y heterocromatina. Cerca de la mitad del genoma humano está formada por secuencias derivadas de elementos transponibles: retrotransposones LINE, SINE y LTR (incluidos los retrovirus endógenos), elementos SVA y transposones de DNA. En las células somáticas la mayoría permanece silenciada en heterocromatina constitutiva marcada por H3K9me3 y unida a HP1α.

Tau y reactivación de TEs. En el modelo de tauopatía PS19 y en cerebros de pacientes con EA, la agregación de tau desestabiliza esa heterocromatina: desplaza a HP1α de los nucleosomas y erosiona la lámina nuclear. La pérdida de condensación reactiva TEs normalmente silenciados, principalmente LINE y LTR. Sus transcritos de doble cadena adoptan conformación Z (Z-RNA) y activan a ZBP1, lo que desencadena necroptosis neuronal (Liu et al., 2026).

snATAC-seq. La secuenciación de cromatina accesible por transposasa en núcleos únicos (snATAC-seq) mide qué regiones del genoma están abiertas en cada núcleo. El resultado es una matriz dispersa de núcleos × picos de accesibilidad. Como la apertura de la cromatina precede o acompaña a la transcripción, un aumento de accesibilidad sobre TEs puede servir como indicador indirecto de su desrepresión.

Atlas SEA-AD. El Seattle Alzheimer's Disease Brain Cell Atlas (Gabitto et al., 2024) perfiló el MTG de donadores que abarcan todo el espectro neuropatológico de la EA, con anotación celular jerárquica (clase, subclase, supertipo). Cada donador está caracterizado con medidas neuropatológicas (ADNC, Braak, Thal, CERAD) y un Continuous Pseudo-progression Score (CPS) que los ordena en un continuo de progresión. Esto permite ubicar los cambios epigenómicos dentro de tipos celulares específicos.

Herramientas previas. Los datos de célula única se manejan con el ecosistema AnnData/scanpy (Wolf et al., 2018). Para ATAC, SnapATAC2 (Zhang et al., 2024) ofrece reducción de dimensionalidad espectral adaptada a matrices binarias y dispersas. Para comparar condiciones entre donadores, el enfoque recomendado es agregar los conteos por donador (pseudobulk) y modelarlos con DESeq2 o edgeR. Así se evita la pseudorreplicación de los análisis a nivel de núcleo (Squair et al., 2021). Las coordenadas de los TEs se obtienen de la tabla rmsk de UCSC, generada con RepeatMasker y Repbase.

## 2. Planteamiento del problema

<!-- AYUDA: Describan la dificultad, necesidad o vacío de conocimiento que
desean atender, a quién afecta y qué sucede si no se resuelve. No confundan el
problema con la herramienta ni con la solución.
EJEMPLO: La identificación manual de genes de resistencia en varios genomas es
lenta, propensa a errores y difícil de reproducir. -->

La evidencia de que tau reactiva TEs proviene sobre todo de modelos animales y de análisis transcriptómicos o de tejido completo. No se sabe si la desrepresión de TEs deja una huella detectable en la accesibilidad de la cromatina humana, ni en qué tipos celulares ocurre, ni en qué punto de la progresión de la enfermedad aparece. Los datos de snATAC-seq de SEA-AD podrían responder esta pregunta, pero su tamaño (un objeto .h5ad de ~18 GB con cientos de miles de núcleos), la estructura multidonador del diseño y la necesidad de cruzarlos con una anotación de repeticiones externa vuelven difícil un análisis manual, correcto estadísticamente y reproducible.

## 3. Justificación

<!-- AYUDA: Expliquen por qué vale la pena realizar el proyecto, cuál es su
relevancia biológica, científica, técnica o social y quién podría beneficiarse.
EJEMPLO: Un flujo automatizado reducirá errores y permitirá repetir el análisis
con los mismos datos, parámetros y versiones del software. -->

Desde el punto de vista biológico, ubicar la posible apertura de TEs en subclases neuronales concretas y a lo largo del continuo de progresión ayudaría a evaluar si el mecanismo tau → TEs → Z-RNA–ZBP1 descrito en modelos tiene un correlato en tejido humano. También ayudaría a identificar qué poblaciones serían más vulnerables. Desde el punto de vista técnico, un flujo automatizado y versionado permite repetir el análisis con los mismos datos (versiones fechadas de SEA-AD y de rmsk), parámetros y dependencias (environment.yml). Además, maneja el volumen de datos en un servidor sin cargar la matriz completa en memoria. Pueden beneficiarse los grupos que estudian la epigenómica de la neurodegeneración y quienes quieran reutilizar el flujo con otras regiones cerebrales o modalidades de SEA-AD.

## 4. Objetivo general

<!-- AYUDA: Expresen el resultado global mediante un verbo en infinitivo. Debe
ser alcanzable durante el semestre.
EJEMPLO: Desarrollar un flujo reproducible en Python para identificar y comparar
genes de resistencia en un conjunto de genomas de E. coli. -->

Caracterizar los cambios de accesibilidad de la cromatina en poblaciones neuronales a lo largo de la progresión de la enfermedad de Alzheimer, y evaluar su relación con regiones que contienen elementos transponibles, como una aproximación indirecta a su posible reactivación.


## 5. Preguntas de investigación

<!-- AYUDA: Formulen preguntas biológicas o computacionales que puedan
responderse con los datos y métodos disponibles. Indiquen qué evidencia sería
necesaria.
EJEMPLO: ¿Qué genes de resistencia aparecen en cada genoma? Evidencia: una tabla
de presencia y ausencia obtenida de las anotaciones. -->

### Pregunta 1

**Pregunta:** [Escriban la pregunta.]  
**Evidencia necesaria:** [Datos o resultados que permitirán responderla.]

### Pregunta 2

**Pregunta:** [Escriban la pregunta.]  
**Evidencia necesaria:** [Datos o resultados que permitirán responderla.]

## 6. Alcance y limitaciones

<!-- AYUDA: Delimiten organismos, muestras, datos, análisis y resultado esperado. Indiquen
qué no se abordará y las restricciones de tiempo, cómputo, acceso o calidad.
EJEMPLO: Se analizarán como máximo 20 genomas completos de RefSeq. No se
utilizarán datos clínicos ni se realizará validación experimental. -->

### Incluye

- Datos humanos de snATAC-seq del giro temporal medio (MTG) de SEA-AD, en el objeto ya filtrado por los autores (final-nuclei, versión 2024-12-06).
- Los cuatro estadios de cambio neuropatológico de EA: Not AD, Low, Intermediate y High.
- Análisis detallado de tres subclases: L2/3 IT (excitadoras), Sst y Vip (inhibidoras). La agrupación y el pseudobulk se generan para todas las subclases.
- Anotación de TEs de la tabla rmsk de UCSC (hg38), restringida a cromosomas principales.
- Agregación pseudobulk por donador y subclase, y accesibilidad diferencial entre estadios.
- Figuras de control de calidad, composición celular y UMAP.]

### No incluye

- Otras regiones cerebrales de SEA-AD (p. ej. corteza prefrontal dorsolateral) ni otras modalidades (snRNA-seq, multioma, transcriptómica espacial).
- Llamado de picos de novo ni realineamiento de lecturas: se usan los picos y núcleos definidos por SEA-AD.
- Cuantificación de la expresión de TEs ni detección de inserciones polimórficas por donador.
- Validación experimental.

### Limitaciones conocidas

- La anotación rmsk de UCSC usa RepeatMasker 3.0.1 y Repbase de 2010, y no se ha recalculado. La clasificación de subfamilias jóvenes (L1HS, AluY, SVA) es menos precisa que con Dfam.
- Las copias recientes de TEs son casi idénticas entre sí; las lecturas cortas que caen en ellas suelen descartarse al alinear. Esos TEs pueden parecer cerrados aunque no lo estén.
- rmsk describe el genoma de referencia, no las inserciones propias de cada donador.
- Cada fila de rmsk es un fragmento, no necesariamente un elemento completo.
- La accesibilidad es un indicador indirecto de la reactivación de TEs, no una medida de su transcripción.
- Con el filtro de ≥10 núcleos por donador y subclase (MIN_NUCLEI), algunas combinaciones de subclase y estadio pueden quedar con pocas réplicas.
- El objeto .h5ad pesa ~18 GB. El flujo requiere un servidor (el paso pseudobulk está calibrado para ~64 GB de RAM).
- No hay pruebas automatizadas todavía.

## 7. Propuesta de solución

<!-- AYUDA: Describan el producto o estrategia que podría resolver el problema.
Es una propuesta inicial y puede cambiar. Expliquen sus componentes, no sólo las
tecnologías.
EJEMPLO: Un programa modular recibirá identificadores, descargará archivos,
extraerá genes, almacenará resultados y generará visualizaciones. -->

Un flujo modular de scripts numerados en src/ que se ejecutan en orden desde la raíz del repositorio. Cada uno se puede correr desde la terminal o importar como módulo mediante su función main(). Componentes:

1. Descarga (00_descargar_datos.py): obtiene el objeto .h5ad de snATAC-seq, el CSV de metadatos de todos los núcleos y la tabla rmsk, con versiones fijadas en el código.
2. Integración y validación de metadatos (01_analisis_metadatos.py): une los metadatos del .h5ad y del CSV, verifica la coincidencia de identificadores y la unicidad del perfil clínico por donador, y define los grupos de severidad. También genera una tabla por donador y exporta las coordenadas de los picos en formato BED.
3. Descripción de la cohorte (01b_graficos_datos_filtrados.py): grafica la composición de subclases por severidad y de supertipos por donador.
4. Partición de la matriz (02_group_atac_by_condition_and_subclass.py): escribe un .h5ad por combinación de estadio ADNC y subclase, en paralelo.
5. Agregación pseudobulk (03_pseudobulk_by_donor.py): suma los conteos por donador y subclase leyendo la matriz por bloques, con un límite de memoria. Genera las tablas de muestras (colData) para el análisis diferencial y figuras de control de calidad.

### 7.1 Resultado o producto esperado

<!-- AYUDA: Indiquen el entregable concreto: programa, paquete, flujo de análisis,
base de datos, visualizaciones u otro producto.
EJEMPLO: Repositorio ejecutable con scripts, datos de prueba, documentación,
tabla comparativa y figuras regenerables. -->

Un repositorio ejecutable con licencia MIT que incluya: los scripts del flujo; el entorno Conda (environment.yml); la documentación (README, este reporte, CITATION.cff); tablas de metadatos por núcleo y por donador; matrices pseudobulk por subclase listas para análisis diferencial; una tabla de picos anotados con TEs; resultados de accesibilidad diferencial; y figuras regenerables (composición celular, control de calidad del pseudobulk, UMAP).


## 8. Datos

### 8.1 Fuentes de datos

<!-- AYUDA: Incluyan institución, base de datos, URL, identificador, versión o
fecha de consulta y condiciones de uso. No todos los proyectos usarán NCBI.
EJEMPLO: NCBI RefSeq, GCF_000005845.2, consultado el dd/mm/aaaa. -->

| Fuente | Identificador o versión | URL | Fecha de consulta | Licencia o condiciones |
|:--|:--|:--|:--|:--|
| [SEA-AD (Allen Institute), bucket público de AWS S3 sea-ad-single-cell-profiling] | [SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad (versión 2024-12-06)] | [https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad] | [15/10/2026] | [Términos de uso del Allen Institute (https://alleninstitute.org/legal/terms-use/); citar a Gabitto et al. (2024)] |
| [SEA-AD (Allen Institute), bucket público de AWS S3 sea-ad-single-cell-profiling] | [SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv (versión 2024-12-06)] | [https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv] | [15/10/2026] | [Términos de uso del Allen Institute (https://alleninstitute.org/legal/terms-use/); citar a Gabitto et al. (2024)] |
| [UCSC Genome Browser, tabla rmsk de hg38] | [RepeatMasker 3.0.1 + Repbase 2010-03-02 (última actualización en UCSC: 18/10/2022)] | [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz] | [15/10/2026] | [Uso libre para fines académicos según UCSC; citar RepeatMasker y Repbase] |

### 8.2 Características de los datos

<!-- AYUDA: Describan organismos, muestras, variables, formatos, versiones, tamaño y otros
atributos necesarios para interpretar los datos.
EJEMPLO: Archivos FASTA y GFF3 de 20 genomas completos de E. coli. -->

- Organismo y tejido: Homo sapiens, giro temporal medio (MTG), tejido post mortem de donadores del atlas SEA-AD que abarcan todo el espectro de cambio neuropatológico de EA, incluidos donadores de referencia neurotípicos.
- snATAC-seq (.h5ad, ~18 GB): objeto AnnData. X es una matriz dispersa de conteos núcleos × picos; obs contiene los metadatos por núcleo (donador, anotación celular jerárquica, variables neuropatológicas y clínicas, métricas de calidad); var contiene los picos, cuyo identificador codifica sus coordenadas. Los núcleos ya fueron filtrados por los autores (final-nuclei). Coordenadas en GRCh38, con nombres de cromosoma chr1, chr2…
- Metadatos (.csv): una fila por núcleo, incluidos los que eliminó el control de calidad de SEA-AD (all-nuclei). En el flujo solo se conservan los núcleos presentes en el .h5ad.
- RepeatMasker (rmsk.txt.gz, ~152 MB comprimido): 5,683,690 repeticiones en 17 columnas separadas por tabuladores. Coordenadas en base 0, semiabiertas, sin ordenar. La versión de hgdownload no tiene encabezado. Con las clases LINE, SINE, LTR, DNA y Retroposon en chr1–chr22, chrX y chrY quedan 4,505,027 fragmentos de TEs.


### 8.3 Organización de los datos

<!-- AYUDA: Muestren la estructura prevista. No suban datos sensibles, tokens,
contraseñas ni archivos grandes. Usen .gitignore y documenten cómo obtener lo
que no se guarde en Git.
EJEMPLO: data/raw conserva originales y data/processed los derivados. -->

data/, logs/ y results/ están en .gitignore: los datos y resultados se generan localmente con los scripts y no se versionan.

```text
NeuroEpigenomics_AD/
├── data/                                   # No versionado
│   ├── SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad    # 00
│   ├── rmsk.txt.gz                                       # 00
│   ├── metadata/
│   │   └── SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv  # 00
│   └── processed/
│       ├── <Condición>_<Subclase>.h5ad                   # 02
│       ├── metadata_test.csv, peaks.txt                  # 02
│       └── pseudobulk/
│           ├── pseudobulk_<Subclase>.h5ad                # 03
│           ├── colData_<Subclase>.csv                    # 03
│           └── pseudobulk_summary.csv                    # 03
├── docs/
│   ├── reporte-proyecto.md
│   └── CITATION.cff
├── logs/                                   # No versionado (salida de nohup)
├── notebooks/
├── results/                                # No versionado
│   ├── 01_metadata_output/                               # 01
│   ├── 01b_graficas/                                     # 01b
│   └── figures/                                          # 03 (pseudobulk/) y 04
├── src/                                    # Scripts 00–04
├── tests/
├── environment.yml
├── LICENSE
└── README.md
```

### 8.4 Diccionario o formato de los datos

<!-- AYUDA: Describan campos o columnas relevantes. Incluyan fragmentos pequeños
cuando ayuden a comprender el formato, pero no archivos completos.
EJEMPLO: En GFF3, seqid identifica la secuencia; type indica gene, CDS, etc. -->

| Archivo o conjunto | Campo/columna | Tipo | Descripción | Valores o unidades |
|:--|:--|:--|:--|:--|
| [Archivo] | [Campo] | [Tipo] | [Descripción] | [Valores] |

## 9. Metodología

<!-- AYUDA: Esta sección evolucionará. Primero describan el plan y después
actualícenla con lo que realmente ejecutaron, incluidos parámetros y decisiones. -->

### 9.1 Etapas del análisis o desarrollo

<!-- AYUDA: Describan la secuencia desde la obtención de datos hasta la validación
de resultados. Relacionen cada etapa con una pregunta u objetivo.
EJEMPLO: descarga, validación, transformación, análisis, visualización y pruebas. -->

1. Obtención de datos (00_descargar_datos.py). Descarga con wget -c (o con urllib si wget no está instalado) y omite los archivos que ya existen. Objetivo general.
2. Integración y validación de metadatos (01_analisis_metadatos.py). Lee el .h5ad en modo backed='r'; comprueba la intersección de identificadores con el CSV; une los metadatos priorizando la anotación del .h5ad; define Severity_Group; clasifica las métricas de calidad como de núcleo o de biblioteca; construye una tabla por donador (falla si un donador tiene más de un perfil clínico); y exporta atac_peaks.bed. 
3. Descripción de la cohorte (01b_graficos_datos_filtrados.py). Genera barras apiladas de núcleos por subclase y severidad. También grafica la proporción de supertipos por donador ordenada por CPS para L2/3 IT, Sst y Vip: los supertipos con menos de 3% se agrupan en "Otros" y se marcan los donadores con menos de 50 núcleos. 
4. Partición por estadio y subclase (02_group_atac_by_condition_and_subclass.py). Escribe un .h5ad comprimido por combinación ADNC × subclase, con 6 procesos en paralelo (joblib). 
5. Agregación pseudobulk (03_pseudobulk_by_donor.py). Suma los conteos por (subclase, donador) leyendo la matriz en bloques de 20,000 núcleos, mediante una matriz indicadora dispersa. Agrupa las subclases en lotes para que el acumulador no supere 8 GB y descarta los donadores con menos de 10 núcleos. Exporta las matrices y colData, más tres figuras de control de calidad. 
6. Exploración no supervisada (04_umap_atac_by_subclass.py). Por subclase y en conjunto: elimina núcleos y picos sin conteos, selecciona 50,000 picos (snap.pp.select_features), calcula el embedding espectral (snap.tl.spectral), el grafo de 30 vecinos y el UMAP (semilla 2). Colorea por conteos totales, supertipo, subclase y condición. 

### 9.2 Herramientas y tecnologías

<!-- AYUDA: Registren lenguajes, bibliotecas y programas con sus versiones y
propósito. No incluyan credenciales.
EJEMPLO: Python 3.x; Biopython para leer formatos biológicos; Seaborn para
visualización. -->

| Herramienta | Versión | Propósito |
|:--|:--|:--|
| [Herramienta] | [Versión] | [Uso] |


### 9.3 Estrategia de validación

<!-- AYUDA: Expliquen cómo comprobarán código y resultados: pruebas unitarias,
datos conocidos, comparación con otra herramienta o revisión manual.
EJEMPLO: Se compararán cinco anotaciones conocidas y se probarán entradas
válidas, identificadores inexistentes y archivos incompletos. -->

[Describan aquí la validación.]

## 10. Plan de trabajo


### 10.1 Distribución de responsabilidades

<!-- AYUDA: Definan responsabilidades iniciales sin aislar a cada integrante.
Toda contribución importante debe ser revisada mediante Pull Request por otra
persona.
EJEMPLO: Ana desarrolla la descarga y revisa el módulo de visualización. -->

| Integrante | Responsabilidad principal | Responsabilidad de revisión |
|:--|:--|:--|
| [Nombre] | [Responsabilidad] | [Qué o a quién revisará] |

### 10.2 Riesgos y alternativas

<!-- AYUDA: Identifiquen situaciones que podrían impedir o retrasar el proyecto
y definan una alternativa.
EJEMPLO: Los datos requieren demasiado almacenamiento; alternativa: reducir el
número de genomas usando criterios documentados. -->

| Riesgo | Probabilidad | Impacto | Prevención o alternativa |
|:--|:--|:--|:--|
| [Riesgo] | Baja/Media/Alta | Bajo/Medio/Alto | [Acción] |

## 11. Resultados

<!-- AYUDA: Presenten resultados vinculados con preguntas y objetivos. Incluyan
tablas o figuras con títulos, leyendas y archivos de origen. Describan aquí lo
obtenido; interprétenlo en Discusión.
EJEMPLO: Tabla de presencia y ausencia generada por src/compare_genes.py.
PRIMERA SESIÓN: dejen esta sección como pendiente. -->

> Estado: pendiente. Se completará durante el desarrollo.



## 12. Discusión

<!-- AYUDA: Interpreten los resultados, expliquen si responden las preguntas,
compárenlos con los antecedentes y señalen limitaciones. No repitan únicamente
los valores.
EJEMPLO: La distribución observada sugiere..., aunque la interpretación está
limitada por la calidad de las anotaciones. -->

> Estado: pendiente. Se completará después de obtener resultados.

## 13. Conclusiones

<!-- AYUDA: Sinteticen qué se aprendió, qué preguntas se respondieron y si se
alcanzaron los objetivos. Incluyan aportes, limitaciones y trabajo futuro. No
introduzcan resultados nuevos.
EJEMPLO: El flujo permitió identificar..., pero será necesario incorporar... -->

> Estado: pendiente. Se completará al finalizar el proyecto.


## 14. Disponibilidad, licencia y citación

<!-- AYUDA: Indiquen dónde está el código, bajo qué licencia puede reutilizarse
y cómo citarlo. Relacionen esta sección con LICENSE, CITATION.cff, codemeta.json,
release final y, cuando corresponda, un DOI.
EJEMPLO: Código en GitHub bajo MIT; cita disponible en CITATION.cff. -->

**Código:** [https://github.com/rafadmtz/NeuroEpigenomics_AD]  
**Datos:** [No se redistribuyen en el repositorio. Se obtienen con python src/00_descargar_datos.py desde el bucket público de SEA-AD en AWS S3 (sea-ad-single-cell-profiling, versión 2024-12-06) y desde UCSC (rmsk de hg38); ver sección 8.1.]  
**Licencia del código:** [MIT (© 2026 rafadmtz; ver LICENSE). Los datos de SEA-AD y de UCSC/RepeatMasker conservan sus propias condiciones de uso.]  
**Cómo citar:** [docs/CITATION.cff (CFF 1.2.0; fecha de liberación 19/10/2026):
Díaz Martínez, R., Martínez Enciso, S. y González Morales, S. N. (2026). NeuroEpigenomics_AD: Epigenomic dynamics during Alzheimer's disease progression. Repositorio de GitHub. https://github.com/rafadmtz/NeuroEpigenomics_AD
]  
**Versión o release:** [URL]

## 15. Referencias

<!-- AYUDA: Registren publicaciones, datos, software y documentos consultados en
un formato uniforme. Incluyan DOI, URL o identificadores persistentes. Toda cita
del texto debe aparecer aquí.
EJEMPLO: Blattner, F. R. et al. (1997). The complete genome sequence of
Escherichia coli K-12. Science, 277(5331), 1453–1462.
https://doi.org/10.1126/science.277.5331.1453 -->

1. [Referencia 1]
2. [Referencia 2]


---

<!-- ORIENTACIÓN PARA LAS DOS PRIMERAS SESIONES:
Completen Información general, Resumen provisional, secciones 1 a 7, fuentes de
datos preliminares y plan de trabajo. Metodología, Resultados, Discusión,
Conclusiones Y Disponibilidad evolucionarán durante el
semestre. Sustituyan las indicaciones entre corchetes por contenido del equipo. -->