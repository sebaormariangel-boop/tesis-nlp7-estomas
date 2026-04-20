# Meta-análisis de transcriptómica a nivel de célula única en hojas de *Arabidopsis thaliana* para identificar circuitos regulatorios de señalización por sequía y nitrógeno

Código y tablas suplementarias externas asociadas a mi tesis de pregrado sobre transcriptómica de células de guarda en *Arabidopsis thaliana*, con énfasis en la integración de datasets single-cell y single-nucleus, análisis de coexpresión y comparación con datos bulk RNA-seq.

## Resumen

Este repositorio reúne scripts, tablas suplementarias externas y archivos auxiliares utilizados en el desarrollo de la tesis. El trabajo se centra en la identificación de programas transcripcionales y circuitos regulatorios asociados a señalización por sequía, ABA y nitrógeno, con especial atención a células de guarda.

Los análisis integran:

- datos públicos de transcriptómica single-cell y single-nucleus de hoja de *Arabidopsis thaliana*,
- datos bulk RNA-seq de hojas completas y de fracciones estomáticas/no estomáticas separadas por FACS,
- análisis de expresión diferencial,
- anotación celular,
- proyección de firmas transcriptómicas,
- análisis de coexpresión mediante hdWGCNA,
- comparación entre módulos y conjuntos génicos.

## Contenido del repositorio

- `scripts/`: scripts de preprocesamiento, análisis e integración de datos.
- `supplementary_tables/`: tablas suplementarias externas en formato `.csv` o `.xlsx`.

## Tablas suplementarias externas

Las tablas de gran tamaño que no fueron incorporadas directamente en el documento principal o en el PDF de material suplementario se incluyen aquí como archivos suplementarios externos. En el manuscrito, estos archivos se citan como **Tabla Suplementaria X** según corresponda.

Ejemplos:

- **Tabla Suplementaria 6.** 
- **Tabla Suplementaria 7.**
- **Tabla Suplementaria 8.** 
- **Tabla Suplementaria 9.**

## Descripción general del análisis

De manera resumida, el flujo de trabajo incluyó:

1. recopilación y preprocesamiento de datasets públicos single-cell y single-nucleus;
2. control de calidad, remoción de dobletes y evaluación de firma de digestión;
3. integración, reducción dimensional, clustering y anotación celular;
4. validación de anotaciones mediante firmas transcriptómicas independientes;
5. construcción de redes de coexpresión con hdWGCNA;
6. análisis funcional de módulos y genes hub;
7. comparación con resultados de bulk RNA-seq y enriquecimiento entre conjuntos de genes.

## Reproducibilidad

Los análisis fueron realizados principalmente en R. Las versiones exactas de paquetes, parámetros y criterios de filtrado se describen en el manuscrito.

## Correspondencia con el manuscrito

Este repositorio acompaña la tesis titulada:

**Meta-análisis de transcriptómica a nivel de célula única en hojas de *Arabidopsis thaliana* para identificar circuitos regulatorios de señalización por sequía y nitrógeno**

Las tablas y análisis aquí contenidos corresponden a la versión digital del manuscrito y su material suplementario.

## Autor

Sebastián Ortiz  
Ingeniería en Biotecnología Molecular  
Universidad de Chile

## Nota

Este repositorio fue creado con fines académicos para respaldar la reproducibilidad, trazabilidad y consulta de resultados suplementarios asociados a la tesis.
