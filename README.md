### README: Análisis y Modificación de Árboles Filogenéticos

Este documento presenta un conjunto de ejercicios prácticos centrados en el análisis, modificación y generación de árboles filogenéticos utilizando herramientas bioinformáticas. A continuación, se describen los objetivos, métodos y resultados obtenidos en cada ejercicio.

---

## **Ejercicio 1: Análisis de un Árbol Filogenético Existente**

### Descripción:
Se analizó un árbol filogenético en formato **Newick**, construido a partir de secuencias de la subunidad beta 1 de la hemoglobina. Este árbol refleja las relaciones evolutivas entre especies de peces, aves y mamíferos.

### Detalles del análisis:
- **Número de terminales:** 22 (especies individuales).
- **Número de nodos o clados:** 42.
- **Longitudes de ramas:** Varían entre 0.0 y 1.13866, reflejando diferentes niveles de divergencia evolutiva.
- **Especies más cercanas:** KX241132.1 y KX241142.1 (distancia 0.0).
- **Especie más distante:** MZ593243.1, con una distancia acumulada de 29.5649.

### Visualización:
El árbol permite observar agrupamientos claros entre linajes específicos y divergencias evolutivas significativas. Se identificaron clados compactos y linajes altamente divergentes.

### Código principal:
- **Análisis:** Determinación de terminales, nodos, distancias y patrones de agrupamiento.
- **Visualización:** Personalización del árbol para resaltar relaciones clave.

---

## **Ejercicio 2: Modificación de un Árbol Filogenético**

### Descripción:
Se modificó el árbol filogenético del ejercicio anterior para alterar nombres de especies, longitudes de ramas y colores de algunas conexiones. Estas modificaciones destacan relaciones clave entre especies.

### Cambios realizados:
1. **Cambio de nombres:**
   - KX241132.1 → *Eriocnemis luciani*.
   - KX241142.1 → *Haplophaedia aureliae*.
   - MZ593243.1 → *Arenicola marina*.
   
2. **Ajuste de longitudes:**
   - Las ramas que conectan *Eriocnemis luciani* y *Haplophaedia aureliae* se establecieron en 0.5.
   - La rama de *Arenicola marina* se ajustó a 0.9.

3. **Colores asignados:**
   - Azul para *Eriocnemis luciani* y *Haplophaedia aureliae*.
   - Rojo para *Arenicola marina*.

### Resultados:
El árbol modificado resalta visualmente las especies más cercanas y la más distante, facilitando su interpretación.

### Código principal:
- Funciones para cambiar nombres, ajustar longitudes y asignar colores.
- Guardado del árbol modificado en formato Newick para futuras referencias.

---

## **Ejercicio 3: Generación de Árboles Filogenéticos desde Secuencias Aleatorias**

### Descripción:
Se generaron árboles filogenéticos a partir de secuencias de nucleótidos creadas de manera aleatoria. Estas secuencias se almacenaron en un archivo **FASTA** para su posterior análisis.

### Pasos realizados:
1. **Generación de secuencias:**
   - Se crearon 5 secuencias de 5 nucleótidos de forma aleatoria.
   - Se almacenaron en el archivo `./fasta_files/secuencias.fasta`.

2. **Métodos empleados:**
   - **Máxima parsimonia:** Construcción de un árbol filogenético basado en la minimización de eventos evolutivos.
   - **Métodos de distancia:**
     - **UPGMA:** Generación de árboles basados en distancias promedio entre grupos.
     - **Neighbor Joining (NJ):** Método de agrupamiento progresivo para generar árboles filogenéticos.

3. **Cálculo de distancias:**
   - Distancias calculadas con medidas como "identity", "blastn", y "trans".

### Comparación de resultados:
- Los árboles generados mediante máxima parsimonia y métodos de distancia mostraron patrones similares en agrupamientos cercanos, pero diferencias en ramas más largas.
- UPGMA mostró una estructura más balanceada, mientras que NJ destacó mejor las diferencias evolutivas.

### Código principal:
- Generación de secuencias y archivo FASTA.
- Construcción y visualización de árboles mediante parsimonia y métodos de distancia.

---

## **Conclusiones Generales**

1. El análisis del árbol filogenético inicial permitió identificar patrones de divergencia y similitud evolutiva, destacando linajes claramente diferenciados y agrupamientos estrechos.
2. Las modificaciones realizadas en el árbol resaltaron visualmente las relaciones clave, facilitando su interpretación y análisis.
3. Los árboles generados a partir de secuencias aleatorias permitieron explorar los resultados de diferentes métodos de construcción (parsimonia y distancia). Los métodos mostraron fortalezas en diferentes aspectos del análisis evolutivo, destacando la importancia de elegir el enfoque adecuado según el contexto del estudio.

---

## **Archivos Generados**

1. **Árbol inicial:** `hbb_tree.txt`.
2. **Árbol modificado:** `modified_hbb_tree.txt`.
3. **Secuencias FASTA:** `./fasta_files/secuencias.fasta`.
4. **Resultados visuales:** Gráficos de los árboles generados para cada método.

Este conjunto de ejercicios ilustra cómo analizar, modificar y generar árboles filogenéticos, destacando las herramientas bioinformáticas y su aplicación en estudios evolutivos.
