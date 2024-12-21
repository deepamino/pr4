# Análisis y Modificación de Árboles Filogenéticos

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

```python
def analyze_tree(newick_file):
    tree = Phylo.read(newick_file, "newick")
    
    num_terminals = len(tree.get_terminals())
    num_clades = len(list(tree.find_clades()))
    branch_lengths = [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]  # Longitudes de las ramas
    terminal_names = [terminal.name for terminal in tree.get_terminals()]
    
    clades_with_bootstrap = [(clade.confidence, clade) for clade in tree.find_clades() if clade.confidence is not None]
    best_bootstrap, best_clade = (None, None)
    if clades_with_bootstrap:
        best_bootstrap, best_clade = max(clades_with_bootstrap, key=lambda x: x[0])
    
    min_distance, closest_pair = float('inf'), None
    for (name1, name2) in combinations(terminal_names, 2):
        distance = tree.distance(name1, name2)
        if distance < min_distance:
            min_distance = distance
            closest_pair = (name1, name2)
    
    max_distance, most_distant_species = float('-inf'), None
    for name in terminal_names:
        total_distance = sum(tree.distance(name, other) for other in terminal_names if other != name)
        if total_distance > max_distance:
            max_distance = total_distance
            most_distant_species = name
    
    print(f"Información del árbol filogenético:")
    print(f"Número de terminales: {num_terminals}")
    print(f"Número de nodos/clados totales: {num_clades}")
    print(f"Longitudes de las ramas: {branch_lengths}")
    print(f"Nombres de terminales: {terminal_names}")

    if best_clade:
        print(f"Clado con el mejor bootstrap: {best_bootstrap}")

    if closest_pair:
        print(f"Especies más parecidas: {closest_pair} con distancia {min_distance}")
        
    if most_distant_species:
        print(f"Especie más distante al resto: {most_distant_species} con distancia total {max_distance}")
    
    return tree, closest_pair, most_distant_species
```

El resultado de ejecutar esta función sobre el árbol generado a partir de la hemoglobina fue el siguiente:

```python
newick_file = "./data/hbb_tree.txt"
tree, closest_pair, most_distant_species = analyze_tree(newick_file)
```

```text
Información del árbol filogenético:
Número de terminales: 22
Número de nodos/clados totales: 42
Longitudes de las ramas: [0.205701, 1.13866, 0.031261, 0.136879, 0.133087, 0.230584, 0.190835, 0.0953875, 0.162102, 0.28055, 0.0278707, 0.0186001, 0.0572805, 0.241457, 0.0286879, 0.0473548, 0.0106191, 0.00200754, 0.0021201, 0.00865119, 0.00195866, 0.0066246, 0.00852267, 5e-09, 0.0, 0.0, 0.0230908, 0.0215049, 0.00236303, 0.0103989, 0.012181, 0.0269614, 0.00456229, 0.00824827, 0.0273313, 0.00829855, 0.00242122, 0.023679, 0.00222634, 0.00210352, 0.00853381]
Nombres de terminales: ['AB364477.1', 'MZ593243.1', 'BT074827.1', 'BT082972.1', 'MT164172.1', 'NM_001160555.2', 'BT059665.1', 'OL804561.1', 'KX241110.1', 'KX241171.1', 'KX241173.1', 'KX241147.1', 'KX241132.1', 'KX241142.1', 'KX241189.1', 'KX241204.1', 'KX241222.1', 'KX241216.1', 'KX241252.1', 'KX241297.1', 'KX241301.1', 'KX241309.1']
Especies más parecidas: ('KX241132.1', 'KX241142.1') con distancia 0.0
Especie más distante al resto: MZ593243.1 con distancia total 29.56490544
```

- **Visualización:** Personalización del árbol para resaltar relaciones clave.

```python
def plot_tree(tree, color_dict=None, branch_width=2):
    if color_dict:
        for clade in tree.find_clades():
            if clade.name in color_dict:
                clade.color = color_dict[clade.name]

    def color_func(clade):
        return getattr(clade, "color", "black")

    plt.figure(figsize=(15, 8))
    ax = plt.gca()

    Phylo.draw(tree, axes=ax, branch_labels=None, label_colors=color_func, do_show=False)

    for line in ax.findobj(match=lambda obj: isinstance(obj, plt.Line2D)):
        line.set_linewidth(branch_width)

    plt.title("Árbol Filogenético de la Subunidad Beta de la Hemoglobina")
    plt.show()
```

Esta función mostrará el árbol filogenético mediante el uso de la librería `matplotlib`, cosa que es autogestionada por biopython directamente, destacando en color rojo la especie más alejada genéticamente del resto. El resultado de esta función para el árbol de estudio se puede ver en la Figura 1.

<div align="center">
    <img src="resources/ej1_tree.png" width="70%">
    <p><b>Figura 1.</b> Árbol filogenético generado con la subunidad beta de la hemoglobina en distintos organismos.</p>
</div>

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

```python
def modify_tree(tree):
    for clade in tree.find_clades():
        if clade.name == 'KX241132.1':
            clade.name = 'Eriocnemis luciani (KX241132.1)'

        elif clade.name == 'KX241142.1':
            clade.name = 'Haplophaedia aureliae (KX241142.1))'

        elif clade.name == 'MZ593243.1':
            clade.name = 'Arenicola marina (MZ593243.1)'

    for clade in tree.find_clades():
        if clade.name in ['Eriocnemis luciani (KX241132.1)', 'Haplophaedia aureliae (KX241142.1))']:
            clade.branch_length = 0.5
            clade.color = 'blue'

        elif clade.name == 'Arenicola marina (MZ593243.1)':
            clade.branch_length = 0.9
            clade.color = 'red'

    return tree
```

- Función para visualizar el árbol, similar a la del ejercicio anterior aunque con algunas modificaciones.
- Guardado del árbol modificado en formato Newick para futuras referencias.

 El resultado de este ejercicio se puede ver en la Figura 2, donde se han destacado las especies mencionadas tal y como se ha descrito.

<div align="center">
    <img src="resources/ej2_tree.png" width="70%">
    <p><b>Figura 2.</b> Árbol filogenético modificado, destacando tanto la rama más lejana como las más cercancas entre sí.</p>
</div>

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

```python
def generar_secuencia(longitud):
    nucleotidos = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotidos) for _ in range(longitud))
```

```python
def crear_fasta(nombre_archivo, num_secuencias, longitud):
    with open(nombre_archivo, 'w') as archivo:
        for i in range(1, num_secuencias + 1):
            secuencia = generar_secuencia(longitud)
            archivo.write(f">Secuencia_{i}\n")  # Encabezado FASTA
            archivo.write(f"{secuencia}\n")
```

```python
nombre_archivo = "./fasta_files/secuencias.fasta"
num_secuencias = 5
longitud_secuencia = 5

crear_fasta(nombre_archivo, num_secuencias, longitud_secuencia)

print(f"Archivo {nombre_archivo} generado con {num_secuencias} secuencias de {longitud_secuencia} nucleótidos.")
```

- Construcción y visualización de árboles mediante parsimonia y métodos de distancia.

Algunos de los resultados obtenidoss han sido los siguientes:

En cuanto al método de máxima parsimonia, el resultado se muestra en la FIgura 3.

```python
generator = GeneratorTreeFactory().initialize_generator('ParsimoniaTree')
tree_parsimonia = generator.generate_tree(sequences)
generator.show_tree(tree_parsimonia)
```

<div align="center">
    <img src="resources/ej3_par_tree.png" width="70%">
    <p><b>Figura 3.</b> Árbol filogenético por el método de máxima parsimonia.</p>
</div>

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
