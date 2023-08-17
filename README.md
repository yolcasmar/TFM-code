# Código
Los scripts se encuentran en el directorio *Codigo* numerados según el orden de aparición en la metodología del TFM.

Siendo el 01 el código empleado en el apartado 3.4.1, el 02 en el apartado 3.4.2, el 03 en el 3.4.3, el 04 en 3l 3.4.4 y el 05 en el 3.4.5.

El script *my_generate_priorization_tables_function.R* contiene la modificación de la función original del paquete MultiNicheNet y se usa en el script 05.
# Figuras R
Las figuras que aparecen en los resultados del TFM se encuentran en el directorio *Figuras R*.

Dentro se distribuyen en subdirectorios dependiendo del script que se ha empleado para obtenerlas. Dentro de estos subdirectorios se encuentran más subdirectorios para facilitar la búsqueda de cada Figura.

El directorio *decoupler_collectri* contiene a su vez: *ComplexHeatmap*, con los dos mapas de calor mencionados en el apartado 4.1. del trabajo; *expresion_actividad*, que contiene el UMAP de los tipos celulares, los UMAPs de la expresión y actividad de PGR y ESR1 y los gráficos de violín de la expresión y actividad de PGR y ESR1 según el tipo celular y condición de estudio (también en el apartado 4.1.); *pseudobulk* corresponde a los resultados del apartado 4.2., contiene la red de intracción de los TFs PGR y ESR1 y un subdirectorio con los gráficos de volcán (ordenados en más subdirectorios según el tipo celular).

El directorio *multinichenet* contiene los gráficos circos mencionados en el apartado 4.5. del trabajo y la red de interacción ligando-diana, del apartado 4.6.

El directorio *nichenet_niches*, que corresponde al apartado 4.4., contiene: el subdirectorio *epi-stroma* con los gráficos de los LFC cuando los epitelios se consideran células emisoras y el estroma se considera células receptoras; y el subdirectorio *stroma-epi* con los gráficos de los LFC cuando el estroma se considera emisor y los epitelios receptores.

El directorio *nichenet_seurat*, que corresponde al apartado 4.3., contiene: el subdirectorio *seurat epi-stroma* con el gráfico combinado y el mapa de calor ligando-receptor mencionados en el trabajo (considerando los epitelios como emisores y el estroma como receptor); y el subdirectorio *seurat stroma-epi* con los mismos gráficos pero considerando el estroma como emisor y los epitelios como receptores. 

# Code
Scripts are in *Codigo* directory numbered according to the order of appearance in the TFM methodology.

Being the 01 the code employed in section 3.4.1, the 02 in section 3.4.2, the 03, in 3.4.3, the 04 in 3.4.4 and the 05 in 3.4.5.

The script *my_generate_priorization_tables_function.R* contains the original function modification from the MultiNicheNet package and is used in script 05.
# R Figures
The figures that appear in the results of the TFM are found in the *Figuras R* directory.

They are distributed inside in subdirectories depending on the script that has been used to obtain them. Within these subdirectories there are more subdirectories to facilitate the search for each Figure.

The *decoupler_collectri* directory contains in turn: *ComplexHeatmap*, with the two heatmaps mentioned in section 4.1. from work; *expresion_actividad*, which contains the cell type UMAP, UMAPs of PGR and ESR1 expression and activity and violin plots of PGR and ESR1 expression and activity according to cell type and study condition (also in section 4.1.); *pseudobulk* corresponds to the results of section 4.2., contains the interaction network of the PGR and ESR1 TFs and a subdirectory with the volcano plots (sorted into more subdirectories according the cell type).

The *multinichenet* directory contains the circos plots mentioned in section 4.5. of work and the ligand-target interaction network, from section 4.6.

The *nichenet_niches* directory, which corresponds to section 4.4., contains: the *epi-stroma* subdirectory with the LFC graphs when epithelia are considered sender cells and the stroma is considered receptor cell; and the *stroma-epi* subdirectory with the LFC graphs when the stoma is considered sender and the epithelia are considered receptors.

The *nichenet_seurat* directory, which corresponds to section 4.3., contains: the *seurat epi-stroma* subdirectory with the combined plot and the ligand-receptor heatmap mentioned in the work (considering the epithelia as senders and the stroma as receptor); and the *seurat stroma-epi* subdirectory with the same graphs but considering the stroma as sender and the epithelia as receptors.
