
# Cargar paquetes
## Se cargan los paquetes requeridos
library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(decoupleR)
library(OmnipathR)

# Solo necesarios para el tratamiento de datos y las gráficas
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ComplexHeatmap)
library(ggraph)
library(tidygraph)

# Cargar el set de datos, reducirlo para que su tamaño no sea muy grande, normalizar los datos y encontrar los 1000 genes más variables
seuratObj = LoadH5Seurat("../Data/Single-cell/07_sc_AS_Estonia_Wang_Group_SCT_Main_lineages_clean.H5Seurat")
seuratObj <- DietSeurat(object = seuratObj, assays = c("RNA"), dimreducs = "umap")
seuratObj$Group_final <- gsub(pattern = "WOI WOI Control", replacement = "WOI", x = seuratObj$Group_final)
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj, nfeatures = 1000)

# Hacer un umap de los tipos celulares
cell_type_umap <- DimPlot(seuratObj, group.by = "cell_type_zooming_v5", reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("Tipos celulares")
ggsave(filename = "cell_type_umap", path = "../Results/Images/decoupleR_CollecTRI_images/", plot = cell_type_umap, device = png)

# Cargar la colección de regulones de CollecTRI
net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
saveRDS(net, file = "../Data/decoupleR_CollecTRI_data/net.rds")

# Inferencia de actividad con ULM
# Extraer los conteos normalizados transformados por logaritmo
mat <- seuratObj@assays$RNA@data
saveRDS(mat, file = "../Data/decoupleR_CollecTRI_data/mat.rds")

# Calcular ULM
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)
saveRDS(acts, file = "../Data/decoupleR_CollecTRI_data/acts.rds")

# Visualización
# Extraer los ulm y almacenarlos en el ensayo tfsulm del objeto Seurat
seuratObj[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Cambiar el ensayo por defecto
DefaultAssay(object = seuratObj) <- "tfsulm"

# Escalar los datos
seuratObj <- ScaleData(seuratObj)
seuratObj@assays$tfsulm@data <- seuratObj@assays$tfsulm@scale.data

# Este nuevo ensayo puede usarse para graficar las actividades, mientras que el ensayo de RNA se usa para graficar la expresión
# Se grafica para los TFs PGR y ESR1
p1 <- DimPlot(seuratObj, group.by = "cell_type_zooming_v5", reduction = "umap", label = TRUE, pt.size = 0.5) +
  NoLegend() + ggtitle('Tipos celulares')
p2 <- (FeaturePlot(seuratObj, raster=FALSE, features = c("PGR")) &
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Actividad de PGR')
DefaultAssay(object = seuratObj) <- "RNA"
p3 <- FeaturePlot(seuratObj, raster=FALSE, features = c("PGR")) + ggtitle('Expresión de PGR')
DefaultAssay(object = seuratObj) <- "tfsulm"
p1 | p2 | p3

ggsave(filename = "DimPlotPGR", path = "../Results/Images/decoupleR_CollecTRI_images/", plot = p1, device = png)
ggsave(filename = "ActivityPGR", path = "../Results/Images/decoupleR_CollecTRI_images/", plot = p2, device = png)
ggsave(filename = "ExpressionPGR", path = "../Results/Images/decoupleR_CollecTRI_images/", plot = p3, device = png)

p4 <- DimPlot(seuratObj, group.by = "cell_type_zooming_v5", reduction = "umap", label = TRUE, pt.size = 0.5) +
  NoLegend() + ggtitle('Tipos celulares')
p5 <- (FeaturePlot(seuratObj, raster=FALSE, features = c("ESR1")) &
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Actividad de ESR1')
DefaultAssay(object = seuratObj) <- "RNA"
p6 <- FeaturePlot(seuratObj, raster=FALSE, features = c("ESR1")) + ggtitle('Expresión de ESR1')
DefaultAssay(object = seuratObj) <- "tfsulm"
p4|p5|p6

ggsave(filename = "DimPlotESR1", path = "../Results/Images/decoupleR_CollecTRI_images/", plot = p4, device = png)
ggsave(filename = "ActivityESR1", path = "../Results/Images/decoupleR_CollecTRI_images/", plot = p5, device = png)
ggsave(filename = "ExpressionESR1", path = "../Results/Images/decoupleR_CollecTRI_images/", plot = p6, device = png)

## Violin plots
# Calcular los Violin plots de la actividad de los TFs de interés
seuratObj <- SetIdent(seuratObj, value = "cell_type_zooming_v5")

comparisons = list(c("AS", "WOI"))
y_max = 5

Vln_Stroma_PGR <- VlnPlot(seuratObj, features = "PGR",
        pt.size = 0, 
        group.by = "Group_final", 
        idents = "Stroma",
        raster = FALSE,
        y.max = y_max) + 
  xlab("Estroma") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Lum_PGR <- VlnPlot(seuratObj, features = "PGR",
        pt.size = 0, 
        group.by = "Group_final", 
        idents = "Epi-Lumenal",
        raster = FALSE,
        y.max = y_max) + 
  xlab("Epitelio Lumenal") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Cil_PGR <- VlnPlot(seuratObj, features = "PGR",
        pt.size = 0, 
        group.by = "Group_final", 
        idents = "Epi-Ciliated",
        raster = FALSE,
        y.max = y_max) + 
  xlab("Epitelio Ciliado") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Glan_PGR <- VlnPlot(seuratObj, features = "PGR",
        pt.size = 0, 
        group.by = "Group_final", 
        idents = "Epi-Glandular",
        raster = FALSE,
        y.max = y_max) + 
  xlab("Epitelio Glandular") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_GlanSec_PGR <- VlnPlot(seuratObj, features = "PGR",
        pt.size = 0, 
        group.by = "Group_final", 
        idents = "Epi-Glandular secretory",
        raster = FALSE,
        y.max = y_max) + 
  xlab("Epitelio Glandular secretor") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Stroma_PGR | Vln_Lum_PGR | Vln_Cil_PGR | Vln_Glan_PGR | Vln_GlanSec_PGR

y_max = 4
Vln_Stroma_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                          pt.size = 0, 
                          group.by = "Group_final", 
                          idents = "Stroma",
                          raster = FALSE,
                          y.max = y_max) + 
  xlab("Estroma") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Lum_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                       pt.size = 0, 
                       group.by = "Group_final", 
                       idents = "Epi-Lumenal",
                       raster = FALSE,
                       y.max = y_max) + 
  xlab("Epitelio Lumenal") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Cil_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                       pt.size = 0, 
                       group.by = "Group_final", 
                       idents = "Epi-Ciliated",
                       raster = FALSE,
                       y.max = y_max) + 
  xlab("Epitelio Ciliado") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Glan_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                        pt.size = 0, 
                        group.by = "Group_final", 
                        idents = "Epi-Glandular",
                        raster = FALSE,
                        y.max = y_max) + 
  xlab("Epitelio Glandular") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_GlanSec_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                           pt.size = 0, 
                           group.by = "Group_final", 
                           idents = "Epi-Glandular secretory",
                           raster = FALSE,
                           y.max = y_max) + 
  xlab("Epitelio Glandular secretor") + ylab("Nivel de actividad") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Stroma_ESR1 | Vln_Lum_ESR1 | Vln_Cil_ESR1 | Vln_Glan_ESR1 | Vln_GlanSec_ESR1

# Calcular los Violin plots de la expresión de los TFs de interés
DefaultAssay(object = seuratObj) <- "RNA"
y_max = 4.5

Vln_Stroma_PGR <- VlnPlot(seuratObj, features = "PGR",
                          pt.size = 0, 
                          group.by = "Group_final", 
                          idents = "Stroma",
                          raster = FALSE,
                          y.max = y_max) + 
  xlab("Estroma") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Lum_PGR <- VlnPlot(seuratObj, features = "PGR",
                       pt.size = 0, 
                       group.by = "Group_final", 
                       idents = "Epi-Lumenal",
                       raster = FALSE,
                       y.max = y_max) + 
  xlab("Epitelio Lumenal") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Cil_PGR <- VlnPlot(seuratObj, features = "PGR",
                       pt.size = 0, 
                       group.by = "Group_final", 
                       idents = "Epi-Ciliated",
                       raster = FALSE,
                       y.max = y_max) + 
  xlab("Epitelio Ciliado") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))


Vln_Glan_PGR <- VlnPlot(seuratObj, features = "PGR",
                        pt.size = 0, 
                        group.by = "Group_final", 
                        idents = "Epi-Glandular",
                        raster = FALSE,
                        y.max = y_max) + 
  xlab("Epitelio Glandular") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))


Vln_GlanSec_PGR <- VlnPlot(seuratObj, features = "PGR",
                           pt.size = 0, 
                           group.by = "Group_final", 
                           idents = "Epi-Glandular secretory",
                           raster = FALSE,
                           y.max = y_max) + 
  xlab("Epitelio Glandular secretor") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Stroma_PGR | Vln_Lum_PGR | Vln_Cil_PGR | Vln_Glan_PGR | Vln_GlanSec_PGR

y_max = 4

Vln_Stroma_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                           pt.size = 0, 
                           group.by = "Group_final", 
                           idents = "Stroma",
                           raster = FALSE,
                           y.max = y_max) + 
  xlab("Estroma") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Lum_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                        pt.size = 0, 
                        group.by = "Group_final", 
                        idents = "Epi-Lumenal",
                        raster = FALSE,
                        y.max = y_max) + 
  xlab("Epitelio Lumenal") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Cil_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                        pt.size = 0, 
                        group.by = "Group_final", 
                        idents = "Epi-Ciliated",
                        raster = FALSE,
                        y.max = y_max) + 
  xlab("Epitelio Ciliado") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Glan_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                         pt.size = 0, 
                         group.by = "Group_final", 
                         idents = "Epi-Glandular",
                         raster = FALSE,
                         y.max = y_max) + 
  xlab("Epitelio Glandular") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_GlanSec_ESR1 <- VlnPlot(seuratObj, features = "ESR1",
                            pt.size = 0, 
                            group.by = "Group_final", 
                            idents = "Epi-Glandular secretory",
                            raster = FALSE,
                            y.max = y_max) + 
  xlab("Epitelio Glandular secretor") + ylab("Nivel de expresión") + stat_compare_means(comparisons = comparisons, label = "p.signif") + theme(axis.title = element_text(17))

Vln_Stroma_ESR1 | Vln_Lum_ESR1 | Vln_Cil_ESR1 | Vln_Glan_ESR1 | Vln_GlanSec_ESR1


## Hacer heatmap con el paquete ClomplexHeatmap
# Se crean unos identificadores que separan los tipos celulares por condición
table(seuratObj@meta.data$cell_type_zooming_v5, seuratObj@meta.data$Group_final)
seuratObj@meta.data$celltype_aggregate <- paste(seuratObj@meta.data$cell_type_zooming_v5, seuratObj@meta.data$Group_final, sep = "_") # Adaptar a los datos del usuario
seuratObj@meta.data$celltype_aggregate %>% table() %>% sort(decreasing = TRUE)
celltype_id <- "celltype_aggregate" # nombre de la columna de metadatos del tipo celular de interés
seuratObj <- SetIdent(seuratObj, value = seuratObj[[celltype_id]])

# Ver cual es la actividad media por grupo de los 50 TFs más variables:
n_tfs <- 50
# Extraer las actividades del objeto como un dataframe long
df <- t(as.matrix(seuratObj@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(seuratObj)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Obtener los TFs con las medias más variables entre clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Hacer un subset del dataframe long para los top TFs y transformarlo a una matriz wide y hacer la transformada de la matriz
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

top_acts_mat_t <- t(top_acts_mat)


# Elegir la paleta de colores
palette_length <- 100
my_color <- colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Crear el heatmap
Heatmap(top_acts_mat_t, name = "mat", rect_gp = gpar(col = "white", lwd = 2), column_names_side = "top", column_dend_side = "bottom", column_km = 5, row_km = 4)


## Heatmap solo con los tipos celulares de interés
n_tfs <- 50
# Extraer las actividades del objeto como un dataframe long
df <- t(as.matrix(seuratObj@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(seuratObj)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Obtener los TFs con las medias más variables entre clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Hacer un subset del dataframe long para los top TFs y transformarlo a una matriz wide y hacer la trnasformada de la matriz
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  filter(cluster == "Stroma_AS" | cluster == "Stroma_WOI" | cluster == "Epi-Lumenal_AS" | cluster == "Epi-Lumenal_WOI" | cluster == "Epi-Ciliated_AS" | cluster == "Epi-Ciliated_WOI" | cluster == "Epi-Glandular_AS" | cluster == "Epi-Glandular_WOI" | cluster == "Epi-Glandular secretory_AS" | cluster == "Epi-Glandular secretory_WOI") %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

top_acts_mat_t <- t(top_acts_mat)

# Elegir paleta de colores
palette_length <- 100
my_color <- colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Generar el heatmap
Heatmap(top_acts_mat_t, name = "mat", rect_gp = gpar(col = "white", lwd = 2), column_names_side = "top", column_dend_side = "bottom", column_km = 4)

## Pseudobulk, conteos según tipo celular y muestra (unión de las muestras con la condición)
seuratObj$orig.ident <- gsub("_", ".", seuratObj$orig.ident)
seuratObj$samples <- paste0(seuratObj$Group_final, seuratObj$orig.ident)
DefaultAssay(seuratObj) <- "RNA"

counts <- AggregateExpression(seuratObj, group.by = c("cell_type_zooming_v5", "samples"), assays = "RNA", slot = "counts", return.seurat = FALSE)
View(counts$RNA)
# transponer
counts.t <- t(counts$RNA)
# convertir a dataframe
counts.t <- as.data.frame(counts.t)

# Coger solo el nombre del tipo celular de los nombres de fila
splitRows <- gsub("_.*", "", rownames(counts.t))

# Lista de las matrices de conteos de cada tipo celular
counts.split <- split.data.frame(counts.t, f = factor(splitRows))

# Fijar los nombres de columna y transponer
counts.split.modified <- lapply(counts.split, function(x){
  rownames(x) <- gsub(".*_(.*)", "\\1", rownames(x))
  t(x)
  })

# Correr el análisis DE con las células del Estroma
# 1.Obtener la matriz de conteos
counts_stroma <- counts.split.modified$Stroma

# 2.Generar niveles de metadatos de las muestras (asignar a cada muestra su condición)
colData <- data.frame(samples = colnames(counts_stroma))
colData <- colData %>% mutate(condition = ifelse(grepl("WOI", samples), "WOI", "AS")) %>% column_to_rownames(var = "samples")

# Realizar DeSeq2
# obtener objeto Deseq2 
dds <- DESeqDataSetFromMatrix(countData = counts_stroma, colData = colData, design = ~ condition)

# Filtrar
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# correr DESeq2
dds <- DESeq(dds)

# Comprobar los coeficientes de la comparación
resultsNames(dds)

# Generar un objeto con los resultados
res <- results(dds, name = "condition_WOI_vs_AS")
res

# Extraer las columnas de interés por gen
res_data <- as.data.frame(res)
deg <- res_data %>% 
  select(log2FoldChange, stat, pvalue) %>%
  filter(!is.na(stat)) %>%
  as.matrix()
head(deg)
  
# Red de regulones de CollecTRI
net <- get_collectri(organism='human', split_complexes=FALSE)
net

# Se puede visualizar los TG más diferenciales en cada TF a través de sus p-values para interpretar los resultados 
tf <- 'PGR'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_stroma_PGR <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

tf <- 'ESR1'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_stroma_ESR1 <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

# Unir los df filtrados y crear una columna que indique la condición en la que se expresan los TG según el log2FoldChange
df_filtered_stroma <- rbind(df_filtered_stroma_PGR, df_filtered_stroma_ESR1)
df_filtered_stroma <- df_filtered_stroma %>% mutate(CellType = "Stroma")
df_filtered_stroma_AS <- df_filtered_stroma %>% filter(log2FoldChange < 0) %>% mutate(Group = "AS") 
df_filtered_stroma_WOI <- df_filtered_stroma %>% filter(log2FoldChange > 0) %>% mutate(Group = "WOI") 
df_filtered_stroma <- rbind(df_filtered_stroma_AS, df_filtered_stroma_WOI)
# Correr el análisis DE con las células del Epitelio Luminal
# 1.Obtener la matriz de conteos
counts_lumenal <- counts.split.modified$`Epi-Lumenal`

# 2.Generar niveles de metadatos de las muestras (asignar a cada muestra su condición)
colData <- data.frame(samples = colnames(counts_lumenal))
colData <- colData %>% mutate(condition = ifelse(grepl("WOI", samples), "WOI", "AS")) %>% column_to_rownames(var = "samples")

# Realizar DeSeq2
# obtener objeto Deseq2 
dds <- DESeqDataSetFromMatrix(countData = counts_lumenal, colData = colData, design = ~ condition)

# filtrar
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# correr DESeq2
dds <- DESeq(dds)

# Comprobar los coeficientes de la comparación
resultsNames(dds)

# Generar un objeto con los resultados
res <- results(dds, name = "condition_WOI_vs_AS")
res

# Extraer las columnas de interés por gen
res_data <- as.data.frame(res)
deg <- res_data %>% 
  select(log2FoldChange, stat, pvalue) %>%
  filter(!is.na(stat)) %>%
  as.matrix()
head(deg)

# Se puede visualizar los TG más diferenciales en cada TF a través de sus p-values para interpretar los resultados 
tf <- 'PGR'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_lumenal_PGR <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

tf <- 'ESR1'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_lumenal_ESR1 <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

# Unir los df filtrados y crear una columna que indique la condición en la que se expresan los TG según el log2FoldChange
df_filtered_lumenal <- rbind(df_filtered_lumenal_PGR, df_filtered_lumenal_ESR1)
df_filtered_lumenal <- df_filtered_lumenal %>% mutate(CellType = "Epi-Lumenal")
df_filtered_lumenal_AS <- df_filtered_lumenal %>% filter(log2FoldChange < 0) %>% mutate(Group = "AS")
df_filtered_lumenal_WOI <- df_filtered_lumenal %>% filter(log2FoldChange > 0) %>% mutate(Group = "WOI")
df_filtered_lumenal <- rbind(df_filtered_lumenal_AS, df_filtered_lumenal_WOI)
# Correr el análisis DE con las células del Epitelio Ciliado
# 1.Obtener la matriz de conteos
counts_ciliated <- counts.split.modified$`Epi-Ciliated`

# 2.Generar niveles de metadatos de las muestras (asignar a cada muestra su condición)
colData <- data.frame(samples = colnames(counts_ciliated))
colData <- colData %>% mutate(condition = ifelse(grepl("WOI", samples), "WOI", "AS")) %>% column_to_rownames(var = "samples")

# Realizar DeSeq2
# obtener objeto Deseq2 
dds <- DESeqDataSetFromMatrix(countData = counts_ciliated, colData = colData, design = ~ condition)

# filtrar
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# correr DESeq2
dds <- DESeq(dds)

# Comprobar los coeficientes de la comparación
resultsNames(dds)

# Generar un objeto con los resultados
res <- results(dds, name = "condition_WOI_vs_AS")
res

# Extraer las columnas de interés por gen
res_data <- as.data.frame(res)
deg <- res_data %>% 
  select(log2FoldChange, stat, pvalue) %>%
  filter(!is.na(stat)) %>%
  as.matrix()
head(deg)

# Se puede visualizar los TG más diferenciales en cada TF a través de sus p-values para interpretar los resultados  
tf <- 'PGR'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_ciliated_PGR <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

tf <- 'ESR1'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_ciliated_ESR1 <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

# Unir los df filtrados y crear una columna que indique la condición en la que se expresan los TG según el log2FoldChange
df_filtered_ciliated <- rbind(df_filtered_ciliated_PGR, df_filtered_ciliated_ESR1)
df_filtered_ciliated <- df_filtered_ciliated %>% mutate(CellType = "Epi-Ciliated")
df_filtered_ciliated_AS <- df_filtered_ciliated %>% filter(log2FoldChange < 0) %>% mutate(Group = "AS")
df_filtered_ciliated_WOI <- df_filtered_ciliated %>% filter(log2FoldChange > 0) %>% mutate(Group = "WOI")
df_filtered_ciliated <- rbind(df_filtered_ciliated_AS, df_filtered_ciliated_WOI)
# Correr el análisis DE con las células del Epitelio Glandular
# 1.Obtener la matriz de conteos
counts_glandular <- counts.split.modified$`Epi-Glandular`

# 2.Generar niveles de metadatos de las muestras (asignar a cada muestra su condición)
colData <- data.frame(samples = colnames(counts_glandular))
colData <- colData %>% mutate(condition = ifelse(grepl("WOI", samples), "WOI", "AS")) %>% column_to_rownames(var = "samples")

# Realizar DeSeq2
# obtener objeto Deseq2 
dds <- DESeqDataSetFromMatrix(countData = counts_glandular, colData = colData, design = ~ condition)

# filtrar
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# correr DESeq2
dds <- DESeq(dds)

# Comprobar los coeficientes de la comparación
resultsNames(dds)

# Generar un objeto con los resultados
res <- results(dds, name = "condition_WOI_vs_AS")
res

# Extraer las columnas de interés por gen
res_data <- as.data.frame(res)
deg <- res_data %>% 
  select(log2FoldChange, stat, pvalue) %>%
  filter(!is.na(stat)) %>%
  as.matrix()
head(deg)

# Se puede visualizar los TG más diferenciales en cada TF a través de sus p-values para interpretar los resultados  
tf <- 'PGR'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_glandular_PGR <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

tf <- 'ESR1'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_glandular_ESR1 <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

# Unir los df filtrados y crear una columna que indique la condición en la que se expresan los TG según el log2FoldChange
df_filtered_glandular <- rbind(df_filtered_glandular_PGR, df_filtered_glandular_ESR1)
df_filtered_glandular <- df_filtered_glandular %>% mutate(CellType = "Epi-Glandular")
df_filtered_glandular_AS <- df_filtered_glandular %>% filter(log2FoldChange < 0) %>% mutate(Group = "AS")
df_filtered_glandular_WOI <- df_filtered_glandular %>% filter(log2FoldChange > 0) %>% mutate(Group = "WOI")
df_filtered_glandular <- rbind(df_filtered_glandular_AS, df_filtered_glandular_WOI)

# Correr el análisis DE con las células del Epitelio Glandular secretor
# 1.Obtener la matriz de conteos
counts_glandular.secretory <- counts.split.modified$`Epi-Glandular secretory`

# 2.Generar niveles de metadatos de las muestras (asignar a cada muestra su condición)
colData <- data.frame(samples = colnames(counts_glandular.secretory))
colData <- colData %>% mutate(condition = ifelse(grepl("WOI", samples), "WOI", "AS")) %>% column_to_rownames(var = "samples")

# Realizar DeSeq2
# obtener objeto Deseq2 
dds <- DESeqDataSetFromMatrix(countData = counts_glandular.secretory, colData = colData, design = ~ condition)

# filtrar
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# correr DESeq2
dds <- DESeq(dds)

# Comprobar los coeficientes de la comparación
resultsNames(dds)

# Generar un objeto con los resultados
res <- results(dds, name = "condition_WOI_vs_AS")
res

# Extraer las columnas de interés por gen
res_data <- as.data.frame(res)
deg <- res_data %>% 
  select(log2FoldChange, stat, pvalue) %>%
  filter(!is.na(stat)) %>%
  as.matrix()
head(deg)

# Se puede visualizar los TG más diferenciales en cada TF a través de sus p-values para interpretar los resultados 
tf <- 'PGR'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_glandular.secretory_PGR <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

tf <- 'ESR1'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('log2FoldChange', 'stat', 'pvalue')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

# Obtener un df filtrado por pvalue
df_filtered_glandular.secretory_ESR1 <- df[order(df$pvalue),] %>% filter(pvalue <= 1e-05)

ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  labs(title = tf)

# Unir los df filtrados y crear una columna que indique la condición en la que se expresan los TG según el log2FoldChange
df_filtered_glandular.secretory <- rbind(df_filtered_glandular.secretory_PGR, df_filtered_glandular.secretory_ESR1)
df_filtered_glandular.secretory <- df_filtered_glandular.secretory %>% mutate(CellType = "Epi-Glandular.secretory")
df_filtered_glandular.secretory_AS <- df_filtered_glandular.secretory %>% filter(log2FoldChange < 0) %>% mutate(Group = "AS")
df_filtered_glandular.secretory_WOI <- df_filtered_glandular.secretory %>% filter(log2FoldChange > 0) %>% mutate(Group = "WOI")
df_filtered_glandular.secretory <- rbind(df_filtered_glandular.secretory_AS, df_filtered_glandular.secretory_WOI)

df_filtered_total <- rbind(df_filtered_stroma, df_filtered_lumenal, df_filtered_ciliated, df_filtered_glandular, df_filtered_glandular.secretory)
# Reordenar las columnas usando dplyr
df_filtered_total <- df_filtered_total %>% select(source, ID, mor, color, log2FoldChange, stat, pvalue, CellType, Group)
rownames(df_filtered_total) <- NULL
write.csv(df_filtered_total, "../Results/df_TF/df_filtered_total.csv", row.names = FALSE)

# Filtrar los genes
names <- df_filtered_total %>% select(source, ID) %>% rename(name = ID) %>% distinct()

relations <- df_filtered_total %>% select(source, ID, color, log2FoldChange, CellType, Group) %>% rename(from = source, to = ID) %>% distinct()
relations <- relations %>% mutate(color = ifelse(color == 1, "up", "down")) 
relations <- relations %>% rename(regulation = color)

# Crear el grafo con igraph
set_graph_style(plot_margin = margin(1,1,1,1))
graph <- as_tbl_graph(relations, directed = TRUE, vertices = names)

ggraph(graph, layout = "kk") + geom_edge_link(aes(start_cap = label_rect(node1.name), end_cap = label_rect(node2.name), colour = regulation, width = abs(log2FoldChange)), 
                                              arrow = arrow(type = "closed", length = unit(1, 'mm'))) +  
  geom_node_text(aes(label = name)) + facet_edges(CellType~Group, drop = T) +
  theme_graph(base_family = "Arial", base_size = 15, background = "white", foreground = 'darkorange2', border = TRUE, text_colour = "black", bg_text_colour = "black", fg_text_colour = "black",
  strip_text_family = "Arial", strip_text_size = 20, strip_text_face = "bold", strip_text_colour = "black", plot_margin = margin(30, 30, 30, 30)) +
  scale_edge_color_manual(values = c('down' = 'steelblue', 'up' = 'indianred'))  
