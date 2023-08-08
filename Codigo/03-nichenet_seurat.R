
# Paquetes requeridos
library(nichenetr)
library(Seurat) 
library(tidyverse)
library(circlize)
library(SeuratDisk)
library(OmnipathR)

# Cargar la matriz ligando-target, la red de interacción ligando-receptor y los pesos de las redes de interacción de NicheNet:
ligand_target_matrix = readRDS("../Results/nichenet_results/ligand_target_matrix_weighted.rds")

networks = readRDS("../Results/nichenet_results/networks.rds")
lr_network = networks$lr_network

weighted_networks = readRDS("../Results/nichenet_results/weighted_networks_weighted.rds")

# Cargar los datos de expresión de las células que iteraccionan y reducirlos para que no ocupen tanto espacio:
seuratObj = LoadH5Seurat("../Data/Single-cell/07_sc_AS_Estonia_Wang_Group_SCT_Main_lineages_clean.H5Seurat")
seuratObj <- DietSeurat(object = seuratObj, assays = "RNA", dimreducs = "umap")
seuratObj$Group_final <- gsub(pattern = "WOI WOI Control", replacement = "WOI", x = seuratObj$Group_final)

## Realizar el análisis con NicheNet
# Explicar la expresión diferencial entre dos condiciones
# Definir los tipos celulares emisores y los receptores
sender_celltypes <- c("Epi-Lumenal", "Epi-Glandular secretory", "Epi-Glandular", "Epi-Ciliated")
receiver_celltypes <- c("dStroma")
# Combinar los datos de expresión de Seurat con los datos de NicheNet
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = receiver_celltypes, 
  condition_colname = "Group_final", condition_oi = "AS", condition_reference = "WOI", 
  sender = sender_celltypes, 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)

# Ver que ligandos priorizados son predichos para regular la expresión de genes diferencialmente expresados:
heatmap <- nichenet_output$ligand_target_heatmap
heatmap2 <- heatmap
# Se filtran los genes de interés
heatmap2$data <- heatmap$data %>% filter( x == "PGR" | x == "ESR1")
heatmap2 <- heatmap2 + xlab("Genes diana predichos") + ylab("Ligandos priorizados") + guides(fill = guide_legend(title="Potencial regulatorio")) + theme(axis.title = element_text(size = 14), legend.title = element_text(size = 12))

# Calcular el logFC de la expresión de estos ligandos priorizados en cada tipo celular de los definidos anteriormente como emisores
ligand_differential_expression_heatmap <- nichenet_output$ligand_differential_expression_heatmap
ligand_differential_expression_heatmap <- ligand_differential_expression_heatmap + xlab("LFC en emisores") + ylab("Ligandos priorizados")

# Ver el nivel de expresión de los ligandos para cada tipo celular:
ligand_expression_dotplot <- nichenet_output$ligand_expression_dotplot 
ligand_expression_dotplot$guides$size$title <- "Porcentaje expresado"
ligand_expression_dotplot$guides$colour$title <- "Expresión media"
ligand_expression_dotplot <- ligand_expression_dotplot + xlab("Ligandos priorizados") + ylab("Tipo celular") + theme(axis.title = element_text(size = 14), legend.title = element_text(size = 12))

# Visualización de las relaciones ligando-receptor que se han descrito en la literatura y no sebasan en interacciones predichas (bona fide):
ligand_receptor_heatmap <- nichenet_output$ligand_receptor_heatmap
ligand_receptor_heatmap$labels$fill <- "Potencial de interacción previa"
ligand_receptor_heatmap <- ligand_receptor_heatmap + xlab("Receptores") + ylab("Ligandos")

# Heatmap de la actividad de los ligandos prioritarios que se han visto anteriormente para los genes de interés
ligand_aupr_matrix = nichenet_output$ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(nichenet_output$ligand_activities$test_ligand)

rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

ligands <- rownames(ligand_aupr_matrix)
order_ligands <- rev(ligands)

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Ligandos priorizados","Actividad del ligando", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(capacidad de predicción del gen diana)") + theme(legend.text = element_text(size = 9))
p_ligand_aupr <- p_ligand_aupr + guides(fill = guide_legend(title="AUPR\n(capacidad de predicción del gen diana)")) + theme(axis.title = element_text(size = 14), legend.title = element_text(size = 12))

## Figura combinada del gráfico con la actividad de los ligandos, el nivel de expresión de los ligandos en los tipos celulares sender, su logFC y los ligandos para los genes de interés
figures_without_legend = cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  ligand_expression_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expresión en emisores") + xlab("") + scale_y_discrete(position = "right"),
  ligand_differential_expression_heatmap + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  heatmap2 + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(3, 7, 7, 12))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(ligand_expression_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(ligand_differential_expression_heatmap)),
  ggpubr::as_ggplot(ggpubr::get_legend(heatmap2)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

## Se repiten los mismos pasos invirtiendo las células emisoras y receptoras
# Definir los tipos celulares que envian la señal y los receptores
sender_celltypes <- c("dStroma")
receiver_celltypes <- c("Epi-Lumenal", "Epi-Glandular secretory", "Epi-Glandular", "Epi-Ciliated")

# Combinar los datos de expresión de Seurat con los datos de NicheNet
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = receiver_celltypes, 
  condition_colname = "Group_final", condition_oi = "AS", condition_reference = "WOI", 
  sender = sender_celltypes, 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)

# Ver que ligandos top son predichos para regular la expresión de genes diferencialmente expresados:
heatmap <- nichenet_output$ligand_target_heatmap
heatmap2 <- heatmap
# Se filtran los genes de interés
heatmap2$data <- heatmap$data %>% filter( x == "PGR" | x == "ESR1")
heatmap2 <- heatmap2 + xlab("Genes diana predichos") + ylab("Ligandos priorizados") + guides(fill = guide_legend(title="Potencial regulatorio")) + theme(axis.title = element_text(size = 14), legend.title = element_text(size = 12))

# Calcular el logFC de la expresión de estos ligandos prioritarios en cada tipo celular de los definidos anteriormente como sender
ligand_differential_expression_heatmap <- nichenet_output$ligand_differential_expression_heatmap
ligand_differential_expression_heatmap <- ligand_differential_expression_heatmap + xlab("LFC en emisores") + ylab("Ligandos priorizados")

# Ver el nivel de expresión de los ligandos para cada tipo celular:
ligand_expression_dotplot <- nichenet_output$ligand_expression_dotplot 
ligand_expression_dotplot$guides$size$title <- "Porcentaje expresado"
ligand_expression_dotplot$guides$colour$title <- "Expresión media"
ligand_expression_dotplot <- ligand_expression_dotplot + xlab("Ligandos priorizados") + ylab("Tipo celular") + theme(axis.title = element_text(size = 14), legend.title = element_text(size = 12))

# Visualización de las relaciones ligando-receptor que se han descrito en la literatura y no sebasan en interacciones predichas (bona fide):
ligand_receptor_heatmap <- nichenet_output$ligand_receptor_heatmap
ligand_receptor_heatmap$labels$fill <- "Potencial de interacción previa"
ligand_receptor_heatmap <- ligand_receptor_heatmap + xlab("Receptores") + ylab("Ligandos")

# Heatmap de la actividad de los ligandos prioritarios que se han visto anteriormente para los genes de interés
ligand_aupr_matrix = nichenet_output$ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(nichenet_output$ligand_activities$test_ligand)

rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

ligands <- rownames(ligand_aupr_matrix)
order_ligands <- rev(ligands)

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Ligandos priorizados","Actividad del ligando", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_aupr <- p_ligand_aupr + guides(fill = guide_legend(title="AUPR\n(capacidad de predicción del gen diana)")) + theme(axis.title = element_text(size = 14), legend.title = element_text(size = 12))

## Figura combinada del gráfico con la actividad de los ligandos, el nivel de expresión de los ligandos en los tipos celulares sender, su logFC y los ligandos para los genes de interés
figures_without_legend = cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  ligand_expression_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expresión en emisores") + xlab("") + scale_y_discrete(position = "right"),
  ligand_differential_expression_heatmap + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  heatmap2 + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(3, 7, 7, 12))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(ligand_expression_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(ligand_differential_expression_heatmap)),
  ggpubr::as_ggplot(ggpubr::get_legend(heatmap2)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot