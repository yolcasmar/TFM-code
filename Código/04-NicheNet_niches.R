
# Paquetes necesarios
library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) 
library(SeuratDisk)

## Análisis diferencial con NicheNet entre condiciones de interés
## Cargar los datos de expresión
seuratObj <- LoadH5Seurat("../Data/Single-cell/07_sc_AS_Estonia_Wang_Group_SCT_Main_lineages_clean.H5Seurat")
seuratObj <- DietSeurat(object = seuratObj, assays = c("RNA", "SCT"), dimreducs = "umap")
seuratObj$Group_final <- gsub(pattern = "WOI WOI Control", replacement = "WOI", x = seuratObj$Group_final)

# Dividir los tipos celulares según las condiciones a analizar y crear una nueva columna en los metadatos
table(seuratObj@meta.data$cell_type_zooming_v5, seuratObj@meta.data$Group_final) # tipo celular vs condición
seuratObj@meta.data$celltype_aggregate <- paste(seuratObj@meta.data$cell_type_zooming_v5, seuratObj@meta.data$Group_final, sep = "_")
celltype_id <- "celltype_aggregate" # Nombre de la columna de metadatos del tipo celular de interés
seuratObj <- SetIdent(seuratObj, value = seuratObj[[celltype_id]])

# Cargar la red de interacción ligando-receptor y la matriz ligando-target de NicheNet
ligand_target_matrix <- readRDS("../Results/nichenet_results/ligand_target_matrix_weighted.rds")

network <- readRDS("../Results/nichenet_results/networks.rds")
lr_network <- network$lr_network
lr_network <- lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network <- lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)

## 1. Definir los grupos de interés
niches <- list(
  "AS_niche" = list(
    "sender" = c("Epi-Lumenal_AS", "Epi-Ciliated_AS", "Epi-Glandular_AS", "Epi-Glandular secretory_AS"),
    "receiver" = c("Stroma_AS")),
  "WOI_niche" = list(
    "sender" = c("Epi-Lumenal_WOI",  "Epi-Ciliated_WOI", "Epi-Glandular_WOI", "Epi-Glandular secretory_WOI"),
    "receiver" = c("Stroma_WOI"))
  ) 


## 2. Calcular la expresión diferencial entre los grupos
assay_oi <- "SCT" 
seuratObj <- PrepSCTFindMarkers(seuratObj, assay = "SCT", verbose = FALSE) # Necesario antes de la función calculate_niche_de

DE_sender <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # Solo los ligandos importantes para los tipos celulares sender

DE_receiver <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # Solo los receptores ahora

DE_sender <- DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver <- DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))

# Procesar los resultados de la expresión diferencial:
expression_pct <- 0.10
DE_sender_processed <- process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed <- process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")

# Combinar la expresión diferencial de sender-receiver basado en los pares ligando-receptor:
specificity_score_LR_pairs <- "min_lfc"
DE_sender_receiver <- combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)

## 3. Calcular la expresión diferencial entre las diferentes regiones espaciales
include_spatial_info_sender <- FALSE 
include_spatial_info_receiver <- FALSE 

# Como no hay información espacial:
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
  spatial_info <- tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 

if(include_spatial_info_sender == TRUE){
  sender_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed <- process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  sender_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
} else {
  sender_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  }

if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed <- process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  receiver_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  
} else {
  receiver_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }

## 4. Calcular las actividades de los ligandos e inferir las relaciones ligando-target activas
lfc_cutoff <- 0.15 # recomendado para 10x como min_lfc cutoff 
specificity_score_targets <- "min_lfc"

DE_receiver_targets <- calculate_niche_de_targets(seurat_obj = seuratObj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 

DE_receiver_processed_targets <- process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background <- DE_receiver_processed_targets %>% pull(target) %>% unique()
geneset_niche1 <- DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 <- DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()

geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
length(geneset_niche1)
length(geneset_niche2)

top_n_target <- 250

niche_geneset_list <- list(
  "AS_niche" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "WOI_niche" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
  )

ligand_activities_targets <- get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)

## 5. Calcular la expresión escalada de los ligandos, receptores y targets en los tipos celulares de interés
features_oi <- union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)

dotplot <- suppressWarnings(Seurat::DotPlot(seuratObj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl <- dotplot$data %>% as_tibble()
exprs_tbl <- exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

exprs_tbl_ligand <- exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor <- exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target <- exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

exprs_tbl_ligand <- exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor <- exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

## 6. Fracción de expresión y receptor
exprs_sender_receiver <- lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 


## 7. Priorización de interacciones ligando-receptor y ligando-target
prioritizing_weights <- c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)

output <- list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)

prioritization_tables <- get_prioritization_tables(output, prioritizing_weights)

## 8. Visualización del output de NicheNet

# Expresión diferencial de los ligandos
top_ligand_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche

receiver_oi <- "Stroma_AS" 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

# log_FC mínimo comparado con el otro grupo
lfc_plot <- make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = TRUE, heights = NULL, widths = NULL)
lfc_plot$labels$x <- "LFC del receptor\n en Receptores"
lfc_plot$labels$fill <- "Receptor:\nLFC min vs\notros niches"
lfc_plot[[1]]$labels$y <- "Pares Ligando-Receptor priorizados"
lfc_plot[[1]]$labels$x <- "LFC del ligando\n en Emisores"
lfc_plot[[1]]$labels$fill <- "Ligando:\nLFC min vs\notros niches"
lfc_plot

# Lo mismo para el otro grupo
receiver_oi <- "Stroma_WOI"  
filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(50, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

lfc_plot2 <- make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot2$labels$x <- "LFC del receptor\n en Receptores"
lfc_plot2$labels$fill <- "Receptor:\nLFC min vs\notros niches"
lfc_plot2[[1]]$labels$y <- "Pares Ligando-Receptor priorizados"
lfc_plot2[[1]]$labels$x <- "LFC del ligando\n en Emisores"
lfc_plot2[[1]]$labels$fill <- "Ligando:\nLFC min vs\notros niches"
lfc_plot2

## Pasos para los grupos anteriores invertidos, los tipos celulares sender como receptores y viceversa
# Se siguen los pasos anteriores hasta obtener las tablas con las interacciones prioritarias para cada interacción estroma epitelio
# Los epitelios se calculan uno a uno porque no puede usarse más de un tipo celular receptor para calcular la expresión diferencial

## Stroma --> Epi_Lumenal
## 1. Definir los grupos de interés
niches <- list(
  "AS_niche" = list(
    "sender" = c("Stroma_AS"),
    "receiver" = c("Epi-Lumenal_AS")),
  "WOI_niche" = list(
    "sender" = c("Stroma_WOI"),
    "receiver" = c("Epi-Lumenal_WOI"))
   )

## 2. Calcular la expresión diferencial entre los grupos
assay_oi <- "SCT" 

DE_sender <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # Solo los ligandos importantes para los tipos celulares sender

DE_receiver <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # Solo los receptores ahora

DE_sender <- DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver <- DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))

# Procesar los resultados de la expresión diferencial:
expression_pct <- 0.10
DE_sender_processed <- process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed <- process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")

# Combinar la expresión diferencial de sender-receiver basado en los pares ligando-receptor:
specificity_score_LR_pairs <- "min_lfc"
DE_sender_receiver <- combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)

## 3. Calcular la expresión diferencial entre las diferentes regiones espaciales
include_spatial_info_sender <- FALSE 
include_spatial_info_receiver <- FALSE 

# Como no hay información espacial:
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
  spatial_info <- tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 

if(include_spatial_info_sender == TRUE){
  sender_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed <- process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  sender_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
} else {
  sender_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  }

if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed <- process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  receiver_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
} else {
  receiver_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }

## 4. Calcular las actividades de los ligandos e inferir las relaciones ligando-target activas
lfc_cutoff <- 0.15 # recomendado para 10x como min_lfc cutoff 
specificity_score_targets <- "min_lfc"

DE_receiver_targets <- calculate_niche_de_targets(seurat_obj = seuratObj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 

DE_receiver_processed_targets <- process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background <- DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 <- DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 <- DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()

geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
length(geneset_niche1)
length(geneset_niche2)

top_n_target <- 250

niche_geneset_list <- list(
  "AS_niche" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "WOI_niche" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
  )

ligand_activities_targets <- get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)

## 5. Calcular la expresión escalada de los ligandos, receptores y targets en los tipos celulares de interés
features_oi <- union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)

dotplot <- suppressWarnings(Seurat::DotPlot(seuratObj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl <- dotplot$data %>% as_tibble()
exprs_tbl <- exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

exprs_tbl_ligand <- exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor <- exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target <- exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

exprs_tbl_ligand <- exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor <- exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

## 6. Fracción de expresión y receptor
exprs_sender_receiver <- lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 

## 7. Priorización de interacciones ligando-receptor y ligando-target
prioritizing_weights <- c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)

output <- list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)

prioritization_tables <- get_prioritization_tables(output, prioritizing_weights)
saveRDS(prioritization_tables, file = "../Results/priorization_tables/stroma_epilumenal.rds")

## 8. Visualización del output de NicheNet
# Expresión diferencial de los ligandos
top_ligand_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche

# Se guarda la tabla con las interacciones prioritarias para ambos grupos
receiver_oi <- c("Epi-Lumenal_AS") 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(10, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

saveRDS(prioritized_tbl_oi, file = "../Results/df_stroma_epi/stroma_epilumenal_AS.rds")

receiver_oi <- c("Epi-Lumenal_WOI") 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(10, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

saveRDS(prioritized_tbl_oi, file = "../Results/df_stroma_epi/stroma_epilumenal_WOI.rds")

## Stroma --> Epi_Ciliated
## 1. Definir los grupos de interés
niches <- list(
  "AS_niche" = list(
    "sender" = c("Stroma_AS"),
    "receiver" = c("Epi-Ciliated_AS")),
  "WOI_niche" = list(
    "sender" = c("Stroma_WOI"),
    "receiver" = c("Epi-Ciliated_WOI"))
   )

## 2. Calcular la expresión diferencial entre los grupos
assay_oi <- "SCT" 

DE_sender <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # Solo los ligandos importantes para los tipos celulares sender

DE_receiver <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # Solo los receptores ahora

DE_sender <- DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver <- DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))

# Procesar los resultados de la expresión diferencial:
expression_pct <- 0.10
DE_sender_processed <- process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed <- process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")

# Combinar la expresión diferencial de sender-receiver basado en los pares ligando-receptor:
specificity_score_LR_pairs <- "min_lfc"
DE_sender_receiver <- combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)

## 3. Calcular la expresión diferencial entre las diferentes regiones espaciales
include_spatial_info_sender <- FALSE 
include_spatial_info_receiver <- FALSE 

# Como no hay información espacial:
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
  spatial_info <- tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 

if(include_spatial_info_sender == TRUE){
  sender_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed <- process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  sender_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
} else {
  sender_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  }

if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed <- process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  receiver_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
} else {
  receiver_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }

## 4. Calcular las actividades de los ligandos e inferir las relaciones ligando-target activas
lfc_cutoff <- 0.15 # recomendado para 10x como min_lfc cutoff 
specificity_score_targets <- "min_lfc"

DE_receiver_targets <- calculate_niche_de_targets(seurat_obj = seuratObj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 

DE_receiver_processed_targets <- process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background <- DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 <- DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 <- DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()

geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
length(geneset_niche1)
length(geneset_niche2)

top_n_target <- 250

niche_geneset_list <- list(
  "AS_niche" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "WOI_niche" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
  )

ligand_activities_targets <- get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)

## 5. Calcular la expresión escalada de los ligandos, receptores y targets en los tipos celulares de interés
features_oi <- union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)

dotplot <- suppressWarnings(Seurat::DotPlot(seuratObj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl <- dotplot$data %>% as_tibble()
exprs_tbl <- exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

exprs_tbl_ligand <- exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor <- exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target <- exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

exprs_tbl_ligand <- exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor <- exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

## 6. Fracción de expresión y receptor
exprs_sender_receiver <- lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 


## 7. Priorización de interacciones ligando-receptor y ligando-target
prioritizing_weights <- c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)

output <- list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)

prioritization_tables <- get_prioritization_tables(output, prioritizing_weights)
saveRDS(prioritization_tables, file = "../Results/priorization_tables/stroma_epiciliated.rds")

## 8. Visualización del output de NicheNet
# Expresión diferencial de los ligandos
top_ligand_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche

# Se guarda la tabla con las interacciones prioritarias para ambos grupos
receiver_oi <- c("Epi-Ciliated_AS") 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(10, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

saveRDS(prioritized_tbl_oi, file = "../Results/df_stroma_epi/stroma_epiciliated_AS.rds")

receiver_oi <- c("Epi-Ciliated_WOI") 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(10, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

saveRDS(prioritized_tbl_oi, file = "../Results/df_stroma_epi/stroma_epiciliated_WOI.rds")

## Stroma --> Epi_Glandular
## 1. Definir los grupos de interés
niches <- list(
  "AS_niche" = list(
    "sender" = c("Stroma_AS"),
    "receiver" = c("Epi-Glandular_AS")),
  "WOI_niche" = list(
    "sender" = c("Stroma_WOI"),
    "receiver" = c("Epi-Glandular_WOI"))
   ) 

## 2. Calcular la expresión diferencial entre los grupos
assay_oi <- "SCT" 

DE_sender <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # Solo los ligandos importantes para los tipos celulares sender

DE_receiver <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # Solo los receptores ahora

DE_sender <- DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver <- DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))

# Procesar los resultados de la expresión diferencial:
expression_pct <- 0.10
DE_sender_processed <- process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed <- process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")

# Combinar la expresión diferencial de sender-receiver basado en los pares ligando-receptor:
specificity_score_LR_pairs <- "min_lfc"
DE_sender_receiver <- combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)

## 3. Calcular la expresión diferencial entre las diferentes regiones espaciales
include_spatial_info_sender <- FALSE 
include_spatial_info_receiver <- FALSE 

# Como no hay información espacial:
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
  spatial_info <- tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 

if(include_spatial_info_sender == TRUE){
  sender_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed <- process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  sender_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
} else {
  sender_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  }

if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed <- process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  receiver_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
} else {
  receiver_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }

## 4. Calcular las actividades de los ligandos e inferir las relaciones ligando-target activas
lfc_cutoff <- 0.15 # recomendado para 10x como min_lfc cutoff 
specificity_score_targets <- "min_lfc"

DE_receiver_targets <- calculate_niche_de_targets(seurat_obj = seuratObj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 

DE_receiver_processed_targets <- process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background <- DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 <- DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 <- DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()

geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
length(geneset_niche1)
length(geneset_niche2)

top_n_target = 250

niche_geneset_list <- list(
  "AS_niche" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "WOI_niche" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
  )

ligand_activities_targets <- get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)

## 5. Calcular la expresión escalada de los ligandos, receptores y targets en los tipos celulares de interés
features_oi <- union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)

dotplot <- suppressWarnings(Seurat::DotPlot(seuratObj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl <- dotplot$data %>% as_tibble()
exprs_tbl <- exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

exprs_tbl_ligand <- exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor <- exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target <- exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

exprs_tbl_ligand <- exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor <- exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

## 6. Fracción de expresión y receptor
exprs_sender_receiver <- lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 

## 7. Priorización de interacciones ligando-receptor y ligando-target
prioritizing_weights <- c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)

output <- list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)

prioritization_tables <- get_prioritization_tables(output, prioritizing_weights)
saveRDS(prioritization_tables, file = "../Results/priorization_tables/stroma_epiglandular.rds")

## 8. Visualización del output de NicheNet
# Expresión diferencial de los ligandos
top_ligand_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche

# Se guarda la tabla con las interacciones prioritarias para ambos grupos
receiver_oi <- c("Epi-Glandular_AS") 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(10, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

saveRDS(prioritized_tbl_oi, file = "../Results/df_stroma_epi/stroma_epiglandular_AS.rds")

receiver_oi <- c("Epi-Glandular_WOI") 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(10, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

saveRDS(prioritized_tbl_oi, file = "../Results/df_stroma_epi/stroma_epiglandular_WOI.rds")

## Stroma --> Epi_Glandular secretory
## 1. Definir los grupos de interés
niches <- list(
  "AS_niche" = list(
    "sender" = c("Stroma_AS"),
    "receiver" = c("Epi-Glandular secretory_AS")),
  "WOI_niche" = list(
    "sender" = c("Stroma_WOI"),
    "receiver" = c("Epi-Glandular secretory_WOI"))
  )

## 2. Calcular la expresión diferencial entre los grupos
assay_oi <- "SCT" 

DE_sender <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # Solo los ligandos importantes para los tipos celulares sender

DE_receiver <- calculate_niche_de(seurat_obj = seuratObj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # Solo los receptores ahora

DE_sender <- DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver <- DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))

# Procesar los resultados de la expresión diferencial:
expression_pct <- 0.10
DE_sender_processed <- process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed <- process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")

# Combinar la expresión diferencial de sender-receiver basado en los pares ligando-receptor:
specificity_score_LR_pairs <- "min_lfc"
DE_sender_receiver <- combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)

## 3. Calcular la expresión diferencial entre las diferentes regiones espaciales
include_spatial_info_sender <- FALSE 
include_spatial_info_receiver <- FALSE 

# Como no hay información espacial:
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
  spatial_info <- tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 

if(include_spatial_info_sender == TRUE){
  sender_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed <- process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  sender_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
} else {
  sender_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed <- sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  }

if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed <- process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
  
  receiver_spatial_DE_others <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
} else {
  receiver_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }

## 4. Calcular las actividades de los ligandos e inferir las relaciones ligando-target activas
lfc_cutoff <- 0.15 # recomendado para 10x como min_lfc cutoff 
specificity_score_targets <- "min_lfc"

DE_receiver_targets <- calculate_niche_de_targets(seurat_obj = seuratObj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 

DE_receiver_processed_targets <- process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)

background <- DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 <- DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
geneset_niche2 <- DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()

geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
length(geneset_niche1)
length(geneset_niche2)

top_n_target <- 250

niche_geneset_list <- list(
  "AS_niche" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "WOI_niche" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
  )

ligand_activities_targets <- get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)

## 5. Calcular la expresión escalada de los ligandos, receptores y targets en los tipos celulares de interés
features_oi <- union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)

dotplot <- suppressWarnings(Seurat::DotPlot(seuratObj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
exprs_tbl <- dotplot$data %>% as_tibble()
exprs_tbl <- exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

exprs_tbl_ligand <- exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
exprs_tbl_receptor <- exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
exprs_tbl_target <- exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)

exprs_tbl_ligand <- exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor <- exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

## 6. Fracción de expresión y receptor
exprs_sender_receiver <- lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 

## 7. Priorización de interacciones ligando-receptor y ligando-target
prioritizing_weights <- c("scaled_ligand_score" = 5,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 1,
                         "bona_fide" = 1)

output <- list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)

prioritization_tables <- get_prioritization_tables(output, prioritizing_weights)
saveRDS(prioritization_tables, file = "../Results/priorization_tables/stroma_epiglandular-secretory.rds")

## 8. Visualización del output de NicheNet
# Expresión diferencial de los ligandos
top_ligand_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
top_ligand_receptor_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)

ligand_prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche

# Se guarda la tabla con las interacciones prioritarias para ambos grupos
receiver_oi <- c("Epi-Glandular secretory_AS") 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(10, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

saveRDS(prioritized_tbl_oi, file = "../Results/df_stroma_epi/stroma_epiglandular-secretory_AS.rds")

receiver_oi <- c("Epi-Glandular secretory_WOI") 

filtered_ligands <- ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(10, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

saveRDS(prioritized_tbl_oi, file = "../Results/df_stroma_epi/stroma_epiglandular-secretory_WOI.rds")

## Cargar los dataframes del grupo AS
lumenal_AS <- readRDS("../Results/df_stroma_epi/stroma_epilumenal_AS.rds")
ciliated_AS <- readRDS("../Results/df_stroma_epi/stroma_epiciliated_AS.rds")
glandular_AS <- readRDS("../Results/df_stroma_epi/stroma_epiglandular_AS.rds")
glandular_secretory_AS <- readRDS("../Results/df_stroma_epi/stroma_epiglandular-secretory_AS.rds")

# Combinar los df en uno
combined_AS_df <- rbind(lumenal_AS, ciliated_AS, glandular_AS, glandular_secretory_AS)

# cargar las tablas prioritarias
pt_lumenal <- readRDS("../Results/priorization_tables/stroma_epilumenal.rds")
pt_ciliated <- readRDS("../Results/priorization_tables/stroma_epiciliated.rds")
pt_glandular <- readRDS("../Results/priorization_tables/stroma_epiglandular.rds")
pt_glandular_secretory <- readRDS("../Results/priorization_tables/stroma_epiglandular-secretory.rds")

# Coger solo las de ligando-receptor
pt_lumenal <- pt_lumenal$prioritization_tbl_ligand_receptor
pt_ciliated <- pt_ciliated$prioritization_tbl_ligand_receptor
pt_glandular <- pt_glandular$prioritization_tbl_ligand_receptor
pt_glandular_secretory <- pt_glandular_secretory$prioritization_tbl_ligand_receptor

# combinar las tablas en una
pt_combined <- rbind(pt_lumenal, pt_ciliated, pt_glandular, pt_glandular_secretory)

# log_FC mínimo comparado con el otro grupo, usando el df combinado
receiver_oi <- c("Epi-Lumenal_AS", "Epi-Ciliated_AS", "Epi-Glandular_AS", "Epi-Glandular secretory_AS")
lfc_plot <- make_ligand_receptor_lfc_plot(receiver_oi, combined_AS_df, pt_combined, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot$labels$x <- "LFC del receptor\n en Receptores"
lfc_plot$labels$fill <- "Receptor:\nLFC min vs\notros niches"
lfc_plot[[1]]$labels$y <- "Pares Ligando-Receptor priorizados"
lfc_plot[[1]]$labels$x <- "LFC del ligando\n en Emisores"
lfc_plot[[1]]$labels$fill <- "Ligando:\nLFC min vs\notros niches"
lfc_plot

## Cargar los dataframes del grupo WOI
lumenal_WOI <- readRDS("../Results/df_stroma_epi/stroma_epilumenal_WOI.rds")
ciliated_WOI <- readRDS("../Results/df_stroma_epi/stroma_epiciliated_WOI.rds")
glandular_WOI <- readRDS("../Results/df_stroma_epi/stroma_epiglandular_WOI.rds")
glandular_secretory_WOI <- readRDS("../Results/df_stroma_epi/stroma_epiglandular-secretory_WOI.rds")

# Combinar los df en uno
combined_WOI_df <- rbind(lumenal_WOI, ciliated_WOI, glandular_WOI, glandular_secretory_WOI)

# log_FC mínimo comparado con el otro grupo, usando el df combinado
receiver_oi <- c("Epi-Lumenal_WOI", "Epi-Ciliated_WOI", "Epi-Glandular_WOI", "Epi-Glandular secretory_WOI")
lfc_plot <- make_ligand_receptor_lfc_plot(receiver_oi, combined_WOI_df, pt_combined, plot_legend = FALSE, heights = NULL, widths = NULL)
lfc_plot$labels$x <- "LFC del receptor\n en Receptores"
lfc_plot$labels$fill <- "Receptor:\nLFC min vs\notros niches"
lfc_plot[[1]]$labels$y <- "Pares Ligando-Receptor priorizados"
lfc_plot[[1]]$labels$x <- "LFC del ligando\n en Emisores"
lfc_plot[[1]]$labels$fill <- "Ligando:\nLFC min vs\notros niches"
lfc_plot
