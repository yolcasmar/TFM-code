
# Paquetes necesarios
library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(stringr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)

# Cargar las redes de interacción ligando-receptor y la matriz ligando-objetivo actualizadas de la versión 2.0.0 de NicheNet
organism = "human"
if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS("../Data/multinichenet/ligand_target_matrix_nsga2r_final.rds")
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS("../Data/multinichenet/ligand_target_matrix_nsga2r_final_mouse.rds")
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
  }

# Cargar el objeto Seurat y pasar su ensayo de RNA a formato SingleCellExperiment
# Pasar los alias de los genes a sus símbolos y ponerlos en un formato de escritura que no ocasione errores en los pasos posteriores
seuratObj <- LoadH5Seurat("../Data/Single-cell/07_sc_AS_Estonia_Wang_Group_SCT_Main_lineages_clean.H5Seurat")
seuratObj <- DietSeurat(object = seuratObj, assays = "RNA", dimreducs = "umap")
seuratObj$Group_final <- gsub(pattern = "WOI WOI Control", replacement = "WOI", x = seuratObj$Group_final)
sce <- Seurat::as.SingleCellExperiment(seuratObj, assay = "RNA")
sce <- alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

# Pasar los nombres de las muestras y los tipos celulars a un formato de escritura correcto como en el paso anterior
# En este caso, los grupos a comparar ya están en un formato adecuado, pero si no fuera el caso habría que ajustarlo
sce$orig.ident <- sce$orig.ident %>% make.names()
sce$cell_type_zooming_v5 <- sce$cell_type_zooming_v5 %>% make.names()

# Definir las columnas donde se encuentran las muestras, los tipos celulares y los grupos a comparar en los datos
sample_id <- "orig.ident"
group_id <- "Group_final"
celltype_id <- "cell_type_zooming_v5"
covariates <- NA
batches <- NA

# Definir los tipos celulares a analizar como emisores y receptores (en este caso se analizan todos como emisores y receptores)
# Filtrar los datos para tener en cuenta solo los tipos celulares definidos
senders_oi <- c("Stroma", "Epi.Lumenal", "Epi.Ciliated", "Epi.Glandular", "Epi.Glandular.secretory")
receivers_oi <- c("Stroma", "Epi.Lumenal", "Epi.Ciliated", "Epi.Glandular", "Epi.Glandular.secretory")

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]

# Extraer información sobre la abundancia y la expresión de los tipos celulares emisores y receptores y unir la información de la expresión de los ligandos en las emisoras a los receptores en las receptoras
min_cells <- 10 # Recomendado por los creadores del paquete
abundance_expression_info <- get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, 
                                                           senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)

# Definir los contrastes y covarianzas de interés para el análisis de expresión diferencial
contrasts_oi <- c("'AS-WOI','WOI-AS'")

contrast_tbl <- tibble(
  contrast = c("AS-WOI","WOI-AS"), 
  group = c("AS","WOI"))

DE_info <- get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)

celltype_de <- DE_info$celltype_de$de_output_tidy

# Combinar la información de la expresión diferencial de ligandos-células emisoras y receptores-células receptoras
sender_receiver_de <- multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
  )

## Predecir las actividades de los ligandos de NicheNet y las uniones ligando-oobjetivo basándose en los resultados de la expresión diferencial
# Definir los parámetros para el análisis de la actividad de los ligandos de NicheNet
logFC_threshold <- 0.50
p_val_threshold <- 0.05
fraction_cutoff <- 0.05
p_val_adj <- FALSE
top_n_target <- 250

# El análisis de la actividad de los ligandos puede correrse en paralelo para cada tipo celular receptor 
verbose <- TRUE
cores_system <- 8
n.cores <- min(cores_system, sender_receiver_de$receiver %>% unique() %>% length()) # usa un procesador por tipo celular receptor

# Correr el análisis de actividad de ligando (este paso requiere algo de tiempo)
ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose, 
  n.cores = n.cores
  )))

## Usar la información recogida anteriormente para priorizar todos los pares célula emisora- ligando---célula receptora-receptor
# Definir los pesos prioritarios y preparar objetos agrupados
prioritizing_weights_DE <- c("de_ligand" = 1,
                             "de_receptor" = 1)

prioritizing_weights_activity <- c("activity_scaled" = 2)

prioritizing_weights_expression_specificity <- c("exprs_ligand" = 2,
                                                 "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency <- c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance <- c( "abund_sender" = 0,
                                              "abund_receiver" = 0)

prioritizing_weights <- c(prioritizing_weights_DE, 
                          prioritizing_weights_activity, 
                          prioritizing_weights_expression_specificity,
                          prioritizing_weights_expression_sufficiency, 
                          prioritizing_weights_relative_abundance)

# Hacer un dataframe con las muestras y los grupos
sender_receiver_tbl <- sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined <- SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl <- metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
  }

# Correr la priorización. La función original presentaba un error en el código, por lo que se ha utilizado una función modificada
prioritization_tables <- suppressMessages(my_generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
  ))

## Añadir información al conocimiento anterior y expresar la correlación entre ligando-receptor y la expresión del objetivo
lr_target_prior_cor <- lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, 
                                                     grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

# Guardar todo el output generado para poder usarlo en otra ocasión sin necesidad de volver a correr todos los pasos anteriores
path = "../Results/multinichenetr/"

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
  ) 

multinichenet_output = make_lite_output(multinichenet_output)

save = FALSE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output.rds"))
  }

## Visualización de los resultados del análisis de comunicación intercelular
# Gráfico circos de las uniones prioritarias, las 50 predicciones principales entre los contrastes, emisores y receptores de interés
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 25, rank_per_group = TRUE)

prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)

# Red de regulación intercelular. Las uniones son ligando-objetivo (uniones de regulación génica) y no interacciones proteina-proteina de ligando-receptor
prioritized_tbl_oi = get_top_n_lr_pairs(prioritization_tables, 25, rank_per_group = TRUE)

lr_target_prior_cor_filtered = prioritization_tables$group_prioritization_tbl$group %>% unique() %>% lapply(function(group_oi){
  lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi)
  lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
  lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50))
  lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)
}) %>% bind_rows()

colors_sender["Epi.Glandular.secretory"] = "pink" # El amarillo original no se distingue bien
graph_plot = make_ggraph_ligand_target_links(lr_target_prior_cor_filtered = lr_target_prior_cor_filtered, prioritized_tbl_oi = prioritized_tbl_oi, colors = colors_sender)
graph_plot$plot




############
lr_target_prior_cor_filtered = prioritization_tables$group_prioritization_tbl$group %>% unique() %>% lapply(function(group_oi){
  lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi)
  lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (target == "PGR" | target == "ESR1") & (pearson > 0.50 | spearman > 0.50) )
  lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (target == "PGR" | target == "ESR1") & (pearson > 0.50 | spearman > 0.50) )
  lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)
}) %>% bind_rows()

id_PGR_ESR1 <- (lr_target_prior_cor_filtered %>% filter(target == "PGR" | target == "ESR1")) %>% select(id) 

prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% lr_target_prior_cor_filtered$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) 


prioritized_tbl_oi <- prioritized_tbl_oi %>% filter(id %in% id_PGR_ESR1$id)
lr_target_prior_cor_filtered <- lr_target_prior_cor_filtered %>% filter(id %in% id_PGR_ESR1$id) %>% filter(target == "PGR" | target == "ESR1")

colors_sender["Epi.Glandular.secretory"] = "pink" # El amarillo original no se distingue bien
graph_plot = make_ggraph_ligand_target_links(lr_target_prior_cor_filtered = lr_target_prior_cor_filtered, prioritized_tbl_oi = prioritized_tbl_oi, colors = colors_sender)
graph_plot$plot