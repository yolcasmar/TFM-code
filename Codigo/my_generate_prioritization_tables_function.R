
# Modificación de la función scaling_zscore empleada por la función generate_prioritization_tables
my_scaling_zscore <- function(x) {
  if (typeof(x) == "double") {
    if (length(x) == 1) {
      return(0)
    }
    
    if (sum(!is.na(x)) > 1) {
      return((x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
    } else {
      return((x - mean(x, na.rm = TRUE)))
    }
  } else {
    return(x)
  }
}

# Modificación de la función generate_prioritization_tables, sustituyendo la función scaling_zscore por my_scaling_zscore
my_generate_prioritization_tables <- function (sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, 
                                               contrast_tbl, sender_receiver_tbl, grouping_tbl, prioritizing_weights = c(de_ligand = 1, 
                                                                                                                         de_receptor = 1, activity_scaled = 2, exprs_ligand = 2, 
                                                                                                                         exprs_receptor = 2, frac_exprs_ligand_receptor = 1, abund_sender = 0, 
                                                                                                                         abund_receiver = 0), fraction_cutoff, abundance_data_receiver, 
                                               abundance_data_sender) 
{
  requireNamespace("dplyr")
  receiver_receptor_prioritization = sender_receiver_de %>% 
    dplyr::ungroup() %>% dplyr::select(contrast, receiver, 
                                       receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>% 
    dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor) * 
                    lfc_receptor, p_val_receptor_adapted = -log10(p_val_receptor) * 
                    sign(lfc_receptor))
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% 
    dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, 
                                             ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, 
                                                                                                ties.method = "average", na.last = FALSE)), scaled_p_val_receptor = rank(desc(p_val_receptor), 
                                                                                                                                                                         ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), 
                                                                                                                                                                                                                            ties.method = "average", na.last = FALSE)), scaled_lfc_pval_receptor = rank(lfc_pval_receptor, 
                                                                                                                                                                                                                                                                                                        ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, 
                                                                                                                                                                                                                                                                                                                                                           ties.method = "average", na.last = FALSE)), scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, 
                                                                                                                                                                                                                                                                                                                                                                                                                                            ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  receiver_ligand_activity_prioritization_up = ligand_activities_targets_DEgenes$ligand_activities %>% 
    dplyr::ungroup() %>% dplyr::filter(direction_regulation == 
                                         "up") %>% dplyr::select(contrast, receiver, ligand, activity, 
                                                                 activity_scaled) %>% dplyr::rename(activity_up = activity, 
                                                                                                    activity_scaled_up = activity_scaled) %>% dplyr::distinct() %>% 
    dplyr::mutate(scaled_activity_scaled_up = scale_quantile_adapted(activity_scaled_up, 
                                                                     outlier_cutoff = 0.01), scaled_activity_up = scale_quantile_adapted(activity_up, 
                                                                                                                                         outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_up)
  receiver_ligand_activity_prioritization_down = ligand_activities_targets_DEgenes$ligand_activities %>% 
    dplyr::ungroup() %>% dplyr::filter(direction_regulation == 
                                         "down") %>% dplyr::select(contrast, receiver, ligand, 
                                                                   activity, activity_scaled) %>% dplyr::rename(activity_down = activity, 
                                                                                                                activity_scaled_down = activity_scaled) %>% dplyr::distinct() %>% 
    dplyr::mutate(scaled_activity_scaled_down = scale_quantile_adapted(activity_scaled_down, 
                                                                       outlier_cutoff = 0.01), scaled_activity_down = scale_quantile_adapted(activity_down, 
                                                                                                                                             outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_scaled_down)
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% 
    dplyr::select(contrast, sender, ligand, lfc_ligand, p_val_ligand) %>% 
    dplyr::distinct() %>% dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand) * 
                                          lfc_ligand, p_val_ligand_adapted = -log10(p_val_ligand) * 
                                          sign(lfc_ligand))
  sender_ligand_prioritization = sender_ligand_prioritization %>% 
    dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", 
                                           na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", 
                                                                     na.last = FALSE)), scaled_p_val_ligand = rank(desc(p_val_ligand), 
                                                                                                                   ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), 
                                                                                                                                                                      ties.method = "average", na.last = FALSE)), scaled_lfc_pval_ligand = rank(lfc_pval_ligand, 
                                                                                                                                                                                                                                                ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, 
                                                                                                                                                                                                                                                                                                   ties.method = "average", na.last = FALSE)), scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, 
                                                                                                                                                                                                                                                                                                                                                                                  ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, 
                                                                                                                                                                                                                                                                                                                                                                                                                                     ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_ligand)
  ligand_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% 
    dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% 
    dplyr::ungroup() %>% dplyr::select(group, sender, ligand, 
                                       avg_ligand_group) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% 
    dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand_group)) %>% 
    dplyr::arrange(-scaled_avg_exprs_ligand)
  ligand_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% 
    dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% 
    dplyr::ungroup() %>% dplyr::select(group, sender, ligand, 
                                       fraction_ligand_group) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% 
    dplyr::mutate(scaled_avg_frq_ligand = scale_quantile_adapted(fraction_ligand_group)) %>% 
    dplyr::arrange(-scaled_avg_frq_ligand)
  ligand_celltype_specificity_prioritization_pb = sender_receiver_info$pb_df_group %>% 
    dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% 
    dplyr::ungroup() %>% dplyr::select(group, sender, ligand, 
                                       pb_ligand_group) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>% 
    dplyr::mutate(scaled_pb_ligand = scale_quantile_adapted(pb_ligand_group)) %>% 
    dplyr::arrange(-scaled_pb_ligand)
  receptor_celltype_specificity_prioritization = sender_receiver_info$avg_df_group %>% 
    dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% 
    dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, 
                                       avg_receptor_group) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% 
    dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor_group)) %>% 
    dplyr::arrange(-scaled_avg_exprs_receptor)
  receptor_celltype_specificity_prioritization_frq = sender_receiver_info$frq_df_group %>% 
    dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% 
    dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, 
                                       fraction_receptor_group) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% 
    dplyr::mutate(scaled_avg_frq_receptor = scale_quantile_adapted(fraction_receptor_group)) %>% 
    dplyr::arrange(-scaled_avg_frq_receptor)
  receptor_celltype_specificity_prioritization_pb = sender_receiver_info$pb_df_group %>% 
    dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% 
    dplyr::ungroup() %>% dplyr::select(group, receiver, receptor, 
                                       pb_receptor_group) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>% 
    dplyr::mutate(scaled_pb_receptor = scale_quantile_adapted(pb_receptor_group)) %>% 
    dplyr::arrange(-scaled_pb_receptor)
  ligand_receptor_expressed_prioritization = sender_receiver_info$frq_df %>% 
    dplyr::inner_join(grouping_tbl) %>% dplyr::ungroup() %>% 
    dplyr::select(sample, group, sender, receiver, ligand, 
                  receptor, fraction_ligand, fraction_receptor) %>% 
    dplyr::distinct() %>% dplyr::group_by(ligand, receptor, 
                                          sender, receiver, group) %>% dplyr::summarise(n_samples = n(), 
                                                                                        n_expressing = sum(fraction_ligand > fraction_cutoff & 
                                                                                                             fraction_receptor > fraction_cutoff)) %>% dplyr::mutate(fraction_expressing_ligand_receptor = n_expressing/n_samples) %>% 
    dplyr::arrange(-fraction_expressing_ligand_receptor) %>% 
    dplyr::select(-n_samples, -n_expressing) %>% dplyr::ungroup()
  sender_abundance_prioritization = sender_receiver_info$rel_abundance_df %>% 
    dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% 
    dplyr::ungroup() %>% dplyr::select(group, sender, rel_abundance_scaled_sender) %>% 
    dplyr::distinct() %>% dplyr::arrange(-rel_abundance_scaled_sender)
  receiver_abundance_prioritization = sender_receiver_info$rel_abundance_df %>% 
    dplyr::inner_join(sender_receiver_tbl) %>% dplyr::inner_join(contrast_tbl) %>% 
    dplyr::ungroup() %>% dplyr::select(group, receiver, rel_abundance_scaled_receiver) %>% 
    dplyr::distinct() %>% dplyr::arrange(-rel_abundance_scaled_receiver)
  group_prioritization_tbl = contrast_tbl %>% dplyr::inner_join(sender_receiver_de) %>% 
    dplyr::inner_join(ligand_activities_targets_DEgenes$ligand_activities %>% 
                        dplyr::select(-target, -ligand_target_weight) %>% 
                        dplyr::distinct()) %>% dplyr::mutate(lr_interaction = paste(ligand, 
                                                                                    receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, 
                                                                                                                                       sender, receiver, sep = "_")) %>% dplyr::inner_join(sender_receiver_info$avg_df_group) %>% 
    dplyr::inner_join(sender_receiver_info$frq_df_group) %>% 
    dplyr::inner_join(sender_receiver_info$rel_abundance_df) %>% 
    dplyr::inner_join(sender_ligand_prioritization) %>% dplyr::inner_join(receiver_receptor_prioritization) %>% 
    dplyr::inner_join(receiver_ligand_activity_prioritization_up) %>% 
    dplyr::inner_join(receiver_ligand_activity_prioritization_down) %>% 
    dplyr::inner_join(ligand_celltype_specificity_prioritization) %>% 
    dplyr::inner_join(ligand_celltype_specificity_prioritization_frq) %>% 
    dplyr::inner_join(ligand_celltype_specificity_prioritization_pb) %>% 
    dplyr::inner_join(receptor_celltype_specificity_prioritization) %>% 
    dplyr::inner_join(receptor_celltype_specificity_prioritization_frq) %>% 
    dplyr::inner_join(receptor_celltype_specificity_prioritization_pb) %>% 
    dplyr::inner_join(ligand_receptor_expressed_prioritization) %>% 
    dplyr::inner_join(sender_abundance_prioritization) %>% 
    dplyr::inner_join(receiver_abundance_prioritization) %>% 
    mutate(max_scaled_activity = pmax(scaled_activity_scaled_up, 
                                      scaled_activity_scaled_down), na.rm = TRUE)
  sum_prioritization_weights = 2 * prioritizing_weights["de_ligand"] + 
    2 * prioritizing_weights["de_receptor"] + prioritizing_weights["activity_scaled"] + 
    prioritizing_weights["exprs_ligand"] + prioritizing_weights["exprs_receptor"] + 
    prioritizing_weights["frac_exprs_ligand_receptor"] + 
    prioritizing_weights["abund_sender"] + prioritizing_weights["abund_receiver"]
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::mutate(prioritization_score = ((prioritizing_weights["de_ligand"] * 
                                                                                                   scaled_lfc_ligand) + (prioritizing_weights["de_receptor"] * 
                                                                                                                           scaled_lfc_receptor) + (prioritizing_weights["de_ligand"] * 
                                                                                                                                                     scaled_p_val_ligand_adapted) + (prioritizing_weights["de_receptor"] * 
                                                                                                                                                                                       scaled_p_val_receptor_adapted) + (prioritizing_weights["activity_scaled"] * 
                                                                                                                                                                                                                           max_scaled_activity) + (prioritizing_weights["exprs_ligand"] * 
                                                                                                                                                                                                                                                     scaled_pb_ligand) + (prioritizing_weights["exprs_receptor"] * 
                                                                                                                                                                                                                                                                            scaled_pb_receptor) + (prioritizing_weights["frac_exprs_ligand_receptor"] * 
                                                                                                                                                                                                                                                                                                     fraction_expressing_ligand_receptor) + (prioritizing_weights["abund_sender"] * 
                                                                                                                                                                                                                                                                                                                                               rel_abundance_scaled_sender) + (prioritizing_weights["abund_receiver"] * 
                                                                                                                                                                                                                                                                                                                                                                                 rel_abundance_scaled_receiver)) * (1/sum_prioritization_weights)) %>% 
    dplyr::arrange(-prioritization_score)
  sample_prioritization_tbl = sender_receiver_info$avg_df %>% 
    dplyr::inner_join(sender_receiver_info$frq_df) %>% dplyr::inner_join(sender_receiver_info$pb_df) %>% 
    dplyr::inner_join(grouping_tbl) %>% dplyr::left_join(group_prioritization_tbl %>% 
                                                           dplyr::distinct(group, sender, receiver, ligand, receptor, 
                                                                           prioritization_score))
  sample_prioritization_tbl = sample_prioritization_tbl %>% 
    dplyr::mutate(lr_interaction = paste(ligand, receptor, 
                                         sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, 
                                                                                  sender, receiver, sep = "_"))
  sample_prioritization_tbl = sample_prioritization_tbl %>% 
    dplyr::group_by(id) %>% dplyr::mutate(scaled_LR_prod = my_scaling_zscore(ligand_receptor_prod), 
                                          scaled_LR_frac = my_scaling_zscore(ligand_receptor_fraction_prod), 
                                          scaled_LR_pb_prod = my_scaling_zscore(ligand_receptor_pb_prod)) %>% 
    dplyr::ungroup()
  sample_prioritization_tbl = sample_prioritization_tbl %>% 
    dplyr::left_join(abundance_data_receiver) %>% dplyr::left_join(abundance_data_sender)
  sample_prioritization_tbl$n_cells_sender[is.na(sample_prioritization_tbl$n_cells_sender)] = 0
  sample_prioritization_tbl$n_cells_receiver[is.na(sample_prioritization_tbl$n_cells_receiver)] = 0
  sample_prioritization_tbl$keep_sender[is.na(sample_prioritization_tbl$keep_sender)] = 0
  sample_prioritization_tbl$keep_receiver[is.na(sample_prioritization_tbl$keep_receiver)] = 0
  sample_prioritization_tbl = sample_prioritization_tbl %>% 
    dplyr::mutate(keep_sender_receiver = keep_receiver + 
                    keep_sender)
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 
                                                   0] = "Sender & Receiver absent"
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 
                                                   1 & sample_prioritization_tbl$keep_receiver == 0] = "Receiver absent"
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 
                                                   1 & sample_prioritization_tbl$keep_sender == 0] = "Sender absent"
  sample_prioritization_tbl$keep_sender_receiver[sample_prioritization_tbl$keep_sender_receiver == 
                                                   2] = "Sender & Receiver present"
  sample_prioritization_tbl = sample_prioritization_tbl %>% 
    mutate(keep_sender_receiver = factor(keep_sender_receiver, 
                                         levels = c("Sender & Receiver absent", "Receiver absent", 
                                                    "Sender absent", "Sender & Receiver present")))
  ligand_activities_target_de_tbl = ligand_activities_targets_DEgenes$ligand_activities %>% 
    dplyr::inner_join(ligand_activities_targets_DEgenes$de_genes_df %>% 
                        dplyr::rename(target = gene, p_val_adj = p_adj)) %>% 
    dplyr::select(contrast, receiver, ligand, activity, activity_scaled, 
                  target, ligand_target_weight, logFC, p_val, p_val_adj, 
                  direction_regulation) %>% dplyr::distinct()
  ligand_activities_target_de_tbl = ligand_activities_target_de_tbl %>% 
    dplyr::mutate(direction_regulation = factor(direction_regulation, 
                                                levels = c("up", "down")))
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::mutate(direction_regulation = factor(direction_regulation, 
                                                                                                      levels = c("up", "down")))
  group_prioritization_tbl = group_prioritization_tbl %>% dplyr::inner_join(group_prioritization_tbl %>% 
                                                                              dplyr::distinct(id, group, prioritization_score) %>% 
                                                                              dplyr::group_by(id) %>% dplyr::top_n(1, prioritization_score) %>% 
                                                                              dplyr::mutate(top_group = group) %>% dplyr::distinct(id, 
                                                                                                                                   top_group) %>% dplyr::ungroup())
  return(list(group_prioritization_tbl = group_prioritization_tbl, 
              sample_prioritization_tbl = sample_prioritization_tbl, 
              ligand_activities_target_de_tbl = ligand_activities_target_de_tbl))
}
