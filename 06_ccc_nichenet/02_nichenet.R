#!/usr/bin/env Rscript
## NicheNet ligand-activity prioritization: niche senders -> myofibroblast program/HES1
##
## IMPORTANT: lr_network$from is not a ligand ontology.  The NicheNet v2
## network contains intracellular processing proteins (for example PSEN1 and
## MFNG) on its `from` side.  Candidate ligands are therefore intersected with
## a high-confidence OmniPath intercell ligand annotation before scoring.
suppressMessages({library(nichenetr); library(dplyr); library(readr); library(tibble)})
OUT <- '/data/ssc/ccc_nichenet'
dir.create(file.path(OUT,'nichenet'), showWarnings=FALSE, recursive=TRUE)
log <- function(...) cat('[nichenet]', ..., '\n')

OMNIPATH_URL <- paste0(
  'https://omnipathdb.org/intercell?format=tsv',
  '&categories=ligand&scope=generic&aspect=functional',
  '&source=composite&entity_types=protein&transmitter=true&license=academic'
)
OMNIPATH_FILE <- Sys.getenv(
  'OMNIPATH_INTERCELL',
  file.path(OUT, 'reference', 'omnipath_ligand_composite.tsv')
)
OMNIPATH_MIN_SCORE <- as.integer(Sys.getenv('OMNIPATH_MIN_CONSENSUS_SCORE', '4'))
dir.create(dirname(OMNIPATH_FILE), showWarnings=FALSE, recursive=TRUE)
if (!file.exists(OMNIPATH_FILE)) {
  log('downloading OmniPath composite ligand annotation:', OMNIPATH_URL)
  download.file(OMNIPATH_URL, OMNIPATH_FILE, mode='wb', quiet=FALSE)
}

omni <- read_tsv(OMNIPATH_FILE, show_col_types=FALSE) %>%
  filter(
    category == 'ligand', parent == 'ligand', scope == 'generic',
    aspect == 'functional', source == 'composite', entity_type == 'protein',
    transmitter, consensus_score >= OMNIPATH_MIN_SCORE,
    secreted | plasma_membrane_transmembrane | plasma_membrane_peripheral,
    !is.na(genesymbol), genesymbol != ''
  ) %>%
  arrange(desc(consensus_score)) %>%
  distinct(genesymbol, .keep_all=TRUE)
bona_fide_ligands <- omni$genesymbol
write_csv(omni, file.path(OUT, 'nichenet', 'omnipath_bona_fide_ligands.csv'))
log('OmniPath bona fide ligands:', length(bona_fide_ligands),
    '| minimum consensus score:', OMNIPATH_MIN_SCORE)

ltm <- readRDS('/data/ssc/ref/nichenet_v2/ligand_target_matrix_nsga2r_final.rds')
lr  <- readRDS('/data/ssc/ref/nichenet_v2/lr_network_human_21122021.rds')
log('ligand_target_matrix dim:', paste(dim(ltm), collapse=' x '))

bg     <- readLines(file.path(OUT,'export/receiver_background.txt'))
gset   <- readLines(file.path(OUT,'export/geneset_myofib.txt'))
sender <- readLines(file.path(OUT,'export/sender_expressed.txt'))
bg <- bg[bg!='']; gset <- gset[gset!='']; sender <- sender[sender!='']
log('background genes:', length(bg), '| geneset:', length(gset), '| sender expressed:', length(sender))

## orientation: targets should be rownames, ligands colnames. Auto-fix if transposed.
if (length(intersect(gset, rownames(ltm))) < length(intersect(gset, colnames(ltm)))) {
  log('transposing ligand_target_matrix (targets were in columns)'); ltm <- t(ltm)
}
targets <- rownames(ltm); ligands_all <- colnames(ltm)

ligands_lr   <- unique(lr$from); receptors_lr <- unique(lr$to)
expressed_ligands   <- intersect(ligands_lr, sender)
expressed_receptors <- intersect(receptors_lr, bg)
potential_unfiltered <- lr %>%
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>% unique() %>%
  intersect(ligands_all)
potential <- intersect(potential_unfiltered, bona_fide_ligands)
background <- intersect(bg, targets)
geneset_oi <- intersect(gset, targets)
log('potential ligands before ligand-ontology filter:', length(potential_unfiltered),
    '| after filter:', length(potential),
    '| background in ltm:', length(background),
    '| geneset in ltm:', length(geneset_oi))
log('geneset_oi:', paste(geneset_oi, collapse=', '))

## Candidate-universe audit.  This makes every inclusion/exclusion reproducible
## and prevents ad hoc removal of individual intracellular proteins.
omni_audit <- read_tsv(OMNIPATH_FILE, show_col_types=FALSE) %>%
  filter(category == 'ligand', source == 'composite', !is.na(genesymbol), genesymbol != '') %>%
  arrange(desc(consensus_score)) %>%
  distinct(genesymbol, .keep_all=TRUE) %>%
  select(
    genesymbol, consensus_score, secreted,
    plasma_membrane_transmembrane, plasma_membrane_peripheral
  )
candidate_audit <- tibble(candidate=sort(unique(ligands_lr))) %>%
  mutate(
    sender_expressed = candidate %in% sender,
    has_expressed_receiver_lr = candidate %in% (lr %>%
      filter(to %in% expressed_receptors) %>% pull(from) %>% unique()),
    in_ligand_target_matrix = candidate %in% ligands_all
  ) %>%
  left_join(omni_audit, by=c('candidate'='genesymbol')) %>%
  mutate(
    omnipath_composite_ligand = !is.na(consensus_score),
    omnipath_score_pass = omnipath_composite_ligand &
      consensus_score >= OMNIPATH_MIN_SCORE,
    extracellular_topology = coalesce(secreted, FALSE) |
      coalesce(plasma_membrane_transmembrane, FALSE) |
      coalesce(plasma_membrane_peripheral, FALSE),
    included = sender_expressed & has_expressed_receiver_lr &
      in_ligand_target_matrix & omnipath_score_pass & extracellular_topology,
    exclusion_reason = case_when(
      !sender_expressed ~ 'not_expressed_in_sender_cells',
      !has_expressed_receiver_lr ~ 'no_edge_to_expressed_receiver_receptor',
      !in_ligand_target_matrix ~ 'absent_from_ligand_target_matrix',
      !omnipath_composite_ligand ~ 'not_an_OmniPath_composite_ligand',
      !omnipath_score_pass ~ paste0('OmniPath_consensus_score_below_', OMNIPATH_MIN_SCORE),
      !extracellular_topology ~ 'not_secreted_or_plasma_membrane',
      TRUE ~ 'included'
    )
  )
write_csv(candidate_audit, file.path(OUT, 'nichenet', 'candidate_ligand_audit.csv'))

notch_lig <- c('DLL1','DLL3','DLL4','JAG1','JAG2')
notch_rec <- c('NOTCH1','NOTCH2','NOTCH3','NOTCH4')
## which notch ligands survive the potential-ligand filter, and via which receptor
notch_pairs <- lr %>% filter(from %in% notch_lig & to %in% notch_rec)
notch_pairs$ligand_expressed_in_sender <- notch_pairs$from %in% expressed_ligands
notch_pairs$receptor_expressed_in_receiver <- notch_pairs$to %in% expressed_receptors
write_csv(notch_pairs, file.path(OUT,'nichenet/notch_lr_in_prior_network.csv'))
log('notch ligands in potential set:', paste(intersect(notch_lig, potential), collapse=', '))

la <- predict_ligand_activities(geneset=geneset_oi, background_expressed_genes=background,
        ligand_target_matrix=ltm, potential_ligands=potential)
la <- la %>% arrange(desc(aupr_corrected)) %>% mutate(rank=row_number(), n_ligands=nrow(la))
write_csv(la, file.path(OUT,'nichenet/ligand_activities.csv'))
la_notch <- la %>% filter(test_ligand %in% notch_lig)
write_csv(la_notch, file.path(OUT,'nichenet/ligand_activities_notch.csv'))
log('TOP20 ligands:'); print(head(la %>% select(test_ligand, aupr_corrected, rank), 20))
log('NOTCH ligand ranks:'); print(la_notch %>% select(test_ligand, aupr_corrected, rank, n_ligands))

## ligand-target links: top ligands + notch ligands -> program/HES1
best <- la %>% slice_max(aupr_corrected, n=30) %>% pull(test_ligand)
report_ligands <- union(best, intersect(notch_lig, potential))
lt <- lapply(report_ligands, function(lg){
  tryCatch(get_weighted_ligand_target_links(ligand=lg, geneset=geneset_oi,
            ligand_target_matrix=ltm, n=200), error=function(e) NULL)
}) %>% bind_rows()
write_csv(lt, file.path(OUT,'nichenet/ligand_target_links.csv'))
lt_notch <- lt %>% filter(ligand %in% notch_lig)
write_csv(lt_notch, file.path(OUT,'nichenet/notch_ligand_target_links.csv'))
log('notch ligand->target links (to myofib program/HES1):'); print(lt_notch)

st <- file(file.path(OUT,'nichenet/nichenet_status.txt'),'w')
writeLines(c(
    paste0('geneset_used=', length(geneset_oi)),
    paste0('background=', length(background)),
    paste0('potential_ligands_before_ontology_filter=', length(potential_unfiltered)),
    paste0('potential_ligands=', length(potential)),
    paste0('ligand_annotation=OmniPath_intercell_composite'),
    paste0('omnipath_min_consensus_score=', OMNIPATH_MIN_SCORE),
    paste0('ligand_topology=secreted_or_plasma_membrane'),
  paste0('notch_ligands_in_potential=', paste(intersect(notch_lig,potential),collapse=',')),
  paste0('metric=aupr_corrected'),
  paste0('top_ligand=', la$test_ligand[1]),
  paste0('best_notch_rank=', ifelse(nrow(la_notch)>0, min(la_notch$rank), NA)),
  paste0('best_notch_ligand=', ifelse(nrow(la_notch)>0, la_notch$test_ligand[which.min(la_notch$rank)], NA)),
  'NICHENET_OK'), st)
close(st)
log('NICHENET_OK')
