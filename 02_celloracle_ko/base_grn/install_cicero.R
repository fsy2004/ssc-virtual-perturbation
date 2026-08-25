options(Ncpus=6, repos='https://cloud.r-project.org')
BiocManager::install('Gviz', update=FALSE, ask=FALSE)
if(!requireNamespace('cicero',quietly=TRUE)) remotes::install_github('cole-trapnell-lab/cicero-release', ref='monocle3', upgrade='never')
cat('CICERO_DONE Gviz=',requireNamespace('Gviz',quietly=TRUE),' cicero=',requireNamespace('cicero',quietly=TRUE),'\n')
