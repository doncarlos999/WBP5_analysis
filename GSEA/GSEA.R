library(tidyverse)
gsea_jar <- file.path("~","gsea","gsea-3.0.jar")
gene <- "WBP5"
# IMPORTANT:
# If rnk file has quotes gmt file must also have quotes around genes
rnk_files <- Sys.glob(paste0("../DE/"gene,"*all_subs*.rnk"))
gmt_file <- "HOXa_B.gmt"
gmt_files <- gmt_file
working_dir <- getwd()
run_GSEA <- function(gmt_file,rnk_file){  # preranked list version
  #gmt_file <- gmt_files[1]
  #rnk_file <- rnk_files[1]
  analysis_name <- str_remove(rnk_file, "_GSEA.rnk")
  gmt_file_name <- str_remove(gmt_file, ".*/")
  gmt_file_name <- str_remove(gmt_file_name, ".gmt")
  command <- paste("java  -Xmx3G -cp",gsea_jar,
                   "xtools.gsea.GseaPreranked -gmx",
                   file.path(working_dir,gmt_file),
		   "-rnk" ,file.path(working_dir,rnk_file),
                   "-collapse false -nperm 1000 -permute gene_set -scoring_scheme weighted -rpt_label ",
                   paste0(analysis_name,"_",gmt_file_name),
                   "  -num 100 -plot_top_x 20 -rnd_seed 12345 -create_svgs true -set_max 500 -set_min 5 -zip_report false -out" ,
                   working_dir, "-gui false > gsea_output.txt",sep=" ")
  system(command)
}
for (i in seq_along(gmt_files)){
	for (j in seq_along(rnk_files)){
	  cat("\nrunning GSEA on ", gmt_files[i], " with rnk file ",rnk_files[j])
	  run_GSEA(gmt_files[i], rnk_files[j])
		   }
}
system(paste0("rm -fr ",gene,"_GSEA"))
system(paste0("mkdir ",gene,"_GSEA"))
system(paste0("mv ",gene,"*GseaPreranked* ",gene,"_GSEA"))

