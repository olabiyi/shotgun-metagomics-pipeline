library(ggpubr)
library(tools)
library(phyloseq)
library(DescTools)
library(vegan)
library(gtools)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(tidyverse)
library(DESeq2)
library(tidyverse)
library(sva)
library(RColorBrewer)
library(pheatmap)
#library(beepr)
library(vsn)
library(gridExtra)
library(openxlsx)
library(plotly)
library(RColorBrewer)
library(colorspace)
library(BiocParallel)
library(ggrepel)
library(DT)
library(clusterProfiler)
library(edgeR)
library(metagenomeSeq)
library(xlsx) #for writing to excel files
# make sure that you are using the correct architecture of java
#i.e 64 bit java with 64 bit R version or 32 bit java with 32 bit R version
library("rJava") #load rJava
library(venn)
library(patchwork)
register(MulticoreParam(30))
#register(SnowParam(30))
library(GO.db)
library(heatmaply)
library(pathview)
library(FSA)

publication_format <- theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.ticks.length=unit(-0.15, "cm"),
        axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
        axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")), 
        axis.title = element_text(size = 18,face ='bold.italic', color = 'black'), 
        axis.text = element_text(size = 16,face ='bold', color = 'black'),
        legend.position = 'right', legend.title = element_text(size = 15,face ='bold', color = 'black'),
        legend.text = element_text(size = 14,face ='bold', color = 'black'),
        strip.text =  element_text(size = 14,face ='bold', color = 'black'))


custom_palette <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00",
                    "#CAB2D6","#6A3D9A","#FF00FFFF","#B15928","#000000","#FFC0CBFF","#8B864EFF","#F0027F",
                    "#666666","#1B9E77", "#E6AB02","#A6761D","#FFFF00FF","#FFFF99","#00FFFFFF",
                    "#B2182B","#FDDBC7","#D1E5F0","#CC0033","#FF00CC","#330033",
                    "#999933","#FF9933","#FFFAFAFF",colors()) 

# remove white colors
custom_palette <- custom_palette[-c(21:23,
                                    grep(pattern = "white|snow", x = custom_palette, ignore.case = TRUE))]

get_sample_pathway_KO <- function(count_func_profile, sample_name, pathway_kos){
  
  pathway_table <- count_func_profile[pathway_kos,]
  pathway_table <- pathway_table[complete.cases(pathway_table),]
  sample_counts <-  pathway_table[,sample_name] > 0
  pathway_table <- pathway_table[sample_counts,,drop=FALSE]
  return(rownames(pathway_table))
} 

# vector2dataframe <- function(named_vector, col1_name, col2_name){
#   
#   # Example
#   # denitrification <- 
#   #   c(
#   #     'K00370'="nitrate reductase / nitrite oxidoreductase, alpha subunit [EC:1.7.5.1 1.7.99.-] | (GenBank) narG; nitrate reductase 1, alpha subunit",
#   #     'K00371'="nitrate reductase / nitrite oxidoreductase, beta subunit [EC:1.7.5.1 1.7.99.-] | (GenBank) narH; nitrate reductase 1, beta (Fe-S) subunit",
#   #     'K00374'="nitrate reductase gamma subunit [EC:1.7.5.1 1.7.99.-] | (GenBank) narI; nitrate reductase 1, gamma (cytochrome b(NR)) subunit",
#   #     'K02567'="nitrate reductase (cytochrome) [EC:1.9.6.1] | (GenBank) napA; nitrate reductase, periplasmic, large subunit",
#   #     'K02568'="K02568 nitrate reductase (cytochrome), electron transfer subunit | (GenBank) napB; nitrate reductase, small, cytochrome C550 subunit, periplasmic",
#   #     'K00368'="nitrite reductase (NO-forming) [EC:1.7.2.1] | (GenBank) aniA; Copper-containing nitrite reductase precursor",
#   #     'K15864'="nitrite reductase (NO-forming) / hydroxylamine reductase [EC:1.7.2.1 1.7.99.1] | (GenBank) nitrite reductase",
#   #     'K04561'="nitric oxide reductase subunit B [EC:1.7.2.5] | (GenBank) cytochrome c oxidase-like protein, subunit I",
#   #     'K02305'="nitric oxide reductase subunit C | (GenBank) cytochrome c",
#   #     'K00376'="nitrous-oxide reductase [EC:1.7.2.4] | (GenBank) nosZ; TAT-dependent nitrous-oxide reductase"
#   #   )
#   # 
#   # denitrification_Ko2name <- vector2dataframe(denitrification, "ID", "Definition")
#   
#   df <- data.frame(col1=names(named_vector), col2=named_vector)
#   colnames(df) <- c(col1_name,col2_name)
#   return(df)
# }

plot_pca <- function(mat, metadata,  axis1 = "PC1", axis2 = "PC2", color_factor,
                     shape_factor = NULL, sample_column='Sample.ID', distance_method="euclidean"){
  # mat - matrix with samples as rows and column as features
  
  samples <- rownames(metadata)
  mat <- mat[samples,]
  
  mat_pca <- prcomp(mat)
  mat_pca_percentVar <- mat_pca$sdev^2/sum(mat_pca$sdev^2)
  names(mat_pca_percentVar) <- colnames(mat_pca$x)
  
  if(!is.null(shape_factor)){
    
    pca2plot <- data.frame(mat_pca$x[,c(axis1,axis2)],
                           Color=metadata[,color_factor],
                           Shape=metadata[,shape_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- ggplot(pca2plot,aes(x=get(axis1),
                             y=get(axis2),
                             col=Color,
                             shape=Shape)) +
      geom_point(size=5) +
      # geom_label_repel(aes(label=Sample)) +
      # geom_path(aes(color=Type,
      #               group=Type)) +
      xlab(sprintf("%s: %1.2f%% of variance", axis1, mat_pca_percentVar[axis1]*100)) + 
      ylab(sprintf("%s: %1.2f%% of variance", axis2, mat_pca_percentVar[axis2]*100)) + 
      labs(color=color_factor, shape=shape_factor) +
      coord_fixed()
    
  }else{
    
    pca2plot <- data.frame(mat_pca$x[,c(axis1,axis2)],
                           Color=metadata[,color_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    p <- ggplot(pca2plot,aes(x=get(axis1),
                             y=get(axis2),
                             col=Color)) +
      geom_point(size=5) +
      # geom_label_repel(aes(label=Sample)) +
      # geom_path(aes(color=Type,
      #               group=Type)) +
      xlab(sprintf("%s: %1.2f%% of variance", axis1, mat_pca_percentVar[axis1]*100)) + 
      ylab(sprintf("%s: %1.2f%% of variance", axis2, mat_pca_percentVar[axis2]*100)) + 
      labs(color=color_factor, shape=shape_factor) +
      coord_fixed()
    
  }
  
  return(p)
  
}

# A funcion to convert a named vector to a dataframe sorted in alphabetical order 
vector2dataframe <- function(named_vector, column_name){
  # named_vector - A named vector of compacter letters with the groups as the name and letters as the value
  
  sorted_group_names <- sort(names(named_vector))
  df <- data.frame(group=named_vector[sorted_group_names])
  colnames(df) <- column_name 
  rownames(df) <- sorted_group_names
  return(df)
  
}



remove_rare_features <- function(feature_table, cut_off_percent=3/4, top=100, by_sample=TRUE){
  # feature_table - feature table matrix with samples as columns and features as rows
  #cut_off_percent - cut-off fraction  or decimal between 0.001 to 1 
  #of the total number of samples for determining the most abundant features
  # by default removes features that are not present in 3/4 of the total number of samples
  # top - an integer specifying the maximum number of abundant features to retian after filtering
  # by_sample - A boolean set to TRU or FALSE whether features should be removed based on
  #  there presences in less the then cut_off_perecent. If set to FALSE, features will be dropped
  # based on relative abundance at the set threshold /cut_off_percent
  # all also returns a maximum of top (by default 100) abundant features
  # remove rare OTUs/ASVs
  
  if(by_sample){
    # Filter by occurrence in a fraction of samples
    # Define a cut-off for determining what's rare
    cut_off <- cut_off_percent * ncol(feature_table)
    # Get the occurrence for each feature
    feature_occurence <- rowSums(feature_table > 0)
    # Get the names of the abundant features
    abund_features <- names(feature_occurence[feature_occurence >= cut_off])
  }else{
    # filter based on relative abundance
    abund_features <- apply(X = feature_table, MARGIN = 1, FUN = function(x) max(x) > cut_off_percent)
  }
  
  
  abun_features.m <- feature_table[abund_features,]
  # drop excess features
  if(nrow(abun_features.m) > top){
    abund_features <- names(abun_features.m %>% rowMeans() %>% sort(decreasing = TRUE))[1:top]
    # This is to ensure that my edge correlation function does not break
    # As it assumes that the ASV are in order
    abund_features <- DescTools::SortMixed(abund_features)
    abun_features.m <- feature_table[abund_features,]
  }
  return(abun_features.m)
}

create_contrast_df <- function(dds,compar_print_name, all_comparisons){
  
  compar2use <-  all_comparisons[[compar_print_name]]
  res_df <- results(dds, contrast=compar2use, alpha = DESEQ_PADJ_CUTOFF) %>% 
    as.data.frame %>% 
    as.tibble %>% 
    # add a linear fold change column by raising 2 to the powere of the log2folchange vaue
    mutate(linearFC = ifelse(is.na(log2FoldChange),
                             yes = NA,
                             no = ifelse(log2FoldChange>0,
                                         yes = 2^log2FoldChange,
                                         no = -1/(2^log2FoldChange))%>% 
                               signif(digits = 3))) %>% 
    mutate(pvalue = ifelse(test = is.na(pvalue),
                           yes = NA,
                           no = signif(pvalue,digits = 3))) %>% 
    mutate(padj = ifelse(test = is.na(padj),
                         yes = NA,
                         no = signif(padj,digits = 3))) %>% 
    mutate(pass = ifelse(test = abs(as.numeric(linearFC)) >= LINEAR_FC_CUTOFF & 
                           padj <= PADJ_CUTOFF & 
                           !is.na(padj),
                         yes = ifelse(test = as.numeric(linearFC)>0,
                                      yes="up",
                                      no="down"),
                         no = "")) %>%
    dplyr::select(linearFC,pvalue,padj,pass) %>%
    # select(linearFC,pvalue,padj,pass) %>% 
    rename(!!paste0("linearFC.",compar_print_name) := linearFC) %>%
    rename(!!paste0("pvalue.",compar_print_name) := pvalue) %>%
    rename(!!paste0("padj.",compar_print_name) := padj) %>%
    rename(!!paste0("pass.",compar_print_name) := pass) %>%
    as.data.frame(check.names = FALSE)
  
  return(res_df)
}

create_results_dataframe <- function(dds, norm_counts, all_comparisons,
                                     test_name="KO", col_data=no_control_metadata){
  
  #test <- "KO"
  ### Add preamble: gene name and counts
  res_dim2 <- 5 # Columns per contrast in res_df: linearFC, pvalue, padj and pass
  
  res_df_initial <- data.frame(gene = rownames(rowData(dds)),
                               norm_counts, check.names = FALSE)
  
  # Save special header info:
  res_df_compar_head = c(test_name,rep("",nrow(col_data)))  
  # Save widths of data groups:
  res_df_grouping <- c(1,ncol(norm_counts))
  
  ### Add statistical tables
  res_df <-  map_dfc(.x = names(all_comparisons),
                     .f = create_contrast_df, all_comparisons=all_comparisons, dds=dds) %>% 
    as.data.frame(check.names = FALSE, StringsAsFactor=FALSE)
  
  res_df <- cbind(res_df_initial,res_df)
  
  for ( compar_print_name  in names(all_comparisons)){
    
    res_df_compar_head <- c(res_df_compar_head,
                            rep(compar_print_name,
                                res_dim2))
    # Add width of group to res_df_grouping
    res_df_grouping <- c(res_df_grouping,res_dim2)
    
  } 
  
  res_df1=list(res_df=res_df,
               res_df_grouping=res_df_grouping,
               res_df_compar_head=res_df_compar_head)
  
  
  ## Additional columns - Add pass any column
  
  ### Add pass_any field aka significant fields
  # pass_any = TRUE if any contrast pass value is 1
  # Do only if there are more than one 'pass' column (not in Scales)
  if (sum(str_detect(string = names(res_df),
                     pattern = "pass")) > 1 ) {
    res_df$pass_any <-
      apply(X = res_df[,str_detect(string = names(res_df),
                                   pattern = "pass")],
            MARGIN = 1,
            FUN = function(x) any(unlist(x)!=""))
    # Convert TRUE/FALSE to 1/"":
    res_df$pass_any <- ifelse(test = !is.na(res_df$pass_any) & res_df$pass_any==TRUE,
                              yes = 1,
                              no = "")
  } else {
    res_df$pass_any <- ifelse(test = res_df[,str_detect(string = names(res_df),pattern = "pass")]!="",
                              yes = 1,
                              no = "")
  }
  
  res_df_grouping <- c(res_df_grouping,1)
  res_df_compar_head <- c(res_df_compar_head, "")
  
  final <- list(res_df=res_df, res_df1=res_df1)
  return(final)
}

define_contrasts <- function(factor_column, col_data=no_control_metadata){
  
  # Get the unique levels of the factor column
  levels <- unique(as.character(col_data[,factor_column]))
  
  # Get all pairwise combinations of the factor
  combinations <- utils::combn(levels,2)
  
  # Paste the factor column name and every element in each column separated by commas
  comparisons <- apply(X = combinations, MARGIN = 2, FUN = function(column) paste0(factor_column,',',paste(column,collapse = ",")))
  
  comparison_names <-  apply(X = combinations, MARGIN = 2, FUN = function(column) paste(column,collapse = ".vs.")) 
  
  all_comparisons <- str_split(string = comparisons,pattern = ",")
  
  names(all_comparisons) <- comparison_names
  
  return(all_comparisons)
}

run_deseq2_with_vsd_normalization <- function(Design, count_func_profile, metadata){
  
  common.ids <- intersect(colnames(count_func_profile),rownames(metadata))
  metadata <- metadata[common.ids,]
  count_func_profile <- count_func_profile[,common.ids]
  dds0 <- DESeqDataSetFromMatrix(countData = count_func_profile,
                                 colData = metadata,
                                 design =  Design)
  
  dds <- DESeq(dds0,
               parallel = TRUE,
               betaPrior = TRUE,  
               modelMatrixType = "expanded")
  
  # Deseq2 normalization
  norm_counts <- counts(dds,normalized=TRUE)
  
  vsd <- varianceStabilizingTransformation(dds,blind=TRUE)
  
  final <- list(dds0=dds0, dds=dds, norm_counts=norm_counts, vsd=vsd)
  
  return(final)
}

make_volcano_plots <- function(dds, all_comparisons, plot_dir, 
                               plot_MA=FALSE, category="All", 
                               save_volcano=FALSE, 
                               LOG_FC_CUTOFF=0.3785116,
                               DESEQ_PADJ_CUTOFF=0.1){
  
  volcano_plots <- list()
  ### Drawing MA plots and volcano plots of contrasts
  for ( compar_print_name  in names(all_comparisons)){
    compar2use <-  all_comparisons[[compar_print_name]]
    # Get results:
    compar_res <- results(dds,
                          contrast=compar2use,
                          alpha = DESEQ_PADJ_CUTOFF)
    if(plot_MA){
      # Draw MA plot:
      png(filename = sprintf("%s/MAPlot_%s.png",plot_dir,compar_print_name))
      
      plotMA(compar_res,
             alpha = PADJ_CUTOFF,
             main=paste0("MAplot ",compar_print_name),
             ylim = range(compar_res$log2FoldChange,na.rm = T)
      )
      abline(h = c(LOG_FC_CUTOFF,-LOG_FC_CUTOFF),col="green")
      dev.off()
    }
    
    # Draw volcano plot:
    data2plot <- compar_res %>% 
      as_tibble %>% 
      filter(!is.na(padj) & padj != 1) %>%
      mutate(highFC=(abs(log2FoldChange)>LOG_FC_CUTOFF),
             lowPADJ=(padj<=PADJ_CUTOFF),
             Significant=((abs(log2FoldChange)>LOG_FC_CUTOFF) & (padj<=PADJ_CUTOFF)),
             neglog10=-log10(pvalue)) 
    
    volcano_plot <- 
      ggplot(data2plot, aes(x=log2FoldChange, y=neglog10)) +
      geom_point(aes(colour = Significant), 
                 size=2.5) +
      scale_colour_manual(values = c("TRUE"= "red", "FALSE"= "black")) + 
      xlab(expression(log[2]~fold~change)) +
      ylab(expression(-log[10]~pvalue)) +
      geom_vline(xintercept = c(-LOG_FC_CUTOFF, LOG_FC_CUTOFF),color = "blue", lty = 2) +
      publication_format  
    if(save_volcano){
      ggsave(plot = volcano_plot, dpi = 600,
             filename = sprintf("%s/%s_VolcanoPlot_%s.png",plot_dir,category,compar_print_name))
    }
    volcano_plots[[compar_print_name]] <-  volcano_plot + labs(x=NULL, y=NULL) + ggtitle(compar_print_name)
    
  }
  return(volcano_plots)
}

combine_plots <- function(graphs,Title,Xlab,Ylab=NULL,Legend='right', ncols=NULL){
  
# Get none null graphs and remove the legends and x labs
graphs <- map(graphs[lengths(graphs)>0],
              .f = function(x){ x + theme(legend.position = Legend , axis.title.x = element_blank())} )

if(!is.null(ncols)){
  
  p <- wrap_plots(graphs, ncol = ncols,
                  guides = 'collect')
}else{
  
p <- wrap_plots(graphs, ncol = if_else(condition = length(graphs) < 3,
                                       true = 1,false = 2),
                guides = 'collect')
}
if(!is.null(Ylab)){
  
  p <- (wrap_elements(grid::textGrob(Ylab, rot = 90, gp=grid::gpar(fontsize=16,fontface='bold'))) + p) + 
    plot_layout(widths = c(1,30))
  
}

# combine plots
p <- wrap_elements(grid::textGrob(Title, x=unit(0.53,"npc"), gp=grid::gpar(fontsize=18,fontface='bold'))) / 
  p / 
  wrap_elements(grid::textGrob(Xlab, x=unit(0.53,"npc"), gp=grid::gpar(fontsize=16,fontface='bold'))) + 
  plot_layout(heights = c(1,30,1))



p

}

combine_enrichement_dfs <- function(enrich_dfs,GO_df){

# Get only dataframes with enrichment results    
dfs <- Filter(function(x) nrow(x) > 0, enrich_dfs)

# Add comparison and Go category columns and re-arrange the columns
sig_enrich_df <- map2_dfr(.x = dfs ,
                             .y = names(dfs), 
                             .f = function(x,y){
                               
                               if(!is_empty(x) && nrow(x) > 0){
                                 x[,'GO_Category'] <- GO_df[GO_df$go_id %in% x$ID,"ontology"]
                                 x[,'Comparison'] <- y
                                 x %>% dplyr::select(Comparison,ID,GO_Category,everything())
                               }
                               
                               
                             })

return(sig_enrich_df )
}


# A function to return whether a gene/feature between two comparisons is 
# up or down regulated relative to the first group
regulation_table <- function(sig_table, group1, group2, suffix="_mean", high='up', low='down'){
  
  ROWNAMES <- rownames(sig_table)
  comp_mean_table <- sig_table %>% 
    dplyr::select(contains(suffix)) %>% 
    dplyr::select(matches(sprintf("%s%s|%s%s", group1,suffix,group2,suffix))) %>% 
    as.data.frame(check.names=FALSE)
  
  original_colnames <- colnames(comp_mean_table)
  colnames(comp_mean_table) <- c("group1","group2")
  rownames(comp_mean_table) <- ROWNAMES
  
  comp_mean_table <- comp_mean_table %>%
    mutate(regulated=case_when(
      (group1 > group2) ~ high,
      (group1 < group2) ~ low,
      TRUE ~ "equal")
    )
  rownames(comp_mean_table) <- ROWNAMES
  colnames(comp_mean_table)[1:2] <- original_colnames
  return(comp_mean_table)
}


enrichment_analysis <- function(res_df, comparison_names, TERM2GENE=sub_GO2GENE,
                                TERM2NAME=sub_GO2NAME,
                                filter_genes_by= c('up', 'down'), 
                                DE_stat_method="deseq", background=NULL){
  
  if(DE_stat_method != "deseq"){
    regulated <- 'regulated'
  }
  
  if(all(sort(filter_genes_by) == sort(c('down','up')))){
    
    if(DE_stat_method == "deseq"){
    cond <- "pass.{comparison} %in% c('up', 'down') "
    }else{
      cond <- "{regulated} %in% c('up', 'down') "
    }
        #plot_title <- "All significant genes GO enrichment for {gsub(pattern='_|\\\\.',replacement= ' ' , x=comparison)}"
 
   }else if(filter_genes_by == 'up'){
     
     if(DE_stat_method == "deseq"){
    cond <- "pass.{comparison} %in% c('up') "
    #plot_title <- "{toTitleCase(filter_genes_by)} regulated GO enrichment for {comparison}"
     }else{
    cond <- "{regulated} %in% c('up') "    
     }
     
  }else{
    
    if(DE_stat_method == "deseq"){
    cond <- "pass.{comparison} %in% c('down') "
    #plot_title <- "{toTitleCase(filter_genes_by)} regulated GO enrichment for {comparison}"
    }else{
    cond <- "{regulated} %in% c('down') " 
    }
    
  }
  
  if(DE_stat_method == "deseq"){
  plot_title <- "{gsub(pattern='_|\\\\.',replacement= ' ' , x=comparison)  %>% toTitleCase}"
  }else{
    plot_title <- "{gsub(pattern='\\\\-',replacement= ' vs ' , x=comparison)  %>% toTitleCase}"
  }
  
  dfs <- vector(mode = 'list', length = length(comparison_names))
  graphs <- vector(mode = 'list', length = length(comparison_names))
  names(dfs) <- comparison_names
  names(graphs) <- comparison_names
  #####dplyr::rename(comp_col= comparison) %>%   
  for (comparison in comparison_names) {
    
    if(DE_stat_method == "deseq"){
    sig_df <- res_df %>% 
      dplyr::select(gene, grep(x = colnames(res_df),pattern=comparison)) 
    sig_df <- sig_df[complete.cases(sig_df),] %>% filter_(glue::glue(cond))
    sig_genes <- sig_df[, "gene"]
    background <- res_df[,"gene"]
    }else{
     
      if(!is.null(background)){
        sig_df <- res_df %>% 
          filter(comp_col == comparison) %>% 
          dplyr::select(gene, regulated) %>% 
          filter_(glue::glue(cond))

        sig_genes <- sig_df[, "gene"]
      }else{
        stop("You must provide a vector of gene names to serve as background for enrichment test")
      }
      
    }
    
    if(is_empty(sig_genes)){ next }
    # GO enrichment
    enriched_genes <- enricher(gene = as.character(sig_genes),
                               universe = as.character(background), 
                               TERM2GENE = TERM2GENE, 
                               TERM2NAME = TERM2NAME)
    
    #enriched_df <- as.data.frame(enriched_genes)
    dfs[[comparison]] <- as.data.frame(enriched_genes)
    # plot and display only enriched GO/pathway terms
    if(nrow(dfs[[comparison]]) > 0){
      
      graphs[[comparison]] <- clusterProfiler::dotplot(enriched_genes, showCategory=30) + 
                        ggtitle(glue::glue(plot_title))
      
      
    }
  
    }
  
  return(list(dfs=dfs,graphs=graphs))
}


#function to calculate standard error
SE<- function(x){
  sd(x)/sqrt(length(x))
}

# Retrieve a samples character vector order by the levels in group column of a metadata
order_samples <- function(levels_order, metadata, metadata_column='treatment'){
  map(levels_order,.f = function(level) rownames(metadata[metadata[,metadata_column] == level,])) %>% 
    flatten_chr()
}

# run permanova on a matrix using a choice distance metric, euclidean by default.
permanova <- function(mat, y_vect=NULL, dist_method="euclidean", formula=NULL, metadata=NULL){
  # EXAMPLE: permanova(vsd.df,metadata$treatment)
  
  if(!is.null(formula) && !is.null(metadata)){
    d <- vegan::vegdist(mat, method = dist_method)
    adonis2(formula,data = metadata)
  }else{
    d <- vegan::vegdist(mat, method = dist_method)
    adonis(d  ~ y_vect)
  }
  
}

mapGO2KO <- function(KO2GO_path="data/function_tables/KO2GO.txt", output_dir='data/function_tables/') {
  # function to map GO ids to KO ids
  # In order to get a mapping of KO to GO terms, I did the following
  #1. Navigate to https://www.genome.jp/linkdb/ on your browser
  #	This site is very useful for getting mapping between various databases
  #2. Click on the GO button in the seachable databases in DBGET box
  #3. Set "from" textbox to orthology and the "to" textbox to GO
  #4. Click download
  #5. rename the file from orthology_go.list to KO2GO.txt
  #6. use this regular expression "\toriginal$" using sed or in notepad
  # to remove the last column containing the word original
  
  KO2GO <-  read.table(file = KO2GO_path, col.names = c('KO','GO'))
  KO2GO <- apply(KO2GO,2,toupper) %>% as.data.frame()
  GO_ids <- unique(KO2GO$GO)
  names(GO_ids)  <- GO_ids 
  GO2name <- map_dfr(.x = GO_ids, .f = function(id) {
    id_info <- GO.db::GOTERM[[id]]
    if(!is.null(id_info)){
      return(data.frame(id=id,
                        term=AnnotationDbi::Term(id_info), 
                        ontology=AnnotationDbi::Ontology(id_info)
      ))
    }
  })
  
  TERM2NAME <- GO2name[,c('id','term')]
  TERM2GENE <- KO2GO[,c('GO', 'KO')]
  TERM2GENE[,'KO'] <- gsub(pattern = 'KO:', replacement = '', TERM2GENE[,'KO'])
  
  GO2GENE <- clusterProfiler::buildGOmap(TERM2GENE)
  GO2NAME <- clusterProfiler::go2term(unique(GO2GENE$GO))
  GO2ONTOLOGY <- clusterProfiler::go2ont(unique(GO2GENE$GO))
  GO_df <- cbind(GO2NAME,ontology=GO2ONTOLOGY$Ontology)
  
  GO_CATEGORIES <- c(BP='BP',MF='MF',CC='CC') 
  sub_GO_dfs <- map(.x = GO_CATEGORIES, 
                    .f =  function(category){
                      subset(x = GO_df, subset = ontology == category)
                    })
  
  # Subset the genes/KO per category
  sub_GO2GENE <- map(.x = GO_CATEGORIES, 
                    .f =  function(category){
                      subset(x = GO2GENE, subset = GO %in% sub_GO_dfs[[category]]$go_id)
                    })
  # Subset the GO names per category
  sub_GO2NAME <- map(.x = GO_CATEGORIES, 
                     .f =  function(category){
                       subset(x = GO2NAME, subset = go_id %in% sub_GO_dfs[[category]]$go_id)
                     })
  
  
  data2write <- list(GO2GENE=GO2GENE,
                     BP_GO2GENE= sub_GO2GENE$BP,
                     MF_GO2GENE= sub_GO2GENE$MF,
                     CC_GO2GENE= sub_GO2GENE$CC,
                     GO2NAME=GO2NAME,
                     BP_GO2NAME= sub_GO2NAME$BP,
                     MF_GO2NAME= sub_GO2NAME$MF,
                     CC_GO2NAME= sub_GO2NAME$CC,                     
                     GO_df=GO_df)
  
  file_names <- names(data2write)
  
  walk2(.x = data2write, .y = file_names, 
        .f = function(.x,.y){
          
          write.table(x = .x, 
                      file = sprintf('%s/%s.tsv', 
                                     output_dir,.y), 
                      sep = "\t", 
                      row.names = FALSE, col.names = TRUE, quote = FALSE)
          
        })
}

pathway2ko <- function(output_dir){
  # link pathway to KO id
  pathway2ko<- KEGGREST::keggLink("ko",'pathway')
  pathway2ko <- data.frame(pathway=names(pathway2ko), KO=pathway2ko)
  pathway2ko <- pathway2ko %>%
    filter(str_detect(string = pathway,pattern = 'map')) %>%
    mutate(KO=str_replace(KO, pattern = "ko:",replacement =  ''))
  write.table(x = pathway2ko, 
              file = sprintf('%spathway2GENE.tsv', 
                             output_dir), 
              sep = "\t", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}

pathway2name <- function(output_dir){
  
  pathway2NAME_list <- AnnotationDbi::as.list( KEGG.db::KEGGPATHID2NAME)
  pathways <- names(pathway2NAME_list)
  names(pathways) <- pathways
  pathway2NAME<- map_df(.x = pathways, 
                        .f = function(pathway=.x) {
                          data.frame(pathway=paste0("path:map",pathway), name=pathway2NAME_list[[pathway]])
                          })
  write.table(x = pathway2NAME, 
              file = sprintf('%s/pathway2NAME.tsv', 
                             output_dir), 
              sep = "\t", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}
plot_pcoa <- function(vsd, col_data, axis1 = "PC1", axis2 = "PC2", color_factor,
                      shape_factor = NULL, plot_type = "3D", syms2use = NULL, sample_column='sample'){
  
  
  # EXAMPLES
  
  # # Make 2D plot without a shape factor
  # plot_pcoa(vsd, col_data, axis1 = "PC1", axis2 = "PC2", color_factor= 'Treatment', plot_type = "2D")
  # 
  # # Make 2D plot with a shape factor
  # plot_pcoa(vsd, col_data, axis1 = "PC1", axis2 = "PC2", color_factor= 'Treatment', shape_factor = 'Treatment',
  #           plot_type = "2D")
  # 
  # # Make 3D plot without a shape factor
  # plot_pcoa(vsd, col_data, color_factor= 'Treatment', plot_type = "3D")
  # 
  # # Make 3D plot with a shape factor
  # # avalable symbol names can be found here https://plot.ly/r/reference/#scatter-marker-symbol
  # my_symbols <- c(Body_Female='diamond',Body_Immat='triangle-down',Body_Male='square',Body_Young_F='circle',
  #                 Eye_Female= "cross", Eye_Immat= "pentagon-dot", Eye_Male= "octagon-open", Eye_Young_F= "star",
  #                 PL_NFD = "star-diamond", PL_FD = "asterisk")
  # plot_pcoa(vsd, col_data, color_factor= 'Treatment', 
  #           shape_factor = 'Treatment', plot_type = "3D", syms2use = my_symbols)
  
  vsd <- t(assay(vsd)) %>% as.data.frame()
  samples <- col_data[,sample_column] %>% as.character()
  vsd <- vsd[samples, ]
  vsd_pca <- prcomp(vsd)
  #vsd_pca <- prcomp(t(log(counts_data)),center=F, scale=F)
  vsd_pca_percentVar <- vsd_pca$sdev^2/sum(vsd_pca$sdev^2)
  names(vsd_pca_percentVar) <- colnames(vsd_pca$x)
  
  if(plot_type == "3D"){
    
    if(!is.null(shape_factor)){
      pca2plot3d <- data.frame(vsd_pca$x[,c("PC1","PC2","PC3")],
                               Color=col_data[,color_factor],
                               Shape=col_data[,shape_factor],
                               Sample=str_replace(string =  col_data[,sample_column],
                                                  pattern = "^X",
                                                  replacement = ""))
      # pca2plot3d$Shape <- 
      #   paste0(substr(pca2plot3d$Shape,3,4),
      #          substr(pca2plot3d$Shape,1,1)) 
      
      p <-plot_ly(pca2plot3d) %>% 
        add_markers(x = ~PC1, y = ~PC2, z = ~PC3, 
                    colors = rainbow(4),     #c("red","blue"),
                    text=~Sample, 
                    color = ~Color, 
                    symbol=~Shape,
                    marker = list(symbol = syms2use),
                    symbols = syms2use)     %>%
        layout(scene = list(xaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC1",
                                                         vsd_pca_percentVar[axis1]*100)),
                            yaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC2",
                                                         vsd_pca_percentVar[axis1]*100)),
                            zaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC3",
                                                         vsd_pca_percentVar[axis1]*100))))
    }else{
      
      pca2plot3d <- data.frame(vsd_pca$x[,c("PC1","PC2","PC3")],
                               Color=col_data[,color_factor],
                               Sample=str_replace(string = col_data[,sample_column],
                                                  pattern = "^X",
                                                  replacement = ""))
      
      p <-plot_ly(pca2plot3d) %>% 
        add_markers(x = ~PC1, y = ~PC2, z = ~PC3, 
                    colors = rainbow(4),     #c("red","blue"),
                    text=~Sample, 
                    color = ~Color)     %>%
        layout(scene = list(xaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC1",
                                                         vsd_pca_percentVar[axis1]*100)),
                            yaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC2",
                                                         vsd_pca_percentVar[axis1]*100)),
                            zaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC3",
                                                         vsd_pca_percentVar[axis1]*100))))
      
      
    }
    
    
  }else{
    
    if(!is.null(shape_factor)){
      
      pca2plot <- data.frame(vsd_pca$x[,c(axis1,axis2)],
                             Color=col_data[,color_factor],
                             Shape=col_data[,shape_factor],
                             Sample=str_replace(string = col_data[,sample_column],
                                                pattern = "^X",
                                                replacement = ""))
      
      p <- ggplot(pca2plot,aes(x=get(axis1),
                               y=get(axis2),
                               col=Color,
                               shape=Shape)) +
        geom_point(size=5) +
        # geom_label_repel(aes(label=Sample)) +
        # geom_path(aes(color=Type,
        #               group=Type)) +
        xlab(sprintf("%s: %1.2f%% of variance", axis1, vsd_pca_percentVar[axis1]*100)) + 
        ylab(sprintf("%s: %1.2f%% of variance", axis2, vsd_pca_percentVar[axis2]*100)) + 
        labs(color=color_factor, shape=shape_factor) +
        coord_fixed()
      
    }else{
      
      pca2plot <- data.frame(vsd_pca$x[,c(axis1,axis2)],
                             Color=col_data[,color_factor],
                             Sample=str_replace(string = col_data[,sample_column],
                                                pattern = "^X",
                                                replacement = ""))
      p <- ggplot(pca2plot,aes(x=get(axis1),
                               y=get(axis2),
                               col=Color)) +
        geom_point(size=5) +
        # geom_label_repel(aes(label=Sample)) +
        # geom_path(aes(color=Type,
        #               group=Type)) +
        xlab(sprintf("%s: %1.2f%% of variance", axis1, vsd_pca_percentVar[axis1]*100)) + 
        ylab(sprintf("%s: %1.2f%% of variance", axis2, vsd_pca_percentVar[axis2]*100)) + 
        labs(color=color_factor, shape=shape_factor) +
        coord_fixed()
      
    }
    
  }
  
  return(p)
  
}

make_pcoa <- function(ps.rarefied,
                          method = "PCoA", 
                          distance= "bray",
                          metadata, ellipse=FALSE, 
                          color_var='MembraneType',
                          shape_var=NULL,
                          palette = c("red","green","blue","yellow","pink")){
  library(ggpubr)
  library(tools)
  library(phyloseq)
  
  # Add metadata to the phyloseq object
  #note a column in your metadata must be called samples else you'll get an error
  #here I add sample column which is the rownames of the passed metadata
  metadata <- data.frame(sample=rownames(metadata),metadata)
  sample_data(ps.rarefied) <- metadata
  common.ids <- intersect(rownames(metadata), colnames(otu_table(ps.rarefied)))
  metadata <- metadata[common.ids,]
  
  # make pcoa plot
  pcoa <- ordinate(ps.rarefied,method = method, distance = distance)
  
  
  pcoa_gg <- plot_ordination(physeq = ps.rarefied,
                             ordination = pcoa,
                             type = "samples",
                             axes = c(1,2), 
                             color = color_var,
                             shape = shape_var)
  
  if(!is.null(shape_var)){
    
    pcoa_gg <- ggscatter(data = pcoa_gg$data, x = 'Axis.1', y = 'Axis.2',
                         color = color_var, shape = shape_var, size = 6, 
                         ellipse = ellipse, palette = palette)
  }else{
    pcoa_gg <- ggscatter(data = pcoa_gg$data, x = 'Axis.1', y = 'Axis.2',
                         color = color_var,size = 6, ellipse = ellipse, palette = palette)  
  }
  
  
  #plot the axes ticks inwards 
  pcoa_gg <- pcoa_gg + theme(axis.ticks.length=unit(-0.15, "cm"),
                             axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
                             axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")))
  
  
  if(method == "PCoA"){
    Cumul_eig <- ncol(pcoa$values) - 1
    pcoa_gg <- pcoa_gg +
      labs(x= sprintf('PC1 [%s%%]', round(pcoa$values[1, Cumul_eig] * 100, digits = 2)),
           y= sprintf('PC2 [%s%%]', round((pcoa$values[2, Cumul_eig] - pcoa$values[1, Cumul_eig]) * 100,
                                          digits = 2)))
  }
  
  pcoa_gg <- ggpar(pcoa_gg,font.x = c(18,'bold.italic','black'), font.tickslab =  c(16,'bold','black'),
                   font.y = c(18,'bold.italic','black'),
                   legend.title = toTitleCase(gsub(pattern = '\\.', replacement = ' ', x = groups)), 
                   legend = 'right', font.legend = c(14, "bold", "black"))
  legend_titles <- toTitleCase(gsub(pattern = '\\.', replacement = ' ', x = groups))
  
  if(length(groups) == 2){
    pcoa_gg <- pcoa_gg + labs(color=legend_titles[1], shape=legend_titles[2])
  }
  
  return(pcoa_gg)
  
}



################################### Continue from here
make_3d_pcoa <- function(ps.rarefied,
                         method = "PCoA", 
                         distance= "bray",
                         metadata,
                         show_sample_name= FALSE,
                         color_var,
                         shape_var=NULL, 
                         size_var=NULL,
                         syms2use=NULL,
                         sizes2use=NULL,
                         sample_column="Sample_name",
                         palette=rainbow(10), ...){
  
  # ... other options that can be passed to plotly::add_markers function
  library(plotly)
  library(tools)
  library(phyloseq)
  
  # Add metadata to the phyloseq object
  #note a column in your metadata must be called samples else you'll get an error
  #here I add sample column which is the rownames of the passed metadata
  metadata <- data.frame(sample=rownames(metadata),metadata)
  sample_data(ps.rarefied) <- metadata
  
  common.ids <- intersect(rownames(metadata), colnames(otu_table(ps.rarefied)))
  if(is_empty(common.ids)){
    common.ids <- intersect(rownames(metadata), rownames(otu_table(ps.rarefied)))
  }
  metadata <- metadata[common.ids,]
  
  # make pcoa plot
  pcoa <- ordinate(ps.rarefied,method = method, distance = distance)
  
  if(method == "PCoA"){
    Cumul_eig <- ncol(pcoa$values) - 1
      x_lab <-  sprintf('PC1 [%s%%]',
                        round(pcoa$values[1, Cumul_eig] * 100, digits = 2))
      y_lab <- sprintf('PC2 [%s%%]', 
                       round((pcoa$values[2, Cumul_eig] - pcoa$values[1, Cumul_eig]) * 100,
                                          digits = 2))
      z_lab <- sprintf('PC3 [%s%%]', 
                       round((pcoa$values[3, Cumul_eig] - pcoa$values[2, Cumul_eig]) * 100,
                             digits = 2))
  }
  
  vectors_df <- as.data.frame(pcoa$vectors)

    
    if(!is.null(shape_factor)){
      pca2plot3d <- data.frame(vectors_df[,c("Axis.1","Axis.2","Axis.3")],
                               Color=metadata[,color_factor],
                               Shape=metadata[,shape_factor],
                               Sample=str_replace(string =  metadata[,sample_column],
                                                  pattern = "^X",
                                                  replacement = ""))
      
      
      
      p <- plot_ly(pca2plot3d) %>% 
        add_markers(x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, 
                    colors = palette,
                    text=~Sample, 
                    color = ~Color, 
                    symbol=~Shape, ...)     %>%
        layout(scene = list(xaxis = list(title = x_lab),
                            yaxis = list(title = y_lab),
                            zaxis = list(title = z_lab)))
    }else{
      
      pca2plot3d <- data.frame(vectors_df[,c("Axis.1","Axis.2","Axis.3")],
                               Color=metadata[,color_factor],
                               Sample=str_replace(string = metadata[,sample_column],
                                                  pattern = "^X",
                                                  replacement = ""))
      
      p <-plot_ly(pca2plot3d) %>% 
        add_markers(x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, 
                    colors = palette,
                    text=~Sample, 
                    color = ~Color, ... = )     %>%
        layout(scene = list(xaxis = list(title = x_lab),
                            yaxis = list(title = y_lab),
                            zaxis = list(title = z_lab)))
      
      
    }
    
    
  
  return(p)
  
}


# Plot 3D PCA or pcoa
make3D_plot <- function(mat, metadata, method="PCoA",  distance="bray", axis1 = "Axis.1", axis2 = "Axis.2", axis3 = "Axis.3", color_factor,
                        shape_factor = NULL, size_factor = NULL, sample_column='Sample.ID',
                        syms2use=NULL, colors2use=rainbow(4), sizes2use=NULL){
  # mat - matrix with samples as rows and column as features or phyloseq object for PCoA
  # Example
  # my_symbols <- c(NVNB='square',CAR='circle', NV='triangle-up')
  # 
  # # The size of the points mxust be specied as a ranged and not as explicit numbers
  # sizes <- c(20,60)
  # For PCA
  # make3D_plot(mat = t(clr_table), metadata = no_control_metadata,
  #            axis1 = "PC1", axis2 = "PC2", axis3 = "PC3",
  #            color_factor=color_factor, shape_factor = shape_factor,
  #            size_factor = size_factor, sample_column='Sample.ID',
  #            syms2use=my_symbols, 
  #            colors2use=annotation_colors[[color_factor]], 
  #            sizes2use=sizes)
  #  For PCoA
  # make3D_plot(mat = ps.rarefied, metadata = no_control_metadata,
  #            axis1 = "Axis.1", axis2 = "Axis.2", axis3 = "Axis.3",
  #            color_factor=color_factor, shape_factor = shape_factor,
  #            size_factor = size_factor, sample_column='Sample.ID',
  #            syms2use=my_symbols, 
  #            colors2use=annotation_colors[[color_factor]], 
  #            sizes2use=sizes)
  
  library(plotly)
  library(tools)
  library(phyloseq)
  
  # Define labels for the 3 axis based on the input
  XLAB <- sub("Axis\\.", "PC", axis1)
  YLAB <- sub("Axis\\.", "PC", axis2)
  ZLAB <- sub("Axis\\.", "PC", axis3)
  
  
  if(method == "PCoA"){
    
    # Add metadata to the phyloseq object
    #note a column in your metadata must be called samples else you'll get an error
    #here I add sample column which is the rownames of the passed metadata
    metadata <- data.frame(sample=rownames(metadata),metadata)
    sample_data(mat) <- metadata
    
    common.ids <- intersect(rownames(metadata), colnames(otu_table(mat)))
    if(is_empty(common.ids)){
      common.ids <- intersect(rownames(metadata), rownames(otu_table(mat)))
    }
    metadata <- metadata[common.ids,]
    
    # make pcoa plot
    pcoa <- ordinate(mat,method = method, distance = distance)
    
    vectors_df <- as.data.frame(pcoa$vectors)
    df <- vectors_df[,c(axis1, axis2, axis3)]
    
    
    colnames(df) <- c("PC1", "PC2", "PC3")
    Cumul_eig <- ncol(pcoa$values) - 1
    x_lab <-  sprintf('%s [%s%%]', XLAB,
                      round(pcoa$values[sub("Axis\\.", "", axis1), Cumul_eig] * 100, digits = 2))
    y_lab <- sprintf('%s [%s%%]', YLAB,
                     round((pcoa$values[sub("Axis\\.", "", axis2), Cumul_eig] - pcoa$values[1, Cumul_eig]) * 100,
                           digits = 2))
    z_lab <- sprintf('%s [%s%%]', ZLAB,
                     round((pcoa$values[sub("Axis\\.", "", axis3), Cumul_eig] - pcoa$values[2, Cumul_eig]) * 100,
                           digits = 2))
  }else{
    
    samples <- rownames(metadata)
    mat <- mat[samples,]
    mat_pca <- prcomp(mat)
    mat_pca_percentVar <- mat_pca$sdev^2/sum(mat_pca$sdev^2)
    names(mat_pca_percentVar) <- colnames(mat_pca$x)
    
    df <- mat_pca$x[,c(axis1, axis2, axis3)]
    
    x_lab <-  sprintf("%s: %1.2f%% of variance", XLAB, mat_pca_percentVar[axis1]*100)
    y_lab <- sprintf("%s: %1.2f%% of variance", YLAB, mat_pca_percentVar[axis2]*100)
    z_lab <- sprintf("%s: %1.2f%% of variance", ZLAB, mat_pca_percentVar[axis3]*100)
    
  }
  
  # Color, shape and size factors are provided
  if(!is.null(shape_factor) && !is.null(size_factor)){
    
    pca3plot <- data.frame(df,
                           Color=metadata[,color_factor],
                           Shape=metadata[,shape_factor],
                           Size=as.numeric(metadata[,size_factor]),
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z = ~PC3, 
                 mode = "markers",
                 color = ~Color,
                 colors = colors2use,
                 text= ~Sample, 
                 symbol= ~Shape,
                 symbols = syms2use,
                 marker = list(sizemode = 'diameter'),
                 size= ~Size,
                 sizes= sizes2use) 
    
    # Color and shape factors are provided
  } else if(!is.null(shape_factor)){
    
    pca3plot <- data.frame(df,
                           Color=metadata[,color_factor],
                           Shape=metadata[,shape_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z = ~PC3, 
                 mode = "markers",
                 color = ~Color,
                 colors = colors2use,
                 text= ~Sample, 
                 symbol= ~Shape,
                 symbols = syms2use)
    
    # Only a Color factor is provided
  }else{
    
    pca3plot <- data.frame(df,
                           Color=metadata[,color_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z=~PC3,
                 mode = "markers",
                 colors = colors2use,
                 text=~Sample, 
                 color = ~Color)
    
  }
  
  
  
  p <-  p  %>%
    layout(scene = list(xaxis = list(title = x_lab),
                        yaxis = list(title = y_lab),
                        zaxis = list(title = z_lab) )
           
    )
  
  return(p)
  
}







unique_factor_combination <- function(data, factor_column = 'Treatment'){
  
  # Get the unique levels of the factor column
  levels <- unique(as.character(data[,factor_column]))
  
  # Get all pairwise combinations of the factor
  combinations <- utils::combn(levels,2)
  return(combinations)
}

make_significant_heatmap <- function(res_df, vsd, col_data, group1, group2 = NULL, group3 = NULL,top = 50000,
                                     plot_dir, analysis_name, hmcol = c("red", "black", "green"),
                                     factor_order = NULL, factor_name = NULL,width = 15, height = 10,
                                     annotation_colors = list(treatment=c(Control='blue',B12_enriched='brown'),
                                                              WaterType=c(PW="blue", TWW="purple"), 
                                                              SoilType=c(Clay='chocolate4',Loam='chocolate3',Loamy_sand='chocolate1')),
                                     add_row_annotation=TRUE, sample_order=NULL,
                                     breaksList = NULL, save_plot=TRUE, ...) {
  
  
  # subset the table to contain the same sample names
  common_ids <- intersect(colnames(vsd), rownames(col_data))
  #res_df <- res_df[,common_ids]
  vsd <- vsd[,common_ids]
  col_data <- col_data[common_ids,]
  if(!is.null(factor_name) & !is.null(factor_order)){
    col_data[,factor_name] = factor(x = col_data[,factor_name], levels = factor_order)
  }
  
  # Get genes to plot in descending padj order. use top N genes for plotting
  
  padj_cols <- str_detect(string = names(res_df),
                          pattern = "padj")
  
  genes2plot <-
    res_df %>%
    # Add column with minimum padj:
    cbind(min_padj=apply(X = res_df[,padj_cols, drop=FALSE] ,
                         MARGIN=1,
                         FUN = function(x) ifelse(test = all(is.na(x)),
                                                  yes = NA,
                                                  no = min(x,na.rm = T)))) %>%
    as_tibble %>%
    # Select only rows that passed any padj
    # Select only rows that passed on day 1:
    filter(pass_any=="1") %>% 
    # filter(min_padj<=DESEQ_PADJ_CUTOFF) %>% 
    # Select top 2000 genes by descending padj:
    top_n(n = top, wt = dplyr::desc(min_padj)) %>% 
    dplyr::select(gene) %>%
    unlist
  
  
  # Get to locations of genes to plot:
  t1 <- rownames(vsd) %in% genes2plot
  
  # Extract data from VSD
  expr_data_list <- vsd[t1,]
  
  # Get numeric data:
  mat2plot <- assay(expr_data_list)
  
  # Get column data:
  col_data_expr <- col_data
  
  # Scale by row
  mat2plot <- mat2plot %>% t %>% scale %>% t
  
  # remove gene rows for which no SD was coming
  mat2plot <- mat2plot[apply(mat2plot, MARGIN = 1, FUN = function(x) sd(x) != 0),]
  
  # re-arrange the order of samples for ploting on the heatmap
  if(!is.null(sample_order)){
    mat2plot <- mat2plot[,sample_order]
  }
  
  # Create col_annotation and
  # Sort col_annotation by factors in col_data:
  # add group2 name if provided
  if(is.null(group2)){
    col_annotation <- data.frame(group1=col_data_expr[,c(group1)],
                                 row.names = rownames(col_data_expr)) 
    colnames(col_annotation)[1] <- group1 
  }else if(!is.null(group2) && !is.null(group3)){
    
    col_annotation <- data.frame(group1=col_data_expr[,c(group1)],
                                 group2=col_data_expr[,c(group2)],
                                 group3=col_data_expr[,c(group3)],
                                 row.names = rownames(col_data_expr)) 
    colnames(col_annotation)[1] <- group1 
    colnames(col_annotation)[2] <- group2 
    colnames(col_annotation)[3] <- group3 
  }else{
    col_annotation <- data.frame(group1=col_data_expr[,c(group1)],
                                 group2=col_data_expr[,c(group2)],
                                 row.names = rownames(col_data_expr)) 
    colnames(col_annotation)[1] <- group1 
    colnames(col_annotation)[2] <- group2 
  }
  
  col_annotation$sample <- NULL
  
  
  # manually create color range
  hmcol = c("red", "black", "green")
  # expand the color range
  hmcol = colorRampPalette(hmcol)(255)
  if(save_plot){
  png(filename = sprintf("%s/%s.heatmap.png",plot_dir,analysis_name),
      width = width,
      height = height,
      units = "in",
      res=300)
  }
  # Get row annotation
  row_annotation <- res_df[, str_detect(names(res_df),"pass") & !names(res_df)=="pass_any",drop=F]
  colnames(row_annotation) <- colnames(row_annotation) %>% str_replace("pass.","")
  
  rownames(row_annotation) <- res_df$gene
  row_annotation <- row_annotation[row.names(mat2plot),,drop=F]
  row_annotation[row_annotation==""] <- NA
  
  pheatmap_data <-
    pheatmap(mat2plot,
             clustering_method = "complete",
             clustering_distance_rows = "correlation",
             color = rev(hmcol),
             cluster_cols = FALSE,
             annotation_row =if(add_row_annotation) row_annotation else NA,
             annotation_col = col_annotation,
             annotation_colors = annotation_colors,
             show_rownames = if(nrow(mat2plot) <= 30) TRUE else FALSE,
             width = width, height = height, ...)
  
  if(save_plot){
  dev.off()
  }
  final <- list(pheatmap_data, row_annotation)
  names(final) <- c('pheatmap_data','row_annotation')
  return(final)
  
}


# by mapping file
subset_resdf <- function(resdf, original_mapping_file, subset_mapping_file, 
                         treatment_column="Treatment") {
  
  # EXAMPLE
  #subset_resdf(resdf = res_df, original_mapping_file=col_data, subset_mapping_file=mapAB, treatment_column="Treatment")
  
  original_treatments <- unique(as.character(original_mapping_file[,treatment_column]))  
  subset_treatments <- unique(as.character(subset_mapping_file[,treatment_column]))
  
  treatments_to_exclude <- original_treatments[!original_treatments %in% subset_treatments]   
  
  exclusion_pattern <- paste(treatments_to_exclude, collapse = '|')
  
  treatment_cols <- colnames(resdf)[grepl(pattern = paste(subset_treatments, collapse = '|'), 
                                          x = colnames(resdf)) & 
                                      !grepl(pattern = exclusion_pattern, x = colnames(resdf))]
  
  # Get the annotation column names - there are 10 annotation columns from the end of res_df 
  start <- ncol(resdf)-10
  end <- ncol(resdf)
  annotation_cols <- colnames(resdf)[start:end]
  
  columns_to_select <- c("gene",as.character(subset_mapping_file$SampleID),treatment_cols,annotation_cols ) 
  
  return(columns_to_select)
  
}

# by comparison names
subset_resdf_by_comparison <- function(resdf, col_data,comparison_names, condition){ 
  columns_to_exclude <-  comparison_names[condition]
  comparisons_to_keep <- setdiff(comparison_names, columns_to_exclude)
  pattern_to_keep <- paste0(comparisons_to_keep, collapse = "|")
  treatment_cols <- colnames(resdf)[grepl(pattern = pattern_to_keep, x = colnames(resdf))]
  # Get the annotation column names - there are 10 annotation columns from the end of res_df 
  start <- ncol(resdf)-10
  end <- ncol(resdf)
  annotation_cols <- colnames(resdf)[start:end]
  
  columns_to_select <- c("gene",as.character(col_data$SampleID),treatment_cols,annotation_cols)
  return(columns_to_select)
}


# functions for DESEq analysis
get_results <- function(ds,y,treatment,independentVariable,alpha=0.1,LINEAR_FC_CUTOFF=1.3,PADJ_CUTOFF=0.05){
  gene_names <- rownames(ds)
  results_df <- data.frame()
  if(y != treatment ){
    results_df <- results(ds, contrast=c(independentVariable, treatment,y), 
                          alpha=alpha)%>% 
      as.data.frame %>% 
      as.tibble %>% 
      mutate(linearFC = ifelse(is.na(log2FoldChange),
                               yes = NA,
                               no = ifelse(log2FoldChange>0,
                                           yes = 2^log2FoldChange,
                                           no = -1/(2^log2FoldChange))%>% 
                                 signif(digits = 3))) %>% 
      mutate(pvalue = ifelse(test = is.na(pvalue),
                             yes = NA,
                             no = signif(pvalue,digits = 3))) %>% 
      mutate(padj = ifelse(test = is.na(padj),
                           yes = NA,
                           no = signif(padj,digits = 3))) %>% 
      mutate(pass = ifelse(test = abs(as.numeric(linearFC)) >= LINEAR_FC_CUTOFF & 
                             padj <= PADJ_CUTOFF & 
                             !is.na(padj),
                           yes = ifelse(test = as.numeric(linearFC)>0,
                                        yes="up",
                                        no="down"),
                           no = "")) %>% 
      as.data.frame
    rownames(results_df) <- gene_names
  }
  
  return(results_df)
}

deseq_results<- function(ds,treatment, treatments, independentVariable){
  #ds is the intial output from DESEQ analysis
  # treatment is a level/group in the independentVariable
  # treatments is a character vector of all the levels in independentVariable
  # independentVariable is a string representin the independent variable that was analysed
  #USAGE
  # get_deseq_results(ds=ds,treatment="NS+NF+NW+NC",treatments=treatments, independentVariable='treatment')
  
  # .x is a sub-treatment/level in the treatments variable
  map(.x = treatments, .f = ~get_results(ds,.x,treatment,independentVariable))
  
}

significant_deseq_results <- function(res,alpha=0.05){
  res <- res[order(res$padj, na.last=NA), ]
  res_sig <- res[(res$padj < alpha), ]
  return(res_sig)
}

deseq2dataframe <- function(res,treatment1,treatment2){
  
  current_result <- res[[treatment1]][[treatment2]]
  if(!is.null(current_result)){
    prefix <- paste0(treatment1, ' vs ')
    prefix <- rep(prefix, nrow(current_result))
    result <- data.frame(comparison=paste(prefix,rep(treatment2, nrow(current_result))), 
                         as.data.frame(current_result, stringsAsFactors=FALSE), stringsAsFactors=FALSE)
    
    result <- data.frame(Taxon=rownames(result), result, stringsAsFactors=FALSE)
    if(nrow(result) > 0){
      rownames(result) <- 1:nrow(result)
    }
    return(result)
  }
  
}

concat_deseq_results <- function(treatments,result_list){
  #concat deseq's list of results to one dataframe
  map_df(.x = treatments,.f = function(x){
    map_df(.x = treatments, .f = ~deseq2dataframe(res = result_list, treatment1 = x, treatment2 = .x))
  })
  
}

round_numeric_columns <- function(data,numColStart){
  if(nrow(data) > 0){
    data[,numColStart:ncol(data)] <- apply(X = data[,numColStart:ncol(data)],
                                           MARGIN = 2, FUN = function(x) { round(x,3) })
  }
  
  return(data)
}

deseq_result_dataframe <- function(ds,treatments,independentVariable='treatment'){
  # Function to generate dataframe of deseq results both significant and non significant results
  #ds is the intial output from DESEQ analysis
  # treatments is a character vector of all the levels in independentVariable
  # independentVariable is a string representin the independent variable that was analysed
  #USAGE
  #deseq_result_dataframe(ds=ds,treatments=treatments,independentVariable='treatment')
  
  
  # Get deseq results for all treatment combinations
  res <- map(.x = treatments, 
             .f = ~deseq_results(ds=ds,treatment = .x, treatments = treatments, independentVariable =independentVariable))
  
  # name the results according to treatment names
  names(res) <- treatments
  for (treatment in treatments) {
    
    names(res[[treatment]]) <- treatments
  }
  
  
  # Get the significant results 
  res_sig <- map(.x = res, .f = function(x){
    
    map(.x = x, .f = ~significant_deseq_results(.x))
  })
  
  # Concat deseq's list of results to one dataframe
  significant_results <- concat_deseq_results(treatments,res_sig)
  all_results <- concat_deseq_results(treatments,res)  
  
  # Round numeric columns to 3 decimal places
  significant_results <- round_numeric_columns(significant_results,numColStart=3)
  all_results <- round_numeric_columns(all_results,numColStart=3)
  
  #Arrange by taxon name
  significant_results <- significant_results %>% arrange(Taxon)
  all_results <- all_results %>% arrange(Taxon)
  
  final <- list(all_results, significant_results)
  names(final) <- c('all_results','significant_results')
  
  return(final)
}

# GOseq gene enrichment analysis
get_genes_binary_vector <- function(result_df){
  genes <- as.integer(result_df$pass !="")
  names(genes) <- rownames(result_df)
  return(genes)
} 


subset_prediction_table <- function(prediction_table, genelist, sep=';') {
  # genelist can be a file path or a character vector of KO ids only
  # prediction_table is your predicted table of functions by Tax4Fun
  # sep is a characters string specifying the seperator to be used to ge the KO names
  
  #value - return a subsetted data frame and  list if ginelist is a file 
  # but returns only the subsetted dataframe is genelist is a vector of KO names 
  
  # USAGE 
  # 1. passing a file path with Ko ids and there definitions on seperate lines
  # nitrogen_genes <- subset_prediction_table(prediction_table = function_profile,
  #                                           genelist = '../PhD/nitrogen_genes',
  #                                           sep = ';')
  
  # 2. passing a vector of KO ids
  # 
  
  sep <- paste0(sep,'|')
  if( length(genelist) == 1 && file.exists(genelist)){
    #sep <- paste0(sep,'|')
    # open, read and close the genelist file
    con <- file(description = genelist, open = 'r')
    lines <- readLines(con)
    close(con)
    # Get a vector of KO ids
    matches <- regexpr(pattern = '^(K\\d+)',text=lines)
    ids <- regmatches(lines,matches)
    ids <- unique(ids)
    # KO_list <- vector(mode = 'list', length = length(seq(1,length(nitrogen),2)))
    # # assumes the ids are followed by their description on separate lines
    # ids <- lines[seq(1,length(lines),2)]
    # description <- lines[seq(2,length(lines),2)]
    # names(KO_list) <- ids
    
    #for (i in seq_along(ids)) {
    #  KO_list[[i]] <- description[i]
    #}
    
    pattern <- paste(ids, collapse = sep)
    #final <- list(genes_df,KO_list)
    #names(final) <- c('gene_df','KO_list')
    #return(final)
    
  }else if(genelist != '' && is.character(genelist) && length(genelist) > 0 && all(startsWith(genelist,'K'))){
    #sep <- paste0(sep,'|')
    pattern <- paste(genelist, collapse = sep)
    
  }else{
    message("invalid genelist: can only be  a file path or a character vector of KO ids")
  }
  
  genes_idices <- grep(pattern = pattern,
                       x = colnames(prediction_table))
  
  genes_df <- prediction_table[,genes_idices]
  return(genes_df)
  
}




subset_prediction_table <- function(prediction_table, genelist, sep=';', Tax4Fun=TRUE) {
  # genelist can be a file path or a character vector of KO ids only
  # prediction_table is your predicted table of functions by Tax4Fun
  # sep is a characters string specifying the seperator to be used to ge the KO names
  
  #value - return a subsetted data frame and  list if ginelist is a file 
  # but returns only the subsetted dataframe is genelist is a vector of KO names 
  
  # USAGE 
  # 1. passing a file path with Ko ids and there definitions on seperate lines
  # nitrogen_genes <- subset_prediction_table(prediction_table = function_profile,
  #                                           genelist = '../PhD/nitrogen_genes',
  #                                           sep = ';')
  
  # 2. passing a vector of KO ids
  # 
  
  sep <- paste0(sep,'|')
  if( length(genelist) == 1 && file.exists(genelist)){
    
    # open, read and close the genelist file
    con <- file(description = genelist, open = 'r')
    lines <- readLines(con)
    close(con)
    # Get a vector of KO ids
    matches <- regexpr(pattern = '^(K\\d+)',text=lines)
    ids <- regmatches(lines,matches)
    ids <- unique(ids)
    
    
    pattern <- paste(ids, collapse = sep)
    
    
  }else if(genelist != '' && is.character(genelist) && length(genelist) > 0 && all(startsWith(genelist,'K'))){
    #sep <- paste0(sep,'|')
    pattern <- paste(genelist, collapse = sep)
    
  }else{
    message("invalid genelist: can only be  a file path or a character vector of KO ids")
  }
  
  if(Tax4Fun){
    genes_idices <- grep(pattern = pattern,
                         x = colnames(prediction_table))
    
    genes_df <- prediction_table[,genes_idices]
    
  }else{
    genes_idices <- grep(pattern = pattern,
                         x = rownames(prediction_table)) 
    
    genes_df <- prediction_table[genes_idices,]
    
  }
  
  return(genes_df)
  
}

# Beta diversity
betadiversity <- function(ps.rarefied,
                          method = "PCoA", 
                          distance= "bray",
                          groups = "treatments",
                          metadata, ellipse=FALSE, 
                          color_var='MembraneType',
                          shape_var=NULL,
                          palette = c("red","green","blue","yellow","pink")){
  library(ggpubr)
  library(tools)
  library(phyloseq)
  
  # Add metadata to the phyloseq object
  #note a column in your metadata must be called samples else you'll get an error
  #here I add sample column which is the rownames of the passed metadata
  metadata <- data.frame(sample=rownames(metadata),metadata)
  sample_data(ps.rarefied) <- metadata
  
  # make pcoa plot
  pcoa <- ordinate(ps.rarefied,method = method, distance = distance)
  
  
  pcoa_gg <- plot_ordination(physeq = ps.rarefied,
                             ordination = pcoa,
                             type = "samples",
                             axes = c(1,2), 
                             color = color_var,
                             shape = shape_var)
  
  if(length(groups)== 2){
    
    pcoa_gg <- ggscatter(data = pcoa_gg$data, x = 'Axis.1', y = 'Axis.2',
                         color = color_var, shape = shape_var, size = 6, 
                         ellipse = ellipse, palette = palette)
  }else{
    pcoa_gg <- ggscatter(data = pcoa_gg$data, x = 'Axis.1', y = 'Axis.2',
                         color = color_var,size = 6, ellipse = ellipse, palette = palette)  
  }
  
  
  #plot the axes ticks inwards 
  pcoa_gg <- pcoa_gg + theme(axis.ticks.length=unit(-0.15, "cm"),
                             axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
                             axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")))
  
  
  if(method == "PCoA"){
    Cumul_eig <- ncol(pcoa$values) - 1
    pcoa_gg <- pcoa_gg +
      labs(x= sprintf('PC1 [%s%%]', round(pcoa$values[1, Cumul_eig] * 100, digits = 2)),
           y= sprintf('PC2 [%s%%]', round((pcoa$values[2, Cumul_eig] - pcoa$values[1, Cumul_eig]) * 100,
                                          digits = 2)))
  }
  
  pcoa_gg <- ggpar(pcoa_gg,font.x = c(18,'bold.italic','black'), font.tickslab =  c(16,'bold','black'),
                   font.y = c(18,'bold.italic','black'),
                   legend.title = toTitleCase(gsub(pattern = '\\.', replacement = ' ', x = groups)), 
                   legend = 'right', font.legend = c(14, "bold", "black"))
  legend_titles <- toTitleCase(gsub(pattern = '\\.', replacement = ' ', x = groups))
  
  if(length(groups) == 2){
    pcoa_gg <- pcoa_gg + labs(color=legend_titles[1], shape=legend_titles[2])
  }
  
  #find if there is a significant clustering between treatments using adonis test aka 
  #PERMANOVA
  #calculate and get your distance matrix - here we use bray curtis but check
  # the help page for more options
  bray_dist <- phyloseq::distance(ps.rarefied, method="bray")
  permanova<- vector(mode = 'list', length = length(groups))
  names(permanova) <- groups
  
  for (group in groups){
    permanova[[group]] <- adonis(bray_dist  ~ metadata[,group])
  }
  
  result <- list(pcoa,pcoa_gg,permanova,bray_dist)
  names(result) <- c("ordination","plot","permanova","bray_dist")
  
  return(result)
  
}


make_sample_distance_heatmap <- function(vsd_matrix, plot_title, colors){
  mat <- t(assay(vsd_matrix))
  # Checking sample to sample distances:
  sampleDists <- dist(mat)
  sampleDistMatrix <- as.matrix(sampleDists)
  # png(filename = sprintf("%s/sample_heatmap.png",plot_dir),
  #     width = 14, height = 8, res = 300, units = 'in')
  samples_heatmap  <- pheatmap(sampleDistMatrix,
                               clustering_distance_rows=sampleDists,
                               clustering_distance_cols=sampleDists,
                               col=colors, main = plot_title)
  
  return(samples_heatmap)
  #dev.off()
}


draw_heatmap <- function(t_OTU_table, color, col.order=NULL,
                         row_names_gp=gpar(fontface="bold",cex=0.6),
                         column_names_gp=gpar(fontface="bold"),
                         row_names_side = 'right', 
                         rect_gp = gpar(col="grey"), ...) {
  #t_OTU_table is the raw trasposed OTU table with the microbe names as column names
  #col.order is how you'll like to arrange the treatment names
  #USAGE
  #draw_heatmap(t_OTU_table = (treat_top_25*100), col.order = c("UnIrrigated_clay", "PW-P","TWW-P","PW+P", "TWW+P"))
  
  
  t_OTU_table <- t(t_OTU_table)
  #get the index of the ordered pathogen names
  id<- order(rownames(t_OTU_table))
  #get their names based on the idex
  id<- rownames(t_OTU_table)[id]
  #re-arrange the names alphabetically
  t_OTU_table <- t_OTU_table[id, ]
  if(!is.null(col.order)){
    t_OTU_table <- t_OTU_table[,col.order]
  }
  
  #f6 = colorRamp2(seq(0, max(t_OTU_table), length = 2), c("white", "red"))
  #make heatmap
  heatMap <- Heatmap(t_OTU_table,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     #column_title = 'Treatment',
                     row_names_side = row_names_side,
                     column_title_side = 'bottom',
                     #clustering_distance_columns = function(x) vegdist(x),
                     row_names_gp = row_names_gp,
                     column_names_gp = column_names_gp,
                     col = color,
                     rect_gp = rect_gp,
                     heatmap_legend_param = list(title= 'Relative Abundance(%)'), ...)
  
  
  #display the heatmap
  return(heatMap)
}



#function for comparing between 2 groups - it performs both t-test and wilcox test
tTest_analysis <- function(x,y,Data){
  #x is a string specifying the name of the independent variable
  #y is a string specifying the name of the dependent variable 
  #Data is a dataframe containing the variables to be analysed
  
  #if all the values are the same
  #if(var(Data[,y]) == 0 || abs(max(Data[,y]) - min(Data[,y])) < 0.02){
  if(var(Data[,y]) == 0){
    indicator.var <- NA; shapiro_pval<- NA; tTest_pval <- 1; tTest_statistic <- NA;
    wilcox_pval <- 1; wilcox_statistic <- NA;
	    # MEan and sD
    mean_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y])
    se_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y], FUNC = SE,FUNC_string = 'se')
    # Mean_rank
    rank_mean_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC_string = 'rank_mean')
    rank_se_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC = SE,FUNC_string = 'rank_se')
  }else if ((length(unique(Data[,y])) > (length(Data[,y])/20))){
    indicator.var <-var.test(Data[,y] ~ Data[,x])$p.value
    
    #check that you have unique values for y
    ifelse(test = length(unique(Data[,y])) > (length(Data[,y])/2),
           yes = shapiro_pval<- (shapiro.test(Data[,y]))$p.value, no = shapiro_pval<-0)
    tTest_res <- t.test(Data[,y] ~ Data[,x])
    tTest_pval <- tTest_res$p.value
    
    tTest_statistic <- tTest_res$statistic
    
    wilcox_res <- wilcox.test(Data[,y] ~ Data[,x])
    wilcox_pval<- wilcox_res$p.value
    
    wilcox_statistic<- wilcox_res$statistic
    # MEan and sD
    mean_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y])
    se_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y], FUNC = SE,FUNC_string = 'se')
    # Mean_rank
    rank_mean_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC_string = 'rank_mean')
    rank_se_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC = SE,FUNC_string = 'rank_se')
  }else{
    indicator.var <- NA; shapiro_pval<- NA; tTest_pval <- NA; tTest_statistic <- NA;
    wilcox_pval <- NA; wilcox_statistic <- NA;
    mean_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y])
    se_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y], FUNC = SE,FUNC_string = 'se')
    # Mean_rank
    rank_mean_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC_string = 'rank_mean')
    rank_se_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC = SE,FUNC_string = 'rank_se')
  }
  results <- data.frame(shapiro_pval = shapiro_pval ,
                        var_pval = indicator.var, 
                        tTest_statistic = tTest_statistic, 
                        tTest_pval = tTest_pval,
                        wilcox_statistic = wilcox_statistic,
                        wilcox_pval = wilcox_pval,
                        mean_df, se_df,
                        rank_mean_df,rank_se_df)
  
  return(results)
}

summarise_group <- function(x_vect, y_vect, FUNC=mean, FUNC_string="mean"){
  # A function to return a one row summary dataframe of
  # y-vect grouped by x_vect using a specified function FUNC
  
  df <- data.frame(x=x_vect, y=y_vect)

  res <- aggregate(y~x, data=df, FUN=FUNC)
  rownames(res) <- res$x
  res <- res[,-1, drop=F] %>% t
  colnames(res) <- paste0(colnames(res),sep=sprintf("_%s",FUNC_string))
  return(res)
}

#function for comparing more than 2 groups, it performs a one way anova and kruskal wallis test
anova_analysis <- function(x,y,Data){
  #x is a string specifying the name of the independent variable
  #y is a string specifying the name of the dependent variable 
  #Data is a dataframe containing the variables to be analysed
  if ((length(unique(Data[,y])) > (length(Data[,y])/20))){
    indicator.aov <-aov(Data[,y] ~ Data[,x])
    anova_res <- anova(indicator.aov)
    anova_pval <- anova_res$'Pr(>F)'[1]
    anova_fval <- anova_res$'F value'[1]
    group_DF <- anova_res$Df[1]
    residual_DF <- anova_res$Df[2]
    ifelse(test = length(unique(residuals(indicator.aov))) > (length(residuals(indicator.aov))/2),
           yes = shapiro_pval<- (shapiro.test(residuals(indicator.aov)))$p.value, no = shapiro_pval<-0)
    Levene_pval<- LeveneTest(lm(Data[,y] ~ Data[,x], data = Data))$'Pr(>F)'[1]
    kruskal_res <- kruskal.test(Data[,y] ~ Data[,x], data = Data)
    kruskal_statistic <- kruskal_res$statistic
    kruskal_pval<- kruskal_res$p.value
    # Means and SE
    mean_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y])
    se_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y], FUNC = SE,FUNC_string = 'se')
    # Mean_rank
    rank_mean_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC_string = 'rank_mean')
    rank_se_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC = SE,FUNC_string = 'rank_se')
  }else{
    anova_pval <- NA; anova_fval <- NA; group_DF <- NA; residual_DF <- NA; shapiro_pval <- NA;
    Levene_pval <- NA; kruskal_statistic <- NA; kruskal_pval<- NA;
    mean_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y])
    se_df <- summarise_group(x_vect = Data[,x], y_vect=Data[,y], FUNC = SE,FUNC_string = 'se')
    # Mean_rank
    rank_mean_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]),FUNC_string = 'rank_mean')
    rank_se_df <- summarise_group(x_vect = Data[,x], y_vect=base::rank(Data[,y]), FUNC = SE,FUNC_string = 'rank_se')
  }
  results <- data.frame(group_DF = group_DF ,
                        residual_DF = residual_DF, 
                        anova_fval = anova_fval ,
                        anova_pval = anova_pval,
                        shapiro_pval = shapiro_pval,
                        Levene_pval = Levene_pval , 
                        kruskal_statistic = kruskal_statistic, 
                        kruskal_pval = kruskal_pval, 
                        mean_df, se_df, 
                        rank_mean_df,rank_se_df)
  return(results)
}



multiple_posHoc_test <- function(data, group, colname="Taxon", post_hoc_method ="dunn.test"){
  #funtion to perform posthoc tests on multiple variable in a dataframe
  
  # data - is a dataframe that contains the group and the numeric variables to test
  # group - a a string that specifies the name of the grouping variable
  # colname - is a string that specifies the general name for the numeric variables
  # post_hoc_method - is a string that specifies the posthoc method either "dunn.test" or "tukeyHSD"
  library(FSA)
  numeric_variables <- colnames(data)[-(which(colnames(data) == group))]
  all_comparisons <- data.frame()
  
  
  if(post_hoc_method == "dunn.test"){
    for (variable in numeric_variables) {
      #print(variable)
      res<- dunnTest(data[,variable],data[,group])
      current_comparison <-cbind(rep(x = variable, times=nrow(res$res)),res$res)
      all_comparisons <- rbind(all_comparisons,current_comparison)
    }
  }else if(post_hoc_method == "tukeyHSD"){
    
    for (variable in numeric_variables) {
      aov.model <- aov(formula = reformulate(termlabels = group,response = variable), data = data)
      res <- TukeyHSD(aov.model)
      current_comparison <- cbind(rep(x = variable, times=nrow(res[[group]])),res[[group]])
      colnames(current_comparison)[1] <- colname
      current_comparison <- as.data.frame(current_comparison)
      current_comparison$comparison <- rownames(current_comparison)
      current_comparison <-  current_comparison[,c(colname,"comparison","diff","lwr","upr","p adj")]
      
      all_comparisons <- rbind(all_comparisons,current_comparison)
    }
    
  }else{
    stop("post_hoc_method must be either dunn.test or tukeyHSD")
  }
  rownames(all_comparisons) <- 1:nrow(all_comparisons)
  colnames(all_comparisons)[1] <- colname
  return(all_comparisons)
}


#function for pairwise comparison between groups
compare_groups <- function(x,y,Data, ...){
  library(ggpubr)
  #x is a string specifying the name of the independent variable
  #y is a string specifying the name of the dependent variable 
  #Data is a dataframe containing the variables to be analysed
  if((length(unique(Data[,y])) > (length(Data[,y])/10))){
    Formula <- paste(sprintf("%s ~ ", y), x, sep = "")
    Forma <- as.formula(Formula)
    comparison <- compare_means(Forma,data = Data, ...) 
    
  }else{
    comparison <- data.frame(.y. = c(NA,NA), group1=c(NA,NA), 
                             group2=c(NA,NA), p = c(NA,NA), p.adj=c(NA,NA),
                             p.format=c(NA,NA), p.signif=c(NA,NA), method=c(NA,NA))
    
  }
  
  return(comparison)
}

# function to make unfaceted and ungrouped line, bar or box plots
draw_plot <- function(x,y,xLab, yLab, Data, date,plot_type,indicator = NULL) {
  library(ggpubr)
  # x is a string specifying the independent variable
  # y is a string specifying the dependent variable
  # xLab is a string specifying the label for the x axis
  # yLab is a string specifying the label for the y axis
  # plot_type is a string that specifies the type of plot 'box' or 'line' or 'bar'
  #indicator is an optional string specifying the name to be added to the figure title instead of y
  
  #ggpubr fails to plot if ';' is part of the variable name
  #so I check for a potential problem then correct it here
  
  #check if a number precess a dash sign and replace it with the number and .
  if(grepl(pattern = "\\d+\\-", x = y)){
    old_y <- y
    y <- gsub(pattern ="(\\d+)\\-",replacement = "\\1\\.",x = y)
    # if(grepl(pattern = "\\d+;", x = y)){
    #   old_y <- y
    #   y <- gsub(pattern ="(\\d+);",replacement = "\\1\\.",x = y) 
    # }
    # if(grepl(pattern = "\\-|;|\\s+|\\(|\\)", x = y)){
    #   old_y <- y
    #   y <- gsub(pattern ="\\-|;|\\s+|\\(|\\)",replacement = "_",x = y)
    # }
    #rename the column with the acceptable name
    colnames(Data)[which(colnames(Data) == old_y)] <- y
  }else if(grepl(pattern = "\\-|;|\\s+|\\(|\\)", x = y)){
    old_y <- y
    y <- gsub(pattern ="\\-|;|\\s+|\\(|\\)",replacement = "_",x = y)
    #rename the column with the acceptable name
    colnames(Data)[which(colnames(Data) == old_y)] <- y
  }
  #debugging
  # print("here you go:-")
  # print(y)
  
  if (plot_type == 'box'){
    if (is.null(indicator) == TRUE ){
      plot.gg<- ggboxplot(data = Data,
                          x = x,
                          y = y,
                          xlab = xLab,
                          ylab = yLab,
                          title = sprintf("Diff. in %s %s",y,date),
                          ggtheme = theme_gray(),
                          palette = c('black','grey','white'),
                          add = "mean_se")
    }else{
      plot.gg<- ggboxplot(data = Data,
                          x = x,
                          y = y,
                          xlab = xLab,
                          ylab = yLab,
                          title = sprintf("Diff. in %s %s",indicator,date),
                          ggtheme = theme_gray(),
                          palette = c('black','grey','white'),
                          add = "mean_se")
    }
    plot.gg + 
      font("xlab", size = 16, face = "bold",color = "black")+
      font("ylab", size = 16, face = "bold",color = "black")+
      font("xy.text", size = 14, face = "bold",color = "black")
    
    
  } else if (plot_type == 'line'){
    if (is.null(indicator) == TRUE ){
      plot.gg<- ggline(data = Data,
                       x = x,
                       y = y,
                       xlab = xLab,
                       ylab = yLab,
                       point.size = 2,
                       size = 1.05,
                       title = sprintf("change in %s %s",y,date),
                       add = "mean_se")
    }else{
      plot.gg<- ggline(data = Data,
                       x = x,
                       y = y,
                       xlab = xLab,
                       ylab = yLab,
                       point.size = 2,
                       size = 1.05,
                       title = sprintf("change in %s %s",indicator,date),
                       add = "mean_se")  
    }
    plot.gg + 
      font("xlab", size = 16, face = "bold",color = "black")+
      font("ylab", size = 16, face = "bold",color = "black")+
      font("xy.text", size = 14, face = "bold",color = "black")
    
  } else if (plot_type == 'bar'){
    if (is.null(indicator) == TRUE ){
      plot.gg <- ggbarplot(data = Data,
                           x = x,
                           y = y,
                           xlab = xLab,
                           ylab = yLab,
                           title = sprintf("Diff. in %s %s",y,date),
                           fill = "grey",
                           ggtheme = theme_gray(),
                           palette = c('black','grey','white'),
                           add = "mean_se")
    }else{
      plot.gg <- ggbarplot(data = Data,
                           x = x,
                           y = y,
                           xlab = xLab,
                           ylab = yLab,
                           title = sprintf("Diff. in %s %s",indicator,date),
                           fill = "grey",
                           ggtheme = theme_gray(),
                           palette = c('black','grey','white'),
                           add = "mean_se")
    }
    plot.gg + 
      font("xlab", size = 16, face = "bold",color = "black")+
      font("ylab", size = 16, face = "bold",color = "black")+
      font("xy.text", size = 14, face = "bold",color = "black")
    
  } else {
    print("error you need to specify a valid plot_type c('box', 'line','bar')")
  }
  
}


#functions to perform differntial taxon/OTU abundances testing by running one-way ANOVA, 
#t-tests and their non-parameteric alternatives
non_parametric_diff_abundance <- function(OTU_table,MAtrix,season,group, 
                                          grouping_info,pAjust = 'fdr', 
                                          convert_percent=TRUE) {
  #a function to perform differential abundance testing of OTUs or taxons using 
  #kruskal Wallis or Wilcoxon test
  
  #otu_table - an otu matrix with colnames as OTUs and rownames as sample names 
  #matrix - a string for labeling in our cause c('soil', 'water','fruit_wash', 'fruit_puree')
  #1.OTU_table - normalized otu table matrix after running the  'filter_OTU_table' function
  #2 MAtrix and season - string for labeling the plot title
  #3 grouping_info - a dataframe of independent variables to be analysed derived at running 'filter_OTU_table' function
  #4.group - a string specifiying the independent variable within the mapping file to be analyzed
  #5.pAjust - is the p-value adjustment method either of c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  #   "fdr", "none")
  
  if (convert_percent){
    #convert relative abundace dat to percentages
    OTU_table <- OTU_table*100
  }
  
  
  #set the variables
  #get the variable name or the OTU names
  variables <- colnames(OTU_table)
  
  #get the dataframe of the independent variable under consideration
  independent_variable.df <- subset(x = grouping_info, select=group)
  
  #get the number of factor levels in the independent variable
  group_factor_length <- length(levels(independent_variable.df[,group]))
  
  #get only the sample present in both the OTU table and the metadata
  common.ids <- intersect(rownames(independent_variable.df), rownames(OTU_table))
  OTU_table <- OTU_table[common.ids,,drop=FALSE]
  independent_variable.df <- independent_variable.df[common.ids,,drop=FALSE]
  #column bind the independent variable to teh otu table
  OTU_table.df <- cbind(independent_variable.df,OTU_table)
  
  #initialize the results containing dataframes
  result.df <- data.frame()
  final_result_comparison.df <- data.frame()
  plot.gg <- vector(mode = 'list', length = length(variables))
  
  
  #check if the independent variable has more than 2 factor levels
  if (group_factor_length > 2){
    #run analysis for each variable i.e otu
    for(index in seq_along(variables)) {
      if(grepl(pattern = "\\-|;|\\s+|\\(|\\)", x = variables[index])){
        old_y <- variables[index]
        #correct names that have alternatives by deleting the alternative
        variables[index] <- gsub(pattern ="/.+",replacement = "",x =variables[index])
        variables[index] <- gsub(pattern ="\\-|;|\\s+|\\(|\\)",replacement = "_",x =variables[index])
        #rename the column with the acceptable name
        colnames(OTU_table.df)[which(colnames(OTU_table.df) == old_y)] <- variables[index]
      }
      #check for a problematic column name then correct it
      res <- try(parse(text = variables[index]), silent=T)
      if(inherits(res,'try-error')){
        colnames(OTU_table.df)[index] <- make.names(names = variables[index])
        variables[index] <- make.names(names = variables[index])
      }
      #treatments analysis
      #perform one-way anova on Treatments
      var_anovaResult <- anova_analysis(x = group, y = variables[index], Data = OTU_table.df)
      #compare the anova groups
      try(
        compare_treat <- compare_groups(x = group, y = variables[index], Data = OTU_table.df)
        , silent = T)
      #combine the anova result for the treatment analysis for each variable
      result.df <- rbind(result.df,var_anovaResult)
      
      #combine the treatment comparison results for each variable
      final_result_comparison.df <- rbind(final_result_comparison.df, compare_treat)
      
      
      #draw the treatment plot
      invisible(print(plot.gg[[index]]<-try(draw_plot(x = group, 
                                                      y = variables[index],
                                                      xLab = group,
                                                      yLab = 'Relative abundance(%)', 
                                                      Data = OTU_table.df,
                                                      date = sprintf('%s %s', MAtrix, season), 
                                                      plot_type = 'box'), silent = T)))
    }
    
  }else{
    #if it is a two factor variable then loop or each otu perform 2 sample tests for each variable or OTU
    for(index in seq_along(variables)){
      if(grepl(pattern = "\\-|;|\\s+|\\(|\\)", x = variables[index])){
        old_y <- variables[index]
        variables[index] <- gsub(pattern ="\\-|;|\\s+|\\(|\\)",replacement = "_",x =variables[index])
        #rename the column with the acceptable name
        colnames(OTU_table.df)[which(colnames(OTU_table.df) == old_y)] <- variables[index]
      }
      #check for a problematic column name then correct it
      res <- try(parse(text = variables[index]), silent=T)
      if(inherits(res,'try-error')){
        colnames(OTU_table.df)[index] <- make.names(names = variables[index])
        variables[index] <- make.names(names = variables[index])
      }
      #2 factor sample t test and wilcoxon analysis
      try(
        var_T_test <- tTest_analysis(x = group, y = variables[index], Data = OTU_table.df)
        , silent = T)
      #combine the water type result for each variable
      result.df <- rbind(result.df,var_T_test)
      
      #draw the water type plot
      invisible(print(plot.gg[[index]]<-try(draw_plot(x = group, 
                                                      y = variables[index],
                                                      xLab = group,
                                                      yLab = 'Relative abundance(%)', 
                                                      Data = (OTU_table.df),
                                                      date = sprintf('%s %s', MAtrix, season), 
                                                      plot_type = 'box'), silent = T)))
      
      
    }
  }
  if (group_factor_length > 2){
    #combine the anova result
    result.df$fdr_Anova_pval <- p.adjust(p =result.df[,'anova_pval'] ,pAjust)
    result.df$fdr_Kruskal_pval <- p.adjust(p =result.df[,'kruskal_pval'] ,pAjust)
    final_result.df <- cbind(Parameter = variables, result.df)
    
    
    
    #rename the column names    
    colnames(final_result.df) <- c('Taxon', 'group_DF', 'residual_DF', 'anova_fVal',
                                   'anova_pVal','shapiro_pVal', 'levene_pVal', 'Krus_statistic',
                                   'Kruskal_pVal', sprintf('%s_Anova_pval',pAjust), sprintf('%s_Kruskal_pval',pAjust))
    
    colnames(final_result_comparison.df) <- c('Taxon', 'group_1','group_2','pVal',
                                              'holm_pVal','p.format','p.signif','method')
    
    #name the plots
    names(plot.gg) <- variables
    
    #order based on the fdr adjusted p-values
    #final_result.df <- final_result.df[order(final_result.df$fdr_Kruskal_pval),]
    
    #create a final list containing all the results
    Final_result <- list(final_result.df,final_result_comparison.df,plot.gg)
    
    #give the parts of the final result names
    names(Final_result) <- c('multiple_group_test', 'comparison','graphs')
    
    
  }else{
    #combine the water type result for the current subset
    result.df$fdr_tTest_pval <- p.adjust(p =result.df[,'tTest_pval'] ,pAjust)
    result.df$fdr_Wilcox_pval <- p.adjust(p =result.df[,'wilcox_pval'] ,pAjust)
    final_result.df <- cbind(Parameter = variables, result.df)
    
    #rename the column names 
    colnames(final_result.df) <- c('Taxon', 'shapiro_pVal', 'Var_pVal',
                                   'tTest_statistic','tTest_pVal', 'wilcox_statistic',
                                   'wilcox_pVal',sprintf('%s_tTest_pval',pAjust), sprintf('%s_Wilcox_pval',pAjust))
    #name the plots
    names(plot.gg) <- variables
    
    Final_result <- list(final_result.df,plot.gg)
    
    #give the parts of the final result names
    names(Final_result) <-c('two_group_test','graphs')
  }
  
  return(Final_result)
  
}

#function to perform single or multiway anova with assumptions testing
factorial_anova <- function(formula,Data,y,is_one_predictor=TRUE,...){
  
  # 1. formula is a formula
  # 2. Data is a dataframe containing the variables to be analysed
  # 3. y is a string specifying the response variable
  # 4. is_one_predictor a boolean specifying if you have only one predictor variable or not
  # 4. ... any parameter that can be passed to aov or lm functions
  
  # #############    Usage ##################
  # # For one predictor variable
  # result<- factorial_anova(formula = reformulate(termlabels = "growth_medium*sterility", response = "Acinetobacter"),
  #                          Data = new_data, y = "Acinetobacter", is_one_predictor = TRUE)
  # 
  # # For multiple predictor variables
  # result<- factorial_anova(formula = reformulate(termlabels = "growth_medium+sterility+growth_medium*sterility",
  #                                                response = "Acinetobacter"),
  #                          Data = new_data, y = "Acinetobacter", is_one_predictor = FALSE)
  # # Get the results
  # result$anova_test
  # result$assumptions_test
  
  
  
  if ((length(unique(Data[,y])) > (length(Data[,y])/20))){
    indicator.aov <-aov(formula=formula,data=Data,...)
    tidy_aov <- broom::tidy(indicator.aov)
    
    ifelse(test = length(unique(residuals(indicator.aov))) > (length(residuals(indicator.aov))/2),
           yes = shapiro_pval <- (shapiro.test(residuals(indicator.aov)))$p.value, no = shapiro_pval<-0)
    
    
    
  }else{
    shapiro_pval <- NA
    Levene_pval <- NA
  }
  
  if(is_one_predictor){
    Levene_pval<- DescTools::LeveneTest(lm(formula = formula, data = Data))$'Pr(>F)'[1]
    assumptions_test <- data.frame(shapiro_pval = shapiro_pval,
                                   Levene_pval = Levene_pval )
  }else{
    assumptions_test <- data.frame(shapiro_pval = shapiro_pval)
  }
  
  results <- list(tidy_aov,assumptions_test)
  names(results) <- c("anova_test","assumptions_test")
  return(results)
}

get_treat_top25 <- function(function_profile,metadata, group){
  common.ids <- intersect(rownames(function_profile),rownames(metadata))
  metadata <- metadata[common.ids,]
  treat_function_profile <- function_profile[common.ids,]
  
  treat_function_profile <- cbind(subset(x = metadata, select=group),treat_function_profile)
  
  treat_function_profile <- aggregate(formula= reformulate(termlabels = group, response = '.'),
                                      data = treat_function_profile, FUN = sum)
  
  rownames(treat_function_profile) <- treat_function_profile[,group]
  treat_function_profile <- treat_function_profile[,-which(colnames(treat_function_profile) == group)]
  # convert to relative abundance 
  treat_function_profile <- t(apply(X = treat_function_profile, MARGIN =  1, FUN = function(x) x/sum(x)))
  
  if(ncol(treat_function_profile) >= 25){
    treat_top_25 <- sort(apply(X = treat_function_profile, MARGIN =  2, FUN = max),decreasing = TRUE )[1:25]
    treat_top_25 <- treat_function_profile[,names(treat_top_25 )]
  }else{
    treat_top_25 <- treat_function_profile
  }
  
  # Rename the columns after with names that exclude the KEGG and EC ids 
  colnames(treat_top_25) <- gsub(pattern = '.+;\\s|\\s\\[.+\\]',replacement = '', colnames(treat_top_25))
  
  return(treat_top_25)
}

#function to collapse the samples in an oTU table with a defined function(fun)  based on a group in metadata 
collapse_samples <- function(taxon_table,metadata,group,fun=sum, convertToRelativeAbundance=FALSE){
  #function to collapse the samples in an oTU table with a defined function(fun)  based on a group in metadata 
  #taxon_table - a matrix count table with samples as rows and features/OTUs as columns
  #metadata - a dataframe to containing the group to collapse samples by. Sample names must be the rownames of the metadata
  #group - an independent factor variable within the metadata to collapse the samples by
  #fun - a function without brackets to apply in order to collapse the samples
  # convertToRelativeAbundance - a boolean set to TRUE OR FALSE if the taxon_table shout be converted to relative abundance
  # default is FALSE
  common.ids <- intersect(rownames(taxon_table),rownames(original_metadata))
  metadata <- droplevels(original_metadata[common.ids,,drop=FALSE])
  taxon_table <- taxon_table[common.ids,,drop=FALSE]
  taxon_table <- cbind(subset(x = metadata, select=group),taxon_table)
  
  taxon_table <- aggregate(formula= reformulate(termlabels = group, response = '.'),
                           data = taxon_table, FUN = fun)
  rownames(taxon_table) <- taxon_table[,1]
  taxon_table <- taxon_table[,-1]
  if(convertToRelativeAbundance){
    taxon_table <- t(apply(X = taxon_table, MARGIN = 1, FUN = function(x) x/sum(x)))
  }
  
  final <- list(taxon_table,metadata)
  names(final) <- c("taxon_table","metadata")
  return(final)
}

# subsample feature table to even depth
rarefaction <- function(ps,depth=0.9*min(sample_sums(ps))) {
  #1. ~ ps is a phyloseq object
  #2. ~ depth is the refaction depth e.g 0.9*min(sample_sums(ps)), min(sample_sums(ps))
  #    a number of your choice e.g 19449
  
  #rarefy the otu table
  ps.rarefied <- rarefy_even_depth(physeq = ps, 
                                   sample.size = depth,
                                   rngseed = 1, 
                                   replace = F, 
                                   verbose = FALSE)
  
  
  
  #plot rarefaction curve
  rare_gg <- rarecurve(t(otu_table(ps.rarefied)), step=50, cex=0.5)
  #print(rare_gg)
  result <- list(ps.rarefied, rare_gg)
  names(result) <- c("ps.rarefied","rare_curve")
  return(result)
}


# Alpha diversity
alpha_diversity <- function(ps.rarefied, metadata, group='treatments', 
                            metrics=c("Observed", "Shannon"),
                            group_order=NULL){
  
  #ps.rarefied - rarefied phyloseq object
  # metatdata - metadata table describing each sample
  # group - independent variable within you metadat for the results should be grouped
  # group_order - a vector specifying the or within the group for plotting
  # metrics - a vector specifying the alpha diversity metrices to plot
  
  
  #get all  sorts of alpha diversity measurements 
  diversity.df <- estimate_richness(ps.rarefied,
                                    measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
  
  #get the samples present in both the metadata and diversity table
  common.ids <- intersect(rownames(diversity.df), rownames(metadata))
  diversity.df <- diversity.df[common.ids,]
  sample_names <- rownames(diversity.df)
  metadata <- as.data.frame(metadata[common.ids,])
  rownames(metadata) <- sample_names
  colnames(metadata) <- group
  
  
  diversity.df <- data.frame(metadata,diversity.df)
  
  #mean and standard deviation of the diversity matrices
  diversity_meanSd <- diversity.df %>%  
    select(-c(se.chao1,se.ACE))%>% 
    group_by(!!sym(group))%>% 
    summarise_all(list(mean=mean,std=sd))
  
  #number of samples per treatment group 
  number_of_samples <- diversity.df %>%  
    select(-c(se.chao1,se.ACE)) %>% 
    group_by(!!sym(group))  %>% 
    summarise(N= n())
  
  
  #function to calculate standard error
  SE<- function(x){
    sd(x)/sqrt(length(x))
  }
  
  diversity_se <- diversity.df %>%  
    select(-c(se.chao1,se.ACE)) %>% 
    group_by(!!sym(group))  %>% 
    summarise_all(list(se=SE))
  
  
  
  
  #re-order the levels of the grouping factor if need be
  if(is.null(group_order)){
    sample_data(ps.rarefied) <- metadata
  }else{
    metadata[,group] <- factor(x =metadata[,group], levels = group_order)
    # Add metadata to the phyloseq object
    sample_data(ps.rarefied) <- metadata
  }
  
  
  #plot alpha diversity
  richness_gg <- plot_richness(ps.rarefied,
                               x=group,
                               measures=metrics) + geom_boxplot()
  
  result <- list(diversity.df,metadata, number_of_samples, diversity_meanSd,diversity_se,richness_gg)
  names(result) <- c("diversity.df","metadata", "N", "mean_sd","SE","plot")
  return(result)
}



findGenesInGff<- function(gff_path, genePattern='cob|cib|blu|hem|cysG|sirC',
                          file_ext='_genomic.gff.gz', removeInName ='_'){
  
  clusterProfiler::Gff2GeneTable(gffFile = gff_path)
  
  Name <- BiocGenerics::basename(gff_path) %>% 
    str_replace_all(file_ext, '') %>% 
    str_replace_all(removeInName, ' ')
  
  load('geneTable.rda')
  
  gene_df <- geneTable %>% 
    filter(str_detect(GeneName, genePattern)) %>% 
    arrange(GeneName)
  
  geneNames <- gene_df %>% pull(GeneName) %>% str_c(collapse = ',')
  gene_df <- gene_df  %>% mutate(OrganismName=strrep(Name, times = nrow(gene_df)))
  file.remove('geneTable.rda')
  
  return(list(genes=geneNames, genesDF=gene_df))
  
}


# Calcultae Aitchison distances between samples
# devtools::install_github("robCompositions", "matthias-da")
# library('robCompositions') 
# https://stackoverflow.com/questions/47777626/error-in-code-calculating-aitchison-distance
aDistFun <- function(x, r){
  # Example:
  # X <- iris [ ,1:3]
  # a <- sapply(seq(1, nrow(X)), FUN=aDistFun, x=X )
  return(aDist( x[r-1, ], x[r, ])) 
}

create_dt <- function(table2show,...){
  DT::datatable(table2show,
                rownames = FALSE, # remove row numbers
                filter = "top", # add filter on top of columns
                extensions = "Buttons", # add download buttons
                options = list(
                  autoWidth = TRUE,
                  dom = "Blfrtip", # location of the download buttons
                  buttons = c("copy", "csv", "excel", "pdf", "print"), # download buttons
                  pageLength = 5, # show first 5 entries, default is 10
                  order = list(0, "asc") # order the title column by ascending order
                ),
                escape = FALSE # make URLs clickable) , ...
  )
  
}



process_metaphlan_taxonomy <- function(taxonomy, prefix='\\w__') {
  #function to process a metaphlan2 taxonopmy assigment table
  #1. ~ file_path is a string specifying the taxonomic assignment file name
  #2 prefix ~ is a regular expression specifying the characters to remove
  # from the taxon names  '\\w__'  for greengenes and 'D_\\d__' for SILVA
  
  #taxon_levels <- c("kingdom","phylum","class","order","family","genus","species", "strain")
  
  taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 
  
  #taxonomy[,'species'] <- paste(taxonomy[,'genus'],taxonomy[,'species'])
  # replace NAa with Other and delete the D_num__ prefix from the taxonomy names
  for (rank in colnames(taxonomy)) {
    #delete the taxonomy prefix
    taxonomy[,rank] <- gsub(pattern = prefix, x = taxonomy[, rank],replacement = '')
    indices <- which(is.na(taxonomy[,rank]))
    taxonomy[indices, rank] <- rep(x = "Other", times=length(indices)) 
    #replace empty cell
    indices <- which(taxonomy[,rank] == "")
    taxonomy[indices,rank] <- rep(x = "Other", times=length(indices))
  }
  taxonomy <- apply(X = taxonomy,MARGIN = 2,
                    FUN =  gsub,pattern = "_",replacement = " ") %>% 
    as.data.frame(stringAsfactor=F)
  #taxonomy <- as.data.frame(taxonomy,stringsAsFactors=F)
  return(taxonomy)
}


#function for format a taxonomy assignment table by appending sufix to a known name
format_taxonomy_table <- function(taxonomy=taxonomy.m,stringToReplace="Other", suffix=";Other") {
  
  for (taxa_index in seq_along(taxonomy)) {
    #indices <- which(taxonomy[,taxa_index] == stringToReplace)
    
    indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
    
    taxonomy[indices,taxa_index] <- 
      paste0(taxonomy[indices,taxa_index-1],
             rep(x = suffix, times=length(indices)))
    
  }
  return(taxonomy)
}


fix_names<- function(taxonomy,stringToReplace,suffix){
  #1~ taxonomy is a taxonomy dataframe with taxonomy ranks as column names
  #2~ stringToReplace is a vector of regex strings specifying what to replace
  #3~ suffix is a string specifying the replacement value
  
  
  for(index in seq_along(stringToReplace)){
    taxonomy <- format_taxonomy_table(taxonomy = taxonomy,
                                      stringToReplace=stringToReplace[index], 
                                      suffix=suffix[index])
  }
  return(taxonomy)
}



combine_previous_and_current_cell <- function(df,current_col_index=8,suffixes, suffixes_indices){
  
  stopifnot(length(suffixes) == length(suffixes_indices), is.numeric(current_col_index))
  
  df[suffixes_indices,current_col_index] <-  paste0(df[suffixes_indices,current_col_index-1],suffixes)
  
  return(df)
  
}


# A function to generate taxon level count matrix based on a taxonomy table and an existing feature table
make_feature_table <- function(count_matrix,taxonomy,taxon_level){
  
   # EAMPLE:
   # make_feature_table(count_matrix = feature_counts_matrix,taxonomy = taxonomy_table,taxon_level = "Phylum")
  
  feature_counts_df <- data.frame(taxon_level=taxonomy_table[,taxon_level],
                                  count_matrix, check.names = FALSE, stringsAsFactors = FALSE)
  
  feature_counts_df <- aggregate(formula=.~taxon_level,data = feature_counts_df,FUN = sum)
  rownames(feature_counts_df) <- feature_counts_df[,"taxon_level"]
  return(feature_counts_df[,-1])
}


get_rarefaction_df <- function(sample_index, rare_curve_list) {
      # rare_curve_list - phloseq rarecurve list generated after runnig alpha_diversity function
  x <- attr(rare_curve_list[[sample_index]], which = "Subsample")
  sample_name <- names(x)[length(x)]
  data.frame(sample=sample_name, sequences=x, observed_species=rare_curve_list[[sample_index]])
  
}

