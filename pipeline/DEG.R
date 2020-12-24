suppressMessages({
if (!require("optparse")) {install.packages("optparse")}
})
# construct args
option_list <- list(
  make_option(c("-m", "--meta"), action="store", default=NA, type='character',
              help="group meta information"),
  make_option(c("-e", "--expression"), action="store", default=NA, type='character',
              help="expression file meta. per line: sample file_path"),
  make_option(c("-a", "--annotation"), action="store", default=NA, type='character',
              help="any other annotations. First column should be gene name"),
  make_option(c("-c", "--compare"), action="store", default=NA, type='character',
              help="Group to be compared. Two groups per line"),
  make_option(c("-k", "--kegg"), action="store", default=NA, type='character',
              help="KEGG annotation"),
  make_option(c("-g", "--go"), action="store", default=NA, type='character',
              help="GO annotation"),
  make_option(c("-t", "--gotable"), action="store", default=NA, type='character',
              help="GO table"),
  make_option(c("-o", "--outdir"), action="store", default=NA, type='character',
              help="output directory.") 
)

suppressMessages({
if (!require("BiocManager")) {install.packages("BiocManager")}
if (!require("DESeq2")) {BiocManager::install("DESeq2")}
if (!require("BiocParallel")) {install.packages("BiocParallel")}
if (!require("clusterProfiler")) {BiocManager::install("clusterProfiler")}
if (!require("dplyr")) {install.packages("dplyr")}
if (!require("magrittr")) {install.packages("magrittr")}
if (!require("tximport")) {BiocManager::install("tximport")}
})

opt <- parse_args(OptionParser(option_list=option_list))
# use multiple cores
register(MulticoreParam(4))
# load data
expression_meta <- read.delim(opt$expression, header=FALSE)
expression_files <- expression_meta$V2 %>% `names<-`(expression_meta$V1)
rcm <- tximport(expression_files, type = "rsem", txIn = FALSE, txOut = FALSE)
colmeta <- read.delim(opt$meta, header = FALSE, row.names = 1, stringsAsFactors = TRUE)
colmeta <- colmeta[match(colnames(rcm$counts), rownames(colmeta)),] %>% 
  as.data.frame() %>% 
  `rownames<-`(colnames(rcm$counts)) %>% 
  `colnames<-`("V2")
groupmeta <- read.delim(opt$compare, col.names = c('group1', 'group2'), header = FALSE)
go_annotation <- read.delim(opt$go) %>% 
  mutate(level = case_when(level=="biological_process" ~ "BP",
                           level=="cellular_component" ~ "CC",
                           level=="molecular_function" ~ "MF")
         )
go_info <- read.delim(opt$gotable) %>% 
  mutate(level = case_when(level=="biological_process" ~ "BP",
                           level=="cellular_component" ~ "CC",
                           level=="molecular_function" ~ "MF")
  )
kegg_annotation <- read.delim(opt$kegg)
oannotation <- read.delim(opt$annotation, row.names = 1) %>% mutate(row_name=row.names(.))
# DESeq2
dds <- DESeqDataSetFromTximport(txi = rcm,
                                colData = colmeta,
                                design = ~ V2)
dds <- DESeq(dds)
# enrich and plot function

# get results
## make output directory
if (!dir.exists(opt$outdir)){dir.create(opt$outdir)}

for (idx in 1:nrow(groupmeta)) {
  compare_name = paste0(groupmeta[idx,1], "VS", groupmeta[idx,2])
  dir.create(paste(opt$outdir, compare_name, sep='/'))
  res <- results(dds, contrast = c("V2", groupmeta[idx,1], groupmeta[idx,2]))
  write.table(res, paste(opt$outdir, compare_name, 'gene.txt', sep='/'), sep='\t', quote = FALSE, col.names = NA)
  subset_index <- (res$padj<0.05) & (abs(res$log2FoldChange)>2)
  subset_index[is.na(subset_index)] <- FALSE
  dematrix <- res[subset_index,]
  if (nrow(dematrix) == 0){
    next
  } else {
    dematrix %<>% 
      as.data.frame() %>%  
      mutate(row_name=row.names(.)) %>% 
      left_join(rcm$counts %>% 
                  as.data.frame() %>%
                  dplyr::select( rownames(colmeta)[colmeta$V2 %in% groupmeta[idx,]] ) %>% 
                  mutate(row_name=row.names(.))
                ) %>% 
      left_join(oannotation) %>%
      relocate("row_name") %>%
      `rownames<-`(.$row_name) %>%
      select(-row_name)
    write.table(dematrix, paste(opt$outdir, compare_name, 'DEG.txt', sep='/'), sep='\t', quote = FALSE, col.names = NA)
    ## kegg enrichment
    gene_list <- dematrix$log2FoldChange %>% `names<-`(rownames(dematrix))
    kegg_enrich <- enricher(names(gene_list), 
                            pAdjustMethod = "fdr", 
                            TERM2GENE = kegg_annotation[c(3,1)], 
                            TERM2NAME = kegg_annotation[c(3,4)])
    
    if (!(is.null(kegg_enrich) || nrow(kegg_enrich@result)==0)){
      write.table(kegg_enrich@result, paste(opt$outdir, compare_name, 'KEGG_enrich.txt', sep='/'), sep='\t', quote = FALSE, row.names = FALSE)
      ### dot plot
      pdf(paste(opt$outdir, compare_name, 'dotplot.pdf', sep='/'), width = 8, height = 10)
      dotplot(kegg_enrich, showCategory=30)
      dev.off()
      ### heatplot
      pdf(paste(opt$outdir, compare_name, 'heatplot.pdf', sep='/'), width = 15, height = 3)
      heatplot(kegg_enrich, foldChange=gene_list)
      dev.off()
      ### cnetplot
      pdf(paste(opt$outdir, compare_name, 'cnetplot.pdf', sep='/'), width = 10, height = 10)
      cnetplot(kegg_enrich)
      dev.off()
    }
    
    ## go enrichment
    for ( lev in c('BP','CC','MF')) {
      go_enrich <- enricher(gene_list, 
                            pAdjustMethod = "fdr", 
                            TERM2GENE = go_annotation[go_annotation$level==lev, c(2,1)], 
                            TERM2NAME = go_info[go_info$level==lev, c(1,2)])
      if (!(is.null(go_enrich) || nrow(go_enrich@result)==0)){
        write.table(kegg_enrich@result, paste(opt$outdir, compare_name, paste0(lev, '_GO_enrich.txt'), sep='/'), sep='\t', quote = FALSE, row.names = FALSE)
        ### dot plot
        pdf(paste(opt$outdir, compare_name, paste0(lev,'dotplot.pdf'), sep='/'), width = 8, height = 10)
        dotplot(kegg_enrich, showCategory=30)
        dev.off()
        ### heatplot
        pdf(paste(opt$outdir, compare_name, paste0(lev,'heatplot.pdf'), sep='/'), width = 15, height = 3)
        heatplot(kegg_enrich, foldChange=gene_list)
        dev.off()
        ### cneplot
        pdf(paste(opt$outdir, compare_name, paste0(lev,'cnetplot.pdf'), sep='/'), width = 10, height = 10)
        cnetplot(kegg_enrich)
        dev.off()
      }
    }
  }
}
