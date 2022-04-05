suppressMessages({
  if (!require("optparse")) {install.packages("optparse")}
  if (!require("BiocManager")) {install.packages("BiocManager")}
  if (!require("DESeq2")) {BiocManager::install("DESeq2")}
  if (!require("clusterProfiler")) {BiocManager::install("clusterProfiler")}
  if (!require("dplyr")) {install.packages("dplyr")}
  if (!require("magrittr")) {install.packages("magrittr")}
})

# construct args
option_list <- list(
  make_option(c("-m", "--meta"), action="store", default=NA, type='character',
              help="group information. Two columns:geneid/group"),
  make_option(c("-e", "--expression"), action="store", default=NA, type='character',
              help="read counts file"),
  make_option(c("-a", "--annotation"), action="store", default=NA, type='character',
              help="any other annotations. First column should be gene name"),
  make_option(c("-c", "--compare"), action="store", default=NA, type='character',
              help="Group to be compared. Two groups per line"),
  make_option(c("-k", "--kegg"), action="store", default=NA, type='character',
              help="KEGG annotation. Two columns: geneid/KO num"),
  make_option(c("-g", "--go"), action="store", default=NA, type='character',
              help="GO annotation. Two columns: geneid/GO num"),
  make_option(c("-o", "--outdir"), action="store", default=NA, type='character',
              help="output directory.")
)

opt <- parse_args(
  OptionParser(
    usage = "usage: %prog [options]",
    prog = "PEx",
    option_list=option_list,
    description = "
    All flies (except read counts matrix) should not have headers
    ")
  )

get_this_file <- function(){
  commandArgs() %>% 
    tibble::enframe(name=NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
}

# load data
read_counts <- read.delim(opt$expression, header=TRUE, row.names = 1) %>% round()
colmeta <- read.delim(opt$meta, header=FALSE, row.names = 1, stringsAsFactors=TRUE)
colmeta <- colmeta[match(colnames(read_counts), rownames(colmeta)),] %>% 
  as.data.frame() %>% 
  `rownames<-`(colnames(read_counts)) %>%
  `colnames<-`("V2")
groupmeta <- read.delim(opt$compare, col.names = c('group1', 'group2'), header = FALSE)
script_dir <-dirname(get_this_file())
load(paste0(script_dir, "/DB.rdata"))
go_annotation <- read.delim(opt$go,  header = FALSE) %>% `names<-`(c('geneid', 'GO')) %>% left_join(go_info[c(1,3)])
kegg_annotation <- read.delim(opt$kegg, header = FALSE) %>% `names<-`(c('geneid', 'KO'))
kegg_annotation %<>% left_join(ko2pathway)
if(!is.na(opt$annotation)){
  oannotation <- read.delim(opt$annotation, row.names = 1, header = FALSE) %>% mutate(geneid=row.names(.))
}
# DESeq2
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = colmeta,
                              design = ~ V2)
dds <- DESeq(dds)
# enrich and plot function

# get results
## make output directory
if (!dir.exists(opt$outdir)){dir.create(opt$outdir)}

for (idx in 1:nrow(groupmeta)) {
  compare_name = paste0(groupmeta[idx,1], "VS", groupmeta[idx,2])
  if (!exists(paste(opt$outdir, compare_name, sep='/'))){dir.create(paste(opt$outdir, compare_name, sep='/'))}
  res <- results(dds, contrast = c("V2", groupmeta[idx,1], groupmeta[idx,2]))
  write.table(res, paste(opt$outdir, compare_name, 'gene.txt', sep='/'), sep='\t', quote = FALSE, col.names = NA)
  subset_index <- (res$padj<0.05) & (abs(res$log2FoldChange)>1)
  subset_index[is.na(subset_index)] <- FALSE
  dematrix <- res[subset_index,]
  if (nrow(dematrix) == 0){
    next
  } else {
    dematrix %<>% 
      as.data.frame() %>%  
      mutate(geneid=row.names(.)) %>% 
      left_join(read_counts %>% 
                  dplyr::select( rownames(colmeta)[colmeta$V2 %in% groupmeta[idx,]] ) %>% 
                  mutate(geneid=row.names(.))
      ) %>% 
      relocate("geneid") %>%
      `rownames<-`(.$geneid)
      
    if (exists("oannotation")){
        dematrix %<>% left_join(oannotation)
      }
    
    for ( direction in c('up', 'down')){
      out_dir <- paste(opt$outdir, compare_name, direction, sep = '/')
      if (!dir.exists(out_dir)){dir.create(out_dir)}
      ## kegg enrichment
      if (direction == 'up'){
        degmatrix <- dematrix %>% filter(log2FoldChange>0)
      }else{
        degmatrix <- dematrix %>% filter(log2FoldChange<0)
      }
      write.table(degmatrix, paste(out_dir, 'DEG.txt', sep='/'), sep='\t', quote = FALSE, row.names = FALSE)
      gene_list <- degmatrix$log2FoldChange %>% `names<-`(degmatrix$geneid)
      kegg_enrich <- enricher(names(gene_list), 
                              pAdjustMethod = "fdr", 
                              TERM2GENE = kegg_annotation[c(3,1)], 
                              TERM2NAME = pathway_3levels[c(1,2)])
    
      if (!(is.null(kegg_enrich) || nrow(kegg_enrich@result)==0)){
        write.table(kegg_enrich@result, paste(out_dir, 'KEGG_enrich.txt', sep='/'), sep='\t', quote = FALSE, row.names = FALSE)
        ### dot plot
        if (!sum(kegg_enrich@result$p.adjust<0.05)==0){
        dotplot(kegg_enrich, showCategory=30)
        ggplot2::ggsave(paste(out_dir, 'dotplot.pdf', sep='/'), width = 8, height = 10)
        ### heatplot
        heatplot(kegg_enrich, foldChange=gene_list)
        ggplot2::ggsave(paste(out_dir, 'heatplot.pdf', sep='/'), width = 15, height = 3)
        ### cnetplot
        cnetplot(kegg_enrich)
        ggplot2::ggsave(paste(out_dir, 'cnetplot.pdf', sep='/'), width = 10, height = 10)
        }
      }
    
    ## go enrichment
      for ( lev in c('BP','CC','MF')) {
        go_enrich <- enricher(names(gene_list), 
                              pAdjustMethod = "fdr", 
                              TERM2GENE = go_annotation[go_annotation$level==lev, c(2,1)], 
                              TERM2NAME = go_info[go_info$level==lev, c(1,2)])
        if (!(is.null(go_enrich) || nrow(go_enrich@result)==0)){
          write.table(go_enrich@result, paste(out_dir, paste0(lev, '_GO_enrich.txt'), sep='/'), sep='\t', quote = FALSE, row.names = FALSE)
          if (!sum(go_enrich@result$p.adjust<0.05)==0){
          ### dot plot
          dotplot(go_enrich, showCategory=30)
          ggplot2::ggsave(paste(out_dir, paste0(lev,'dotplot.pdf'), sep='/'), width = 8, height = 10)
          ### heatplot
          heatplot(go_enrich, foldChange=gene_list)
          ggplot2::ggsave(paste(out_dir, paste0(lev,'heatplot.pdf'), sep='/'), width = 15, height = 3)
          ### cneplot
          cnetplot(go_enrich)
          ggplot2::ggsave(paste(out_dir, paste0(lev,'cnetplot.pdf'), sep='/'), width = 10, height = 10)
          }
        }
      }
    }
  }
}