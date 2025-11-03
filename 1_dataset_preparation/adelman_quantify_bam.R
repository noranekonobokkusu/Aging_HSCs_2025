suppressMessages( {
  library(Rsubread)
  library(dplyr)
} )


args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outfile_gene = paste0('gene.', args[2])

# Quantify on the level of genes
t = featureCounts(infile,
                  annot.ext = "/broad/vangalenlab/vangalen/Genomes/GRCh38.221223/GRCh38/genes/genes.gtf.gz",
                  isGTFAnnotationFile = T,
                  GTF.attrType = "gene_name",
                  useMetaFeatures = T,
                  genome = "/broad/vangalenlab/vangalen/Genomes/GRCh38.221223/GRCh38/fasta/genome.fa",
                  isPairedEnd = T,
                  strandSpecific = 0,
                  nthreads = 4)

results = data.frame(t$annotation[,c("GeneID", "Length")], Counts = t$counts, stringsAsFactors=F)
colnames(results)[3] = "Counts"
write.table(x = results, outfile_gene, quote = F, sep = "\t", row.names = F)
