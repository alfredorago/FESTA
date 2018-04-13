# Convert h5 files into tsv
# This script collects the h5v files produced by kallisto and the sample metadata
# Both are used to generate the count table of exons and the exon and sample annotation

# Clean workspace and print date
rm(list=ls())
date()

# Load libraries
library(tximport)
library(stringr)
library(DESeq2)

# Select source data directory
dataSource <- file.path('../Results/20180413/SIRV_only')

# Set output path
outdir = file.path("../Results", format(Sys.Date(), format = "%Y%m%d"), "SIRV_only/tximport")
dir.create(outdir, recursive = T)

# Set path of data files
files = list.dirs(path = dataSource, full.names = F)
files = grep(pattern = '^[A-Z][0-9]{2}$', x = files, value = T)
files = file.path(dataSource, files , 'abundance.h5')

### Create exon to gene reference table from first sample
tpm = read.table(file = file.path(dataSource, 'F58/abundance.tsv'), header = T, row.names = 1)
IDtable = data.frame(exonID = row.names(tpm))
# Create gene names as substrings: SIRV and one numeric character after the start (WIP)
IDtable$geneID = sapply(X = IDtable$exonID, FUN = function(x) {
  if (grepl(pattern = 'SIRV', x = x)) {
    str_extract(string = x, pattern = 'SIRV.?')
  }
}
)

### Import sample annotation and subset to sequenced samples
sampleData = read.csv('../Results/20180322/Metadata_compiler/sample_metadata.csv',
                  header = T, row.names = 1)
sampleData = sampleData[which(row.names(sampleData)%in%str_extract(string = files, pattern = '[A-Z][0-9]{2}')),]
sampleData$stage = as.factor(sampleData$stage)
names(sampleData)[3] = "mix"

### Import counts and convert to integers for DESeq analyses
txTranscript = tximport(files = files, type = 'kallisto', tx2gene = IDtable,
                        abundanceCol = 'tpm', lengthCol = 'length',
                        txOut = T)
deseqTranscript = DESeqDataSetFromTximport(txi = txTranscript,
                                           colData = sampleData,
                                           design = ~ sirv)
### Save counts and assignment table for further analyses
SIRV_data = list(tpm = counts(deseqTranscript), ID = IDtable[,c('geneID','exonID')], annot = sampleData[,'mix'])
save(SIRV_data, file = file.path(outdir, 'SIRV_data.RData'), compress = 'gzip')
