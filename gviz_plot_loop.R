#!/usr/bin/env Rscript
# load library
library("tidyverse")
library("Gviz")
library("GenomicRanges")
library(biomaRt)
source("sub_chr_ensembl.function.R")
source("sub_strand_gviz.function.R")
# import data
tdt <- as_tibble(read.csv("../TDT/asd.290.IBD_cleaned.tdt", sep=""))
as_tibble(tdt) %>%
  arrange(P) %>%
  mutate(CHR=sub_chr_ensembl(CHR))%>%
  mutate(LOG=-log10(P)) %>%
  mutate_at("CHR", ~gsub("^","chr",.)) -> tdt
recomb.rate <- read.delim("~/Tools/Library/Genetic_map/recomb-hg38/genetic_map_GRCh38_merged.tab")
# hg38; ensembl database.
gene.ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        dataset="hsapiens_gene_ensembl")
# loop through significant snp
for (i in 1:4) {
  print(paste("Plotting variant",i,"using gviz"))
  #
  tdt.sel <- tdt %>% dplyr:::slice(i)
  sel.chr <- tdt.sel$CHR
  sel.pos <- tdt.sel$BP
  range <- 5e5
  gen <- "hg38"
  # create chromosome ideogram and genome axis track
  gtrack <- GenomeAxisTrack(labelPos="below")
  itrack <- IdeogramTrack(genome=gen, chromosome= sel.chr)
  # create adjacent variants Pvalue track
  tdt %>%
    filter(CHR==sel.chr, between(BP, sel.pos-range, sel.pos+range)) %>%
    dplyr::select(c(1,3,11)) -> tdt.sel.region
  # convert to suitable format for gviz
  grange.tdt <- makeGRangesFromDataFrame(tdt.sel.region, keep.extra.columns=TRUE,
                                         ignore.strand=TRUE,
                                         seqnames.field="CHR",
                                         start.field="BP",
                                         end.field="BP")
  dtrack <- DataTrack(grange.tdt, name="-log10(P)", genome=gen,
                      baseline=-log10(0.05/nrow(tdt)),
                      col.baseline="red",lty.baseline="dashed")
  displayPars(dtrack) <- list(ylim=c(0,9))
  # create annotation tracks of variants
  tdt %>%
    filter(CHR==sel.chr, between(BP, sel.pos-range, sel.pos+range)) %>%
    dplyr::select(c(1,2,3)) -> sel.variant
  sel.variant[-1,2] <- " "
  sel.variant %>%
    dplyr::rename(chromosome=1, id=2, start=3) %>%
    mutate(end=start,group=c("sel",rep("adj",nrow(sel.variant)-1))) -> sel.variant
  # make transparent blue color
  col0 <- rgb(0, 0, 255, max = 255, alpha = 100, names = "blue40")
  # atrack with feature colored
  atrack <- AnnotationTrack(name="Variants", genome=gen, chromosome=sel.chr,
                            start=sel.variant$start, end=sel.variant$end,
                            id=sel.variant$id,featureAnnotation="id",
                            fontcolor.feature="darkblue",
                            rotation.title=0,showTitle=TRUE,cex.title=0.5,
                            shape="box",stacking="dense", #below is new code
                            feature=rep(c("selected","adjacent"),c(1,nrow(sel.variant)-1)),
                            col="transparent", selected="red", adjacent=col0
                            )
  # create gene region track
  # query needed biomart field
  out.bm.genes.region <- getBM(
    attributes = c('chromosome_name','exon_chrom_start','exon_chrom_end','strand',
                   'gene_biotype',
                   'ensembl_gene_id','ensembl_exon_id','ensembl_transcript_id',
                   'external_gene_name'), 
    filters = c('chromosome_name','start','end'), 
    values = list(gsub("chr","",sel.chr),sel.pos - range, sel.pos + range), 
    mart = gene.ensembl)
  # reformat dataframe for plotting
  out.bm.genes.region %>%
    dplyr::rename(chromosome=1,
                  start=2,
                  end=3,
                  strand=4,
                  feature=5,
                  gene=6,
                  exon=7,
                  transcript=8,
                  symbol=9) %>%
    mutate(strand=sub_strand_gviz(strand)) %>%
    filter(feature=="protein_coding") %>%
    mutate(symbol=paste0(symbol,strand)) %>%
    mutate_at("symbol",~gsub("\\+$","\u2192",.)) %>%
    mutate_at("symbol",~gsub("\\-$","\u2190",.)) -> genes.region
  grtrack <- GeneRegionTrack(genes.region,name="Known.Genes\nEnsembl",
                             genome=gen,chromosome=sel.chr,
                             transcriptAnnotation="symbol",
                             collapseTranscripts="longest",
                             fontfamily.group="Arial Unicode MS")
  # plot recombination rate track
  recomb.rate %>%
    filter(chrom==sel.chr, between(pos, sel.pos-range, sel.pos+range)) -> recomb.rate.sel
  recomb.rate.sel %>%
    arrange(pos) %>%
    mutate(pos_end=c(pos[-1],pos[nrow(recomb.rate.sel)])) %>%
    dplyr::select(1,2,5,3) -> recomb.df
  grange.recomb <- makeGRangesFromDataFrame(recomb.df, keep.extra.columns=TRUE,
                                            ignore.strand=TRUE,
                                            seqnames.field="chrom",
                                            start.field="pos",
                                            end.field="pos_end")
  rrtrack <- DataTrack(grange.recomb, name="Recomb.\nrate\n(cM/Mbp)", genome=gen,
                       type="l")
  displayPars(rrtrack) <- list(ylim=c(0,100))
  # highlight track
  ht <- HighlightTrack(trackList=list(atrack, grtrack,dtrack,rrtrack,gtrack),
                       start=sel.pos-2000, width=4000,
                       chromosome=sel.chr,
                       col="#FFE3E6")
  # plot with custom track size and save in png file
  output=paste0("2variant_",i,".",sel.chr,"_",sel.pos,".png")
  png(file=output,
      bg="white",units="in",width=8, height=4,res=300)
  plotTracks(list(itrack, ht), sizes=c(1,1,3,5,3,2),
             from=sel.pos-range, to=sel.pos+range)
  dev.off()
  # end
  print(output)
}
