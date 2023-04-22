library(GenomicRanges)
library(seqinr)
library(stringr)
library(ggplot2)
library(seqRFLP)

### starting files: MACS2 peak-calling files ###
### get the intersection and union of 3' End seq peaks from two repeats ###
# macs2 callpeak -t 777chip-1_S2_aln.rmdup.fwd.sort.bam --nomodel -f BAMPE -g 4.6e6 -n 0730c.fwd -B -q 0.00001
# macs2 callpeak -t 777chip-1_S2_aln.rmdup.rev.sort.bam --nomodel -f BAMPE -g 4.6e6 -n 0730c.rev -B -q 0.00001

orient=c("fwd")
dir=c("/home/jingjing/Desktop/ChIP analyses2020/")

input1=paste(dir,"0730t0730c.",orient,"_peaks",sep="")
input2=paste(dir,"0225t0730c.",orient,"_peaks",sep="")
c777=paste(dir, "0730c.",orient,"_peaks",sep="")

rep1=read.table(input1,header=T,sep="",stringsAsFactors = FALSE)
rep2=read.table(input2,header=T,sep="",stringsAsFactors = FALSE)
control=read.table(c777,header=T,sep="",stringsAsFactors = FALSE)
peaks1=GRanges(seqnames=rep1$chr,ranges=IRanges(start=rep1$start,end=rep1$end))
peaks2=GRanges(seqnames=rep2$chr,ranges=IRanges(start=rep2$start,end=rep2$end))
peaksc=GRanges(seqnames=control$chr,ranges=IRanges(start=control$start,end=control$end))

### intermediate files ###
tmp1=rep1[queryHits(findOverlaps(peaks1,peaksc,minoverlap=25L)),]
tmp2=control[subjectHits(findOverlaps(peaks1,peaksc,minoverlap=25L)),]
write.table(cbind(tmp1,tmp2),paste(dir, "0730andControlOverlap.",orient,".txt",sep=""),
                                   quote=FALSE,row.names=FALSE,sep="\t")

tmp1=rep2[queryHits(findOverlaps(peaks2,peaksc,minoverlap=25L)),]
tmp2=control[subjectHits(findOverlaps(peaks2,peaksc,minoverlap=25L)),]
write.table(cbind(tmp1,tmp2),paste(dir, "0225andControlOverlap.",orient,".txt",sep=""),
            quote=FALSE,row.names=FALSE,sep="\t")

### exclude non-specific region ###

rep1.trun=rep1[-queryHits(findOverlaps(peaks1,peaksc,minoverlap=25L)),]
rep2.trun=rep2[-queryHits(findOverlaps(peaks2,peaksc,minoverlap=25L)),]
peaks1=GRanges(seqnames=rep1.trun$chr,ranges=IRanges(start=rep1.trun$start,end=rep1.trun$end))
peaks2=GRanges(seqnames=rep2.trun$chr,ranges=IRanges(start=rep2.trun$start,end=rep2.trun$end))

rep1_overlap=queryHits(findOverlaps(peaks1,peaks2))
rep2_overlap=subjectHits(findOverlaps(peaks1,peaks2))
common=cbind(rep1.trun[rep1_overlap,],rep2.trun[rep2_overlap,])

### extract genomic sequence +- 75bp from peak summit ###

refseq=getSequence(read.fasta(file="/home/jingjing/Desktop/NC_000913\ ref\ files/NC_000913.3.fa"))
frag.fa=cbind(name=common[,10],seq=NA)
for (i in 1:nrow(frag.fa)){
  frag.fa[i,2]=paste(getFrag(refseq[[1]],common[i,5]-550,common[i,5]-50),collapse="")
}
write.fasta(sequences = as.list(frag.fa[,2]), names = frag.fa[,1], 
            file.out = paste(dir,"common_summit_up500bp_",orient,".fa",sep=""))


### MEME command line ###
#meme control_summit150bp_rev.fa -dna -mod anr -revcomp -minw 5 -maxw 20 
### align MEME-analyzed motif and peak files ###
meme=read.table(paste(dir,"common_meme_",orient,sep=""),header=T,sep="",stringsAsFactors = FALSE)
meme$count=sapply(meme$Sequence_name,function(x){sum(meme$Sequence_name==x)})
meme.rmdup=meme[-which(duplicated(meme$Sequence_name)=="TRUE"),]
common.motif=cbind(common,SequenceName=NA,motifStart=NA,Upstream=NA,motifSite=NA,Downstream=NA,MotifCount=NA)
for (i in 1:nrow(common.motif))
{
  tmp=common.motif[i,10]
  pos=which(meme.rmdup$Sequence_name==tmp)
  if(length(pos)>0){
    common.motif$SequenceName[i] = meme.rmdup$Sequence_name[pos]
    common.motif$motifStart[i] = meme.rmdup$Start[pos]
    common.motif$Upstream[i] = meme.rmdup$Upstream[pos]
    common.motif$motifSite[i] = meme.rmdup$Site[pos]
    common.motif$Downstream[i] = meme.rmdup$Downstream[pos]
    common.motif$MotifCount[i] = meme.rmdup$count[pos]
  }
}
write.table(common.motif, paste(dir,"common_motif_",orient,sep=""),quote=FALSE,row.names=FALSE,sep="\t")

### focusing on simple-peak signal ###
tmp=common.motif[,c(1,5,10,22,23,24,25,26)]
common.info=tmp[which(tmp$MotifCount==1),]
common.info=cbind(common.info,motifCoordinate.start=common.info$abs_summit-75+common.info$motifStart-1,
                  motifCoordinate.end=common.info$abs_summit-75+common.info$motifStart+6,gene=NA, 
                  strand=NA, fiveGene=NA, fiveGeneStrand=NA, threeGene=NA,threeGeneStrand=NA,
                  DirectRNA.fwd.start=NA, DirectRNA.fwd.end=NA, DirectRNA.rev.start=NA, DirectRNA.rev.end=NA)

### gene.pos and DirectRNA peak pos ###
dir2=c("/home/jingjing/Documents/20201109 Termseq/")

DirectRNA=list(fwd=NA, rev=NA)
for (i in 1:2){
  file=paste(dir2,"DirectRNA_seq_",names(DirectRNA)[i],sep="")
  peak=read.table(file,header=F,sep="",stringsAsFactors = FALSE, skip=1)[,c(1:3)]
  colnames(peak)=c("chr","start","end")
  DirectRNA[i]=GRanges(seqnames=rep("NC_000913.3",nrow(peak)),ranges=IRanges(start=peak$start,end=peak$end))
}

genes.table=read.table("/media/jingjing/JINGJING/ChIP-seq analyses/MG1655v3_genes_MochiView.txt",header=T,sep="",stringsAsFactors = FALSE)
genes=GRanges(seqnames=genes.table$SEQ_NAME,ranges=IRanges(start=genes.table$START,end=genes.table$END))

query=GRanges(seqnames=common.info$chr,ranges=IRanges(start=common.info$motifCoordinate.start,end=common.info$motifCoordinate.end))
for (i in 1:nrow(common.info)){
  
  ### gene annotate###
  tmp=subjectHits(findOverlaps(query[i],genes))
  if (length(tmp) > 0) 
  {
    common.info$gene[i]=genes.table$FEATURE_NAME[tmp]
    common.info$strand[i]=genes.table$STRAND[tmp]
  }else{
    common.info$gene[i]=c("intergenic")
    common.info$strand[i]=c("intergenic")
  }
  if (orient == "fwd" ) {
    tmp1=follow(query[i],genes)
    common.info$fiveGene[i]=genes.table$FEATURE_NAME[tmp1]
    common.info$fiveGeneStrand[i]=genes.table$STRAND[tmp1]
    tmp2=precede(query[i],genes)
    common.info$threeGene[i]=genes.table$FEATURE_NAME[tmp2]
    common.info$threeGeneStrand[i]=genes.table$STRAND[tmp2]
  }else{
    tmp1=precede(query[i],genes)
    common.info$fiveGene[i]=genes.table$FEATURE_NAME[tmp1]
    common.info$fiveGeneStrand[i]=genes.table$STRAND[tmp1]
    tmp2=follow(query[i],genes)
    common.info$threeGene[i]=genes.table$FEATURE_NAME[tmp2]
    common.info$threeGeneStrand[i]=genes.table$STRAND[tmp2]
    }

  ### DirectRNA-seq signal ###
  tmp=subjectHits(findOverlaps(query[i],DirectRNA$fwd))
  if (length(tmp) > 0){
    common.info$DirectRNA.fwd.start[i]=start(DirectRNA$fwd)[tmp]
    common.info$DirectRNA.fwd.end[i]=end(DirectRNA$fwd)[tmp]
  }
  tmp=subjectHits(findOverlaps(query[i],DirectRNA$rev))
  if (length(tmp) > 0){
    common.info$DirectRNA.rev.start[i]=start(DirectRNA$rev)[tmp]
    common.info$DirectRNA.rev.end[i]=end(DirectRNA$rev)[tmp]
  }
}

write.table(common.info, paste(dir,"common_singlePeakInfo.",orient,sep=""),quote=FALSE,row.names=FALSE,sep="\t")


### read rpoB2 peak files ###
tmp=common.motif[,c(1,2,3,10,26)]
wt_sp=tmp[which(tmp$MotifCount==1),]
wtpeaks=GRanges(seqnames=wt_sp$chr,ranges=IRanges(start=wt_sp$start,end=wt_sp$end))

input1=paste(dir,"greB29t0730c.",orient,"_peaks",sep="")
input2=paste(dir,"greB30t0730c.",orient,"_peaks",sep="")

mut1=read.table(input1,header=T,sep="",stringsAsFactors = FALSE)
mut2=read.table(input2,header=T,sep="",stringsAsFactors = FALSE)

mut1peaks=GRanges(seqnames=mut1$chr,ranges=IRanges(start=mut1$start,end=mut1$end))
mut2peaks=GRanges(seqnames=mut2$chr,ranges=IRanges(start=mut2$start,end=mut2$end))

### intermediate files ###
tmp1=mut1[queryHits(findOverlaps(mut1peaks,wtpeaks,minoverlap=25L)),]
tmp2=mut2[queryHits(findOverlaps(mut2peaks,wtpeaks,minoverlap=25L)),]
write.table(rbind(tmp1,tmp2),paste(dir, "greB_detectedSP.",orient,".txt",sep=""),
            quote=FALSE,row.names=FALSE,sep="\t")

tmp1=rpoB2[queryHits(findOverlaps(rpoB2peaks,peaksc,minoverlap=25L)),]
tmp2=control[subjectHits(findOverlaps(rpoB2peaks,peaksc,minoverlap=25L)),]
write.table(cbind(tmp1,tmp2),paste(dir, "rpoB2andControlOverlap.",orient,".txt",sep=""),
            quote=FALSE,row.names=FALSE,sep="\t")

### exclude non-specific region ###

rpoB1.trun=rpoB1[-queryHits(findOverlaps(rpoB1peaks,peaksc,minoverlap=25L)),]
rpoB2.trun=rpoB2[-queryHits(findOverlaps(rpoB2peaks,peaksc,minoverlap=25L)),]
rpoB1peaks=GRanges(seqnames=rpoB1.trun$chr,ranges=IRanges(start=rpoB1.trun$start,end=rpoB1.trun$end))
rpoB2peaks=GRanges(seqnames=rpoB2.trun$chr,ranges=IRanges(start=rpoB2.trun$start,end=rpoB2.trun$end))

rep1_overlap=queryHits(findOverlaps(rpoB1peaks,rpoB2peaks))
rep2_overlap=subjectHits(findOverlaps(rpoB1peaks,rpoB2peaks))
common.rpoB2=cbind(rpoB1.trun[rep1_overlap,],rpoB2.trun[rep2_overlap,])

B2pos=GRanges(seqnames=common.rpoB2[,1],ranges=IRanges(start=common.rpoB2[,2],end=common.rpoB2[,3]))
common.rpoB2.info=cbind(common.rpoB2, sponPeakName=NA, summit=NA)
for (i in 1:nrow(common.rpoB2)){
  tmp=subjectHits(findOverlaps(B2pos[i],query))
  if (length(tmp)>0){
    common.rpoB2.info$sponPeakName[i]=common.info$name[tmp]
    common.rpoB2.info$summit[i]=common.info$abs_summit[tmp]
  }
}

write.table(common.rpoB2.info,paste(dir, "GamcombinecI.info.",orient,".txt",sep=""),
            quote=FALSE,row.names=FALSE,sep="\t")


#~/Downloads/bowtie-1.2.3-linux-x86_64/bowtie -v 0 -a MG1655v3 -f bowtie_query.fa >motifAllPos --un fail.fa

supple=read.table(c("/home/jingjing/Desktop/suppleTable1"),header=T,sep="",stringsAsFactors = FALSE)
supple.pos=GRanges(seqnames=supple$chr,ranges=IRanges(start=supple$start,end=supple$end))
tmp=supple[queryHits(findOverlaps(supple.pos,query)),]
write.table(tmp,"/home/jingjing/Desktop/suppleTable1",quote=FALSE,row.names=FALSE,sep="\t")
common.exd=cbind(common,motifRev.start=NA,motifRev.end=NA)
tmp=common.exd
query=GRanges(seqnames=tmp[,1],ranges=IRanges(start=tmp[,2],end=tmp[,3]))

### overlap with motifSites ###
file=paste(dir,"MACS2/MotifSite_",orient,".bed",sep="")
motif=read.table(file,header=F,sep="",stringsAsFactors = FALSE)
colnames(motif)=c("chr","motif.start","motif.end","strand","seq")
motifList=GRanges(seqnames=motif$chr,ranges=IRanges(start=motif$motif.start,end=motif$motif.end))

for (i in 1:nrow(common.exd)){
  tmp=subjectHits(findOverlaps(query[i],motifList))
  if (length(tmp) > 0){
    common.exd$motifRev.start[i]=start(motifList)[tmp]
    common.exd$motifRev.end[i]=end(motifList)[tmp]
  }
}
                     

dir2=c("/home/jingjing/Documents/20201109 Termseq/")
motifList=list(fwd=NA,rev=NA)
for (i in 1:2){
  file=paste(dir2,"MACS2/MotifSite_",names(motifList)[i],".bed",sep="")
  motif=read.table(file,header=F,sep="",stringsAsFactors = FALSE)
  colnames(motif)=c("chr","motif.start","motif.end","strand","seq")
  motifList[i]=GRanges(seqnames=motif$chr,ranges=IRanges(start=motif$motif.start,end=motif$motif.end))
} 

union_unique=common
query=GRanges(seqnames=union_unique$chr,ranges=IRanges(start=union_unique$start,end=union_unique$end))
union_unique_ext=cbind(union_unique,motifFwd.start=NA,motifFwd.end=NA, motifRev.start=NA, motifRev.end=NA)
for (i in 1:nrow(union_unique)){
  tmp=subjectHits(findOverlaps(query[i],motifList$fwd))
  if (length(tmp) > 0){
    union_unique_ext$motifFwd.start[i]=start(motifList$fwd)[tmp]
    union_unique_ext$motifFwd.end[i]=end(motifList$fwd)[tmp]
  }
  tmp=subjectHits(findOverlaps(query[i],motifList$rev))
  if (length(tmp) > 0){
    union_unique_ext$motifRev.start[i]=start(motifList$rev)[tmp]
    union_unique_ext$motifRev.end[i]=end(motifList$rev)[tmp]
  }
}

write.table(union_unique_ext, "/media/jingjing/JINGJING/ChIP-seq analyses/rpoB2_union_unique_peaks_motif.txt",
            quote=FALSE,row.names=FALSE,sep="\t")

genes.table=read.table("/media/jingjing/JINGJING/ChIP-seq analyses/MG1655v3_genes_MochiView.txt",header=T,sep="",stringsAsFactors = FALSE)
genes=GRanges(seqnames=genes.table$SEQ_NAME,ranges=IRanges(start=genes.table$START,end=genes.table$END))
query=GRanges(seqnames=rep("MG1655v3",nrow(union_unique)),ranges=IRanges(start=union_unique$abs_summit,end=union_unique$abs_summit+1))
output=cbind(union_unique,genes=NA, strand=NA)

for (i in 1:nrow(union_unique))
{
  tmp=subjectHits(findOverlaps(query[i],genes))
  if (length(tmp) > 0) 
  {
    output$genes[i]=genes.table$FEATURE_NAME[tmp]
    output$strand[i]=genes.table$STRAND[tmp]
  }else{
    output$genes[i]=c("intergenic")
    output$strand[i]=c("intergenic")
  }
}

write.table(output,"/media/jingjing/JINGJING/ChIP-seq analyses/union_unique_peaks.txt",
            quote=FALSE,row.names=FALSE,sep="\t")



### exclude all antisense motif sites located inside the 3'end profile peak ###
### get positions of all profile peaks ###
all_peaks_pos=GRanges(seqnames=all_peaks$chr,ranges=IRanges(start=all_peaks$start,end=all_peaks$end))
all_motif_plus=read.table("/media/jingjing/JINGJING/ChIP-seq analyses/allExIDsite_antisense_motif_200bp_plusgene.bed.txt",
                          header=F,sep="",stringsAsFactors = FALSE)
plus_pos=GRanges(seqnames=all_motif_plus$V1,ranges=IRanges(start=all_motif_plus$V2,end=all_motif_plus$V3))
plus_overlap=subjectHits(findOverlaps(all_peaks_pos,plus_pos))
write.table(all_motif_plus[-plus_overlap,], "/media/jingjing/JINGJING/ChIP-seq analyses/allExIDsite_antisense_motif_200bp_plusgene_nonpeak.bed",
            quote=FALSE,row.names=FALSE,sep="\t")

all_motif_minus=read.table("/media/jingjing/JINGJING/ChIP-seq analyses/allExIDsite_antisense_motif_200bp_minusgene.bed.txt",
                           header=F,sep="",stringsAsFactors = FALSE)
minus_pos=GRanges(seqnames=all_motif_minus$V1,ranges=IRanges(start=all_motif_minus$V2,end=all_motif_minus$V3))
minus_overlap=subjectHits(findOverlaps(all_peaks_pos,minus_pos))
write.table(all_motif_minus[-minus_overlap,], "/media/jingjing/JINGJING/ChIP-seq analyses/allExIDsite_antisense_motif_200bp_minusgene_nonpeak.bed",
            quote=FALSE,row.names=FALSE,sep="\t")

### count reads at stalling ends of extracted NET bam file ###
plus.end=read.table("/home/jingjing/Documents/20201024 NET-seq analysis/plus.sort.trun.bed",
                    header=F,sep="",stringsAsFactors = FALSE)[,3]
plus.freq=as.data.frame(table(plus.end)/22.4)
plus.freq$plus.end=as.numeric(levels(plus.freq$plus.end))
minus.end=read.table("/home/jingjing/Documents/20201024 NET-seq analysis/minus.sort.trun.bed",
                    header=F,sep="",stringsAsFactors = FALSE)[,2]
minus.freq=as.data.frame(table(minus.end)/16.38)
minus.freq$minus.end=as.numeric(levels(minus.freq$minus.end))

ggplot(plus.freq,aes(x=plus.end,y=Freq))+
  geom_point()+
  xlim(1954774,1954799)+
  ylim(0,1)
  
ggplot(minus.freq,aes(x=minus.end,y=Freq))+
  geom_point()+
  xlim(4231135,4231160)+
  ylim(0,5)


### plot MotifLoc ###

### motif length ###
len=7

### motif color primary is polyT ###
motifcolor=c("#d80000","#007e24")

### motif color primary is polyA ###
motifcolor=c("#007e24","#d80000")

temp=read.table(paste(pwd,"/MEMEplot.txt",sep=""),header=T,sep="",stringsAsFactors = FALSE)
peaks.table=temp

peaks.table[,3]=as.factor(peaks.table[,3])
motif=cbind(peaks.table,MotifEnd=peaks.table[,4]+len)

ID=1:23
ggplot(motif, aes(x=MotifStart,y=y,xend=MotifEnd,yend=y,color =Strand))+
  geom_segment()+
  scale_color_manual(values=motifcolor)+
  theme_bw()+
  scale_x_continuous(limits=c(0,200),expand=c(0.01,0.01))+
  scale_y_continuous("ID", labels = as.character(ID), breaks = ID,trans = "reverse", expand=c(0.01,0.01))+
  theme(axis.text=element_text(size=rel(1.4)),
#       axis.text.y=element_blank(),
#       axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


