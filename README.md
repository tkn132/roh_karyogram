# ROH karyogram
Plot ROH density on the chromosome set

```bash
awk '((NR==1) || (($5+0)>0)){print $1,$3,$5}' pilot.hom.summary > pilot.hom.summary.simplified
```



```R
library(GenomicRanges)
library(regioneR)
library(GenomeInfoDb) 
library(karyoploteR)

## read and collapse the roh data.
dat2 = read.table("pilot.hom.summary.simplified", header=T)
head(dat2)

out = data.frame("chr"=NA, "start"=NA, "end"=NA, "nsnps"=NA, "unaff"=NA)

for (i in unique(dat2$CHR)){
sdat = dat2[dat2$CHR == i, ]
group = rle(sdat$UNAFF) #Remove/collapse consecutive duplicate values in sequence
nsnps= group$lengths #number of snps is the lengths of each value
index.start = cumsum(c(1, group$lengths[-length(group$lengths)])) 
index.end = c(index.start[-1] - 1, nrow(sdat))

start = sdat$BP[index.start]
end = sdat$BP[index.end]
head(end-start)
res = data.frame(rep(i,length(start)),start, end, nsnps, group$values)
names(res) = names(out)
out = rbind(out, res)
}
out=out[-1,]
write.csv(out, "../roh/dat/pilot.hom.summary.collapsed.csv", row.names=F, quote=F)
```

```R
library(ggbio)
library(GenomicRanges)
library(regioneR)
library(GenomeInfoDb) 
library(RColorBrewer)

dat = read.csv("dat/pilot.hom.summary.collapsed.csv", header=T)
head(dat)

## convert data into a Grange object
dat$chr = paste0("chr", dat$chr)
head(dat)

dat <- toGRanges(dat, genome="hg38")
head(dat)

## trim all but autosomal chroms
seqlevels(dat) = paste0("chr", 1:22)
seqlengths(dat)

## roh frequency
dat$freq = dat$unaff/371

## color
palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')

## plot karyogram
sdat = dat[dat$unaff>2,]


png('karyogram.unaff3.png', width=1200, height=1200, res=150)
autoplot(sdat, layout="karyogram", aes(fill=freq)) +
  labs(fill="Frequency") + 
  scale_fill_continuous(trans= 'reverse')  + 
  scale_fill_gradientn(colours = palette(100))
dev.off()
```




