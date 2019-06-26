##################
# Tangles 
# by Hannes Becher
##################

# Script allows you to generate tangle plots from files containing read mapping
#  positions and qualities. The input file has to be a gz-compressed text file
#  with three columns: read name, mapping position, and mapping quality (BWA).
#  The file "N1_names.gz" in the folder "data" is an example.
# To create suitable input files, map your data in single-end mode with BWA,
#  sort the resulting files by name, and extract the read names, positions, and
#  mapping qualities with samtools and GNU command line utilities:
#  for i in *bam; do samtools view $i | cut -f 1,4,14 | gzip > ${i%.*}_names.gz; done



# Utility functions ####
# A function to seperate singlet reads from pairs, returns a list of (usually)
#  three elements: file path, singletons, and pairs
sing.pair = function(filepath) {
  # file connections
  print(paste0("Reading file: ", filepath))
  con = gzcon(file(filepath, open="rb"))
  prf <- numeric(0)
  prfq <- numeric(0)
  prr <- numeric(0)
  prrq <- numeric(0)
  si <- numeric(0)
  siq <- numeric(0)
  # read first line
  old = strsplit(readLines(con, n=1), "\t")[[1]]
  purged=F
  lc = 1
  ls = 0
  ll = 0
  # loop over remaining lines
  while ( TRUE ) {
    if(lc %% 10000 == 0) print(paste0("Lines read: ", lc))
    lc <- lc + 1
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) { # enf of file, finishing steps
      if(purged==F){ # there is a singleton in old, save it
        ls <- ls + 1
        si[ls]  <- old[2]
        siq[ls] <- substr(old[3],6,9)       
      } # the alternative is purged=T, the pair would already have been written in that case
      break
    }
    # to do normally
    nw = strsplit(line, "\t")[[1]]
    if(purged==F){
      if(nw[1] == old[1]){# pair found
        ll <- ll + 1
        prf[ll]  <- old[2]
        prfq[ll] <- substr(old[3],6,9)
        prr[ll]  <- nw[2]
        prrq[ll] <- substr(nw[3],6,9)
        purged = T
      } else { # old was a single read
        ls <- ls + 1
        si[ls]  <- old[2]
        siq[ls] <- substr(old[3],6,9)
        
      }
      
    } else { # purged last time round
      purged = F
    }
    # this happens in every iteration
    old = nw 
  }
  print(paste0("Single reads found: ", length(si)))
  print(paste0("Pairs: ", length(prf)))
  close(con)
  si=data.frame(si=as.numeric(si), si.q=as.numeric(siq))
  print("Done.")
  print("")
  if(length(prf)>0){
    pr=data.frame(pr.f=as.numeric(prf), pr.fq=as.numeric(prfq), pr.r=as.numeric(prr), pr.rq=as.numeric(prrq))
    prd=abs(pr[,1]-pr[,3])
    return(list(file=filepath, si=si, pr=data.frame(pr, pr.d=prd)))
  } else {
    
    return(list(file=filepath, si=si))
  }
  
}

# A function to make a tangle plot, takes as input the output of sing.pair().
#  Adjust ll to match the length of your mapping reference. The arguments "nam"
#  and "cl" are passed to the function plot() as "main" and "col".
plot.is <- function(x, ll=16008, nam="", cl = 1){
  plot(c(-1.1,1.1), c(-1.1,1.1), asp=1, type="n", main = nam, frame.plot = F, xlab="", ylab="", axes = F)
  # print("Background plotted.")
  y <- x$pr[,c(1,3)]
  xx <- y/ll *2 * pi
  
  sc <- cbind(sin(xx[,1]), cos(xx[,1]), sin(xx[,2]), cos(xx[,2]))
  for(i in 1:nrow(sc)){
    
    points(c(sc[i,1], sc[i,3]), c(sc[i,2], sc[i,4]), type = "l", col = cl)
  }
}


#  Get data ####

setwd("~/git_repos/tangles/data/")


# A tangle plot example ####
# sing.pair() groups reads by whethere they are paired or singletons
n1 <- sing.pair("N1_names.gz")
# There are three elements in the list: the file name, a data.frame of singleton reads, and
#  a data/frame of paired reads. The data frams have columns for mapping posistion (si/pr.f/pr.r),
#  mapping score (si.q/pr.fq/pr.fr), and the distance between the paired reads mapped (pr.d).
str(n1)


# remove short-insert pairs and artefacts
n1$pr <- n1$pr[n1$pr[,5] > 1499 & n1$pr[,5] < 14501, ]


# Plot ####
# Use transparence when specyfying the colour of tangles,
#  overlapping regions will appear darker.
plot.is(n1, nam="My first Tangle Plot", cl = "#0000FF60")



# multiple samples ####
# this requires you to first extract the read positions from the BAM files
#  supplied in the folder "data". These files were generated with BWA and were
#  sorted by read name with samtools. In the directory with the BAM files,
#  run (in one line):
# for i in *bam; do samtools sort  -n -@ 4 -o ${i%.*}_s.bam $i; done
#  and then:
# for i in *s.bam; do samtools view $i | cut -f 1,4,14 | gzip > ${i%_*}_names.gz; done

#setwd("~/git_repos/tangles/data/")
files <- dir(pattern = "names.gz") # get file names

# should print this: "N1_names.gz" "N2_names.gz" "N3_names.gz" "N4_names.gz" "N5_names.gz" "N6_names.gz"
files

# read files looping over file names with 'lapply'
dat <- lapply(files, function(x) {
  sing.pair(x)
})

# remove pairs with insert size < 1500 and > 14500
dat <- lapply(dat, function(x) {
  a <- x
  a$pr <- a$pr[a$pr[,5] > 1499 & a$pr[,5] < 14501,]
  a
})


# Reproduce the figure in the paper.
par(mfrow=c(2,3))

# 6 plots
lapply(dat, function(x){
  plot.is(x, nam = strsplit(x$file,"_")[[1]][[1]], cl="#0000FF40")
  # plot an outline (optional)
  pos <- seq(0, 2*pi, 2*pi/1000) 
  lines(sin(pos),cos(pos), lwd = 2)
  # # add some numbers (optional)
  nums <- 16
  pos20 <- seq(0, 2*pi, 2*pi/nums)[1:nums] + 2*pi/nums/2
  text(sin(pos20)*1.1,cos(pos20)*1.1, labels=1:nums)
})
par(mfrow=c(1,1))


# Example: Plotting with the Circlize package ####
library(circlize)

# get annotations
gff <- read.table("../dupl11_cons.fasta_mitos.gff")
l <- levels(gff[,3])
circos.par("track.height" = 0.05, start.degree=90, gap.degree=0)
circos.initialize(c("a", "a"),c(1, 16008))
circos.track("a", ylim = c(0,2), bg.border=0)
labs <- sapply(substr(gff[,9],6,20), function(x) strsplit(x,"[(]")[[1]][1])
# plot gene names excluding tRNAs
circos.text(apply(gff[,4:5],1,mean)[gff[,3] !="tRNA"], 3, labs[gff[,3] !="tRNA"], "a",1)

apply(gff[gff[,3] == l[1],4:5], 1, function(x) circos.lines(x, c(0,0), lwd=2, col = "red"))
apply(gff[gff[,3] == l[2],4:5], 1, function(x) circos.lines(x, c(0,0), lwd=2, col = "green"))
apply(gff[gff[,3] == l[3],4:5], 1, function(x) circos.lines(x, c(0,0), lwd=2, col = "blue"))

apply(dat[[1]]$pr, 1, function(x){
  circos.link("a", x[1], "a", x[3], w=0, col="#0000FF60")
})
# dev.off() # reset to remove circlize setting



# Example: Matrices, distences, trees from tangles####

# First the mapping positiions are binned. The binning can be adjusted with
# the seq() functions. Here, we have bins of 250 bp, going a bit further than
# the mito genome's length (16500 instead of 16008).

# using pairs
pair.bins500 <- lapply(dat, function(x){
  a <- cbind(.bincode(x$pr[,1], seq(0, 16500,250)),
             .bincode(x$pr[,3], seq(0, 16500,250)))
  
  t(apply(a, 1, function(y) if(y[1] > y[2]) c(y[2], y[1]) else y))
  
})


# matrices
pair.counts500 <- as.data.frame(sapply(pair.bins500, function(x) {
  m <- max(unlist(pair.bins500))
  mat <- matrix(0, m, m)
  for(i in 1:nrow(x)){
    mat[x[i, 1], x[i, 2]] <- mat[x[i, 1], x[i, 2]] + 1
  }
  mat[upper.tri(mat)]
  }))
names(pair.counts500) <- sapply(dat, function(x) strsplit(x$file, "_")[[1]][[1]])
pair.counts.norm500 <- apply(pair.counts500, 2, function(x) x/sum(x, na.rm = T))

plot(hclust(dist(t(pair.counts.norm500), method="manhattan")),
     main="Cluster dendrogram based on tangle patterns")





# test area ####
str(pair.bins500)
str(pair.counts500)


rM1 <- sing.pair("/media/hannes/2nd/Euphrasia/rDNA/M1rDNA_names.gz")
rpodi <- sing.pair("/media/hannes/2nd/Podi_rDNA/novoP/PodiN2_names.gz")
rSim <- list()
rSim$pr <- data.frame(pr.f = sample.int(10000,1372, replace = T),
                      pr.fq = sample.int(10000,1372, replace = T),
                      pr.r = sample.int(10000,1372, replace = T))

str(rM1)
str(rpodi)
str(rSim)

# remove short-insert pairs and artefacts
rM1$pr <- rM1$pr[rM1$pr[,5] > 1499 & rM1$pr[,5] < 8500, ]
rpodi$pr <- rpodi$pr[rpodi$pr[,5] > 1499 & rpodi$pr[,5] < 23500, ]
str(rM1)
str(rpodi)
hist(rM1$pr[,5])
plot.is(rM1, ll=10000, nam="M1 rDNA", cl = "#0000FF20")
plot.is(rpodi, ll=25000, nam="podi rDNA", cl = "#0000FF60")
hist(rpodi$si$si, breaks=100)

plot.is(rSim, ll=10000, cl = "#0000FF20")
