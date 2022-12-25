library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPQC)

bamfile = "sortbams/SRR502329.bam";
bedfile = "macs/STAT1_30m_IFNa_summits.bed";
bedfile2 = "macs/STAT1_6h_IFNa_summits.bed";

# dict1 = {"SRR502329": "STAT1_30m_IFNa","SRR502327":"STAT1_6h_IFNa","SRR502228":"INP_30m_IFNa","SRR502225":"INP_6h_IFNa"}

reads = readGAlignments(bamfile); # Returns GAlignments object
#bed <- as.data.frame(read.table(bedfile,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")) 
peaks <- import.bed(bedfile, genome = "hg19")
peaks2 <- import.bed(bedfile2, genome = "hg19")
# Peaks are loaded.


bam_views <- BamViews(file, bamRanges = GRanges("chr20", IRanges(start = 10000, end = 110000))) #We may use BamViews 
# to define regions of interest in our reads data (for ex, a particular region on a particular chromosome)

reads <- readGAlignments(bam_views) # Import again

bams <- BamViews(bamfile, bamRanges = peaks) # Use peaks to define views into BAM files
reads <- readGAlignments(bamfile)


# Data Plotting

# Determined by where a lot of things were mapped for chromosome Y
starting =   13121902 
ending = 13801106
ideogram <- IdeogramTrack("chrY", "hg19")
axis <- GenomeAxisTrack()
tx <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, chromosome = "chrY", start = starting, end = ending, name = "Genes")

# Creating data tracks - these are pairs (supposedly numeric, genome interval), with other useful things we care about defined.
# I am guessing it does filtering based on the other arguments by itself.
peaksDT <- DataTrack(peaks, chromosome = "chrY", start = starting, end = ending, name = "Peaks", type = "h")
peaksDT2 <-  DataTrack(peaks2, chromosome = "chrY", start = starting, end = ending, name = "Peaks", type = "h")
#covDT <- DataTrack(coverage(reads)[["chrY"]], chromosome = "chrY", start = starting, end = ending, name = "Coverage", type = "h") - Figure how to plot coverage
blacklist_go <- AnnotationTrack(blacklist, name = "Blacklist", id = blacklist@elementMetadata$name)
plotTracks(list(ideogram, axis, peaksDT2, peaksDT, tx, blacklist_go), chromosome = "chrY", from = starting, to = ending, featureAnnotation = "id")

# Can we use reads data - perhaps to see how reads map?






# Removing blacklisted peaks
blacklisted <- findOverlaps(peaks, blacklist.hg19, type="within")
clean_peaks <- peaks[-from(blacklisted)]







# Visualising mapq scores
reads <- readGAlignments(bamfile, param = ScanBamParam(what="mapq")) # Load reads with mapping qualities by explicitly mentioning it.

table(mcols(reads)$mapq) # This reveals that I don't have mapq scores available. Oh well.










# Quality Control with ChIPQC - there's so much more you can do. Check out their manual!
BlackListFile <- ("blacklist/ENCFF001TDO.bed")
exp <- ChIPQCsample(bamfile,peaks = bedfile,annotation = "hg19",blacklist = BlackListFile, verbose = FALSE)
QCmetrics(exp) # Print the metrics

plotFrip(exp) # fraction read in peaks (inside peaks vs. outside peaks)
plotFribl(exp); # Plot fraction reads in blacklist (inside black and outside black)





# Computing coverage
reads_gr <- granges(reads) # Extract genomic ranges of reads object.

frag_length <- fragmentlength(exp) # Obtain average fragment length.

reads_ext <- resize(reads_gr, width = frag_length) # Extend reads and compute coverage
cover_ext <- coverage(reads_ext)

# Coverage
# We split the genome up into 200 bp bins
bins <- tileGenome(seqinfo(reads), tilewidth = 200, cut.last.tile.in.chrom = TRUE)








# Understanding data tracks - a basic example
itrack <- IdeogramTrack(genome = "hg19", chromosome = "chrY")
gtrack <- GenomeAxisTrack()
set.seed(255)
lim <- c(starting, ending)
coords <- sort(c(lim[1], 
                 sample(seq(from = lim[1], to = lim[2]), 99), 
                 lim[2]))
dat <- runif(99, min = -10, max = 10)
dtrack <- DataTrack(data = dat, start = coords[-length(coords)],
                    end = coords[-1], chromosome = "chrY", genome = "hg19", 
                    name = "Uniform")
plotTracks(list(itrack, gtrack, dtrack), 
           from = lim[1], to = lim[2])
