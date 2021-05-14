
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(DECIPHER)
library(stringr)
library(tidyverse)
library(R4RNA)

# Remove variable from memory
rm(list=ls())

########################################################################

# Script: Type III crRNA-guide design
  # Purpose: Find crRNA sequences in that are conserved in SARS-CoV-2 genomes
    # crRNAs should be able to target diverse SARS-CoV-2 genotypes
    # crRNAs should not cross-react with human mRNAs or other common respiratory flora

########################################################################


# Set working directory (path..to../type-III_diagnostic)
setwd("/type-III_diagnostic")

############### Run First Time Only ###############

msa <- readDNAStringSet("sub_algn.afa")

msa <- msa[width(msa) == max(width(msa))]

ref <- as.character(msa[grep("Wuhan",names(msa))])
ref <- str_split(ref, "")
table(ref[[1]] == "-")
ref <- as.character(ref[[1]])
ref <- data.frame(ref)
ref$algn <- seq(1,nrow(ref))
x <- !ref$ref == "-"
ref$pos[x] <- seq(1, 29903) 

ref <- ref %>% filter(!is.na(pos))

# Remove columns from alignment that aren't represented by a base in Ref
x <- ref$algn
start <- c(1, which(diff(x) != 1 & diff(x) != 0) + 1) # Find overlapping range of these continguous values
end <-  c(start - 1, length(x))

algn.trm <- narrow(msa, start = x[start[1]], end = x[end[2]])
p <- length(x[start])
for (i in seq(2, p)){
  
  algn.trm <- xscat(algn.trm, narrow(msa, start = x[start[i]], end = x[end[i+1]]))
  
}

names(algn.trm) <- names(msa)

#Import SARS-CoV-2 alignment
new <- readDNAStringSet("Alignment_w_new_seqs.fasta")

algn.trm <- algn.trm[1:2]
algn <- append(new, algn.trm)

writeXStringSet(algn, "Alignment_w_SARS_MERS.fasta")

############### Run from Here After First Time ###############

algn <- readDNAStringSet("Alignment_w_SARS_MERS.fasta")

sub.algn <- algn[c(7, 45652, 45651)]


###########################################################################
### Identify Guides

size.window <- 40      # Size of crRNA guide window to scan
five.prime.handle <- "AUUGCGAC"  # 5' handle of crRNA (derived from repeat)
mismatches.to.SARS.MERS <- 2    # (Min tolerated in target seed)
mismatches.in.handle <- 1       # (Min tolerated in crRNA 5' handle)
size.seed <- 18       # Size of seed region of cRNA (nearest to 5' end of crRNA handle)

###########################################################################

cross.reactive <- data.frame(matrix(nrow = 0, ncol = 5), stringsAsFactors = FALSE) ; colnames(cross.reactive) <- c("position", "sequence", "sars.seed.mm", "mers.seed.mm", "three.prime")

# Find guides with less than XX% conservation with SARS-CoV or MERS
for(i in seq(1, width(sub.algn[1]) - size.window  + 7) ){
  
   sars <- narrow(sub.algn[c(1,2)], start = i, end = size.window  -1 + i)
   sars.seed <- narrow(sars, start = width(sars)  - size.seed + 1, end = width(sars))
   sars <- narrow(sars, start = 1, end = width(sars) - size.seed)
   sars <- alignmentConservation(sars)
   sars.seed <- alignmentConservation(sars.seed)
   
   mers <- narrow(sub.algn[c(1,3)], start = i, end = size.window  -1 + i)
   mers.seed <- narrow(mers, start = width(mers)  - size.seed + 1 , end = width(mers))
   mers <- narrow(mers, start = 1, end = width(mers)- size.seed)
   mers <- alignmentConservation(mers)
   mers.seed <- alignmentConservation(mers.seed)
  
   seed.cutoff <- 1 - mismatches.to.SARS.MERS/size.seed
   
   if(sars.seed < seed.cutoff  & mers.seed < seed.cutoff) {
     
     x <- data.frame(matrix(nrow = 1, ncol = 5), stringsAsFactors = FALSE) ; colnames(x) <- colnames(cross.reactive)
     x$position <- i
     x$sequence <- as.character(narrow(sub.algn[1], start = i, end = size.window  -1 + i))
     x$sars.seed.mm <- round(size.seed - (sars.seed * size.seed), 0)
     x$mers.seed.mm <- round(size.seed - (mers.seed * size.seed), 0)
     x$three.prime <- as.character(narrow(sub.algn[1], start = size.window + i, end = size.window  + 7 + i))
     cross.reactive <- rbind(cross.reactive, x)
     
   }
}
  
# Find guides that aren't complementary to the 5' handle
keepers <- c()
for(i in seq(1, nrow(cross.reactive))) {
  
  g <- reverse(complement(DNAStringSet(RNAString(five.prime.handle))))
  
  x <- append( narrow(g, start = 1, end = 4), DNAStringSet(substr(cross.reactive$three.prime[i], 1, 4)))
  y <- append( narrow(g, start = 1, end = 4), DNAStringSet(substr(cross.reactive$three.prime[i], 5, 8)))
  
  if(alignmentConservation(x) == 0 & alignmentConservation(y) <= 0.5) {
    keepers <- append(keepers, i)
  }
  
}

cross.reactive.handle <- cross.reactive[keepers,]

## Cross reference guides with SARS-CoV-2 Genomes
# Find positions of variability in SARS-CoV-2 Genomes (run only first time)
this <- "No"

if(this == "Yes"){
ncov.positions <- data.frame(matrix(nrow = width (algn[7]), ncol = 6))
colnames(ncov.positions) <- c("A", "C", "G", "T", "-", "Other")
x <- algn[seq(1,length(algn)-2)]

for(i in seq(1, nrow(ncov.positions))) {

    
   ncov.positions$A[i] <- sum(letterFrequency(narrow(x, start = i, end = i), "A"))
   ncov.positions$C[i] <- sum(letterFrequency(narrow(x, start = i, end = i), "C"))
   ncov.positions$G[i] <- sum(letterFrequency(narrow(x, start = i, end = i), "G"))
   ncov.positions$T[i] <- sum(letterFrequency(narrow(x, start = i, end = i), "T"))
   ncov.positions$`-`[i] <- sum(letterFrequency(narrow(x, start = i, end = i), "-"))
   ncov.positions$Other[i] <- length(x) - sum( ncov.positions$A, ncov.positions$C, ncov.positions$G, ncov.positions$T, ncov.positions$`-` )
   
}

rm(x)

write.table(ncov.positions, file = "SARS-CoV-2_positional_variability.txt", sep = "\t", row.names = FALSE, quote = FALSE)}

if(this == "No") {ncov.positions <- read_delim("SARS-CoV-2_positional_variability.txt", delim = "\t")}



# Check guides for variability across their target positions

cross.reactive.handle$max.seed.variability <-0
cross.reactive.handle$avg.seed.variability <- 0
cross.reactive.handle$max.stem.variability <- 0
cross.reactive.handle$avg.stem.variability <- 0

for(i in seq(1, nrow(cross.reactive.handle))) {

x <- cross.reactive.handle$position[i]

y <- ncov.positions[seq(x, x + size.window - 1),1:4]
y$total <- (y$A + y$C + y$G + y$T)
y$percent <- round(as.numeric(apply(X=y[,1:4], MARGIN=1, FUN=max))/as.numeric(y$total)*100 , 2)

y.seed <- y[23:40,]
y <- y[1:22,]

cross.reactive.handle$max.seed.variability[i] <- 100 - min(y.seed$percent)
cross.reactive.handle$avg.seed.variability[i] <- round(100 - sum(y.seed$percent)/nrow(y.seed), 2)

cross.reactive.handle$max.stem.variability[i] <- 100 - min(y$percent)
cross.reactive.handle$avg.stem.variability[i] <- round(100 - sum(y$percent)/nrow(y), 2)

}


write.table(cross.reactive.handle, "Guides_found.txt", quote = FALSE, sep = "\t", row.names = FALSE)



###########################################################################
### Blast against common pathogens/microbial flora for cross reactivity

library(rBLAST)
#makeblastdb("Cross_reactivity/Cross_reactive.fna", dbtype = "nucl")   #Run first time (or if you add an organism)

###########################################################################

seq <- DNAStringSet(cross.reactive.handle$sequence)
names(seq) <- cross.reactive.handle$position

bl <- blast(db = "Cross_reactivity/Cross_reactive.fna", type = "blastn")
cl <- predict(bl, seq)

cross.ref <- data.frame(names(readDNAStringSet("Cross_reactivity/Cross_reactive.fna")), stringsAsFactors = FALSE)
colnames(cross.ref) <- "name"
cross.ref$SubjectID <- word(cross.ref$name)

cl$SubjectID <- as.character(cl$SubjectID)
cl <- cl %>% left_join(cross.ref, by = "SubjectID")

master <- cross.reactive.handle %>% filter(!position %in% cl$QueryID)



###########################################################################
### Add target annotation

ncov.features <- ape::read.gff("SARS-CoV-2.gff")

###########################################################################

genes <- ncov.features %>% filter(type == "gene")
genes$attributes <- word(genes$attributes, start = 3, end = 3, sep = "\\;")
genes <- genes %>% select(start, end, attributes)
genes$attributes <- str_replace(genes$attributes, "Name=", "")

target <- data.frame(matrix(ncol = 1, nrow = nrow(master)))
colnames(target) <- "target"

for (i in seq(1, nrow(master))){

  x <- which(master$position[i] >= genes$start & master$position[i] <= genes$end)
  if(isEmpty(x)) {target$target[i] <- "Non-coding"}
  if(length(x) >= 1) {target$target[i] <- genes$attributes[x]}

}

master <- cbind(master, target)

master$guide.sequence.5to3 <- as.character(reverse(complement(DNAStringSet(master$sequence))))

write.table(master, file = "Verified_guides.txt", quote = FALSE, sep = "\t", row.names = FALSE)


###########################################################################
### Visualize targets

#

###########################################################################

graphing <- master
graphing <- graphing[order(graphing$max.seed.variability, graphing$max.stem.variability),]
graphing <- graphing[1:15,]
graphing$number <- seq(1, nrow(graphing))

graphing$target <- factor(graphing$target, genes$attributes)

ggplot() +
  geom_segment(aes(x = 1, y = 0, xend = width(sub.algn[1]), yend = 0), color = "black") +
  geom_segment(aes(x = genes$start, y = 0, xend = genes$end, yend = 0, color = genes$attributes)) +
  geom_segment(aes(x = graphing$position, y = graphing$number, xend = graphing$position + size.window - 1, yend = graphing$number), color = "black", size = 2) +
  geom_segment(aes(x = graphing$position + size.window, y = graphing$number, xend = graphing$position + size.window + size.seed - 1, yend = graphing$number), color = "red", size = 2) +
  #geom_label(aes(label = genes$attributes, x = (genes$start + genes$end) / 2, y = 0, fill = genes$attributes), colour = "white", fontface = "bold", vjust = 0, nudge_y = 0.05) +
  geom_label(aes(label = graphing$position, x = (graphing$position + (graphing$position + size.window - 1)) / 2, y = graphing$number, fill = graphing$target), colour = "white", fontface = "bold", vjust = 0, nudge_y = 0.05) +
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())


###########################################################################

### Get more info on top targets

library()

###########################################################################

guide.info(graphing$position[2])


guide.info <- function(hit.position) {
  
  print(paste0(""))
  print(paste0("Info for guide at position ", hit.position, " targeting " , master$target[master$position == hit.position], " gene"))
  print(paste0(""))
  print(paste0("Target sequence [+ 3' sequence] = ", master$sequence[master$position == hit.position], " + ",
               master$three.prime[master$position == hit.position], "  " ))
  print(paste0("Guide  sequence [+ 5' handle]   = ", as.character(complement(DNAString(master$sequence[master$position == hit.position]))), " + ",
               as.character(reverse(five.prime.handle)), "  " ))
  print(paste0(""))
  print(paste0("Alignment to SARS-CoV:"))
  print(DNAMultipleAlignment(narrow(sub.algn[1:2], start = hit.position, end = hit.position + size.window -1)))
  print(paste0(""))
  print(paste0(""))
  print(paste0("Guide sequence (5' -- 3') = ", 
               as.character(reverse(complement(DNAString(master$sequence[master$position == hit.position])))), "  " ))
  print(paste0(""))
  BrowseSeqs(narrow(sub.algn[1:2], start = hit.position, end = hit.position + size.window -1))
}


pairwiseAlignment(narrow(sub.algn[1], start = hit.position, end = hit.position + size.window -1), narrow(sub.algn[2], start = hit.position, end = hit.position + size.window -1)) 




