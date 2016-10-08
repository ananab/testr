#############################################################################
# Sequence matching (SeqMatch) - test data
#############################################################################
# Purpose: Match amino acid sequences identified by mass spectrometry to all 
# amino acid sequences produced by humans.
# Packages: seqinr, dplyr, tidyr, ggplot2, "assertthat", "stringr"
lapply(c("dplyr", "seqinr", "tidyr", "ggplot2","assertthat","stringr"), library, character.only = T)
#############################################################################

# read csv file from test data
ori_lkup <- read.csv("data/testData.csv")

# select only relevant columns. the name "leading protein" is converted to "leading.protein" by R default
lkup <- ori_lkup %>% 
  select(peptide,leading.protein) %>% 
  # filter out proteins with no sequence information, doesn't matter for test but will become important for real data
  filter(peptide!="") %>% 
  
  # create new column indicating whether the N-terminus is blocked before tryptic digestion
  mutate(Ncapped = 0)

lkup$Ncapped[contains("n",FALSE,lkup$peptide)] <- 1

# variable that stores the peptide sequence in character
pepseq <- as.character(lkup$peptide)

# temp variable for formatting sequence
tmp <- sub(".?\\.","",pepseq)
tmp <- sub("n\\[.*\\]","",tmp)
tmp <- sub("\\..*","",tmp)
# add the formatted sequence as a new column named Lseq
lkup <- lkup %>% 
  mutate(Lseq=tmp) %>% 
  rename(Acc_id = leading.protein)

write.csv(lkup, "results/lkup.csv")

################################################################################
##### Load and clean the ref dataset.
orig_ref <- read.fasta("data/uniprot-all.fasta", seqtype = "AA", as.string = T)
ref <- data.frame(Acc_id = getName(orig_ref), 
                  Rseq = rapply(getSequence(orig_ref, as.string = T), c)) %>%
  ### Clean Acc_id.
  mutate(Acc_id = substr(as.character(Acc_id), 4, 9)) %>%
  ### Clean Rseq.
  mutate(Rseq = toupper(Rseq)) %>%
  ### Remove duplicates.
  distinct(Acc_id, Rseq)

### Save the clean version.
write.csv(ref,"results/ref.csv")

##### Merge the lkup and ref datasets by Acc_id. Do checks.
lkup <- read.csv("results/lkup.csv")
ref <- read.csv("results/ref.csv")
combo <- merge(lkup, ref, by = "Acc_id") %>% 
  ### Keep obs with unique values in Lseq and Rseq only.
  group_by(Rseq) %>%
  filter(!duplicated(Lseq)) 

combo$X.x = NULL
combo$X.y = NULL
### Check if all Acc_ids from lkup are in combo. Must be TRUE.
assert_that(length(intersect(combo$Acc_id, lkup$Acc_id)) == length(unique(lkup$Acc_id))) 

### Check if all Acc_ids from ref are in combo. Must be FALSE.
assert_that(length(intersect(combo$Acc_id, as.character(ref$Acc_id))) != length(unique(ref$Acc_id)))
### Set Match to 1 if Lseq is a substring of Rseq; 0 otherwise.
combo$Lseq <- as.character(combo$Lseq)
combo$Rseq <- as.character(combo$Rseq)

combo <- combo %>% 
  mutate(Match = mapply(grepl,Lseq,Rseq)) %>%
  ### Keep only matched obs.
  filter(Match == TRUE) %>% 
  ### Count character length of Lseq and Rseq.
  mutate(Llen = nchar(Lseq), 
         Rlen = nchar(Rseq),
         ### Get start and end indices of matched Lseq.
         Lstart = mapply(regexpr,Lseq,Rseq,fixed = TRUE),
         Lend = Lstart + Llen -1) %>% 
  ### Remove case with Lstart < 0.
  filter(Lstart>=0) %>% 
  ### Sort by Acc_id, Lstart and Lend.
  arrange(Acc_id,Lstart,Lend) 

### Check if all obs are matched. Must be TRUE.
assert_that(sum(combo$Match) == nrow(combo))
### See which obs are not matched.
which(combo$Match == FALSE)


### Check if each Acc_id is associated with a unique Rseq. Must be TRUE.
assert_that(length(unique(combo$Acc_id)) == length(unique(combo$Rseq)))

ACID_Rseq <- combo %>% 
  select(Acc_id,Rseq) %>% 
  distinct(.keep_all = TRUE)

accid <- sapply(as.numeric(rownames(ACID_Rseq)),
                function(i){
                  replicate(nchar(ACID_Rseq$Rseq[i]),
                            expr = ACID_Rseq$Acc_id[i])})

df <- data.frame(ACID <- c(unlist(accid)))

df <- df %>% 
  mutate(ACID = levels(unlist(accid))[ACID....c.unlist.accid..],
         index = rownames(df)) %>% 
  arrange(ACID,index) %>% 
  group_by(ACID) %>%
  mutate(index = row_number(),
         ACID....c.unlist.accid..=NULL,
         match = rep(nrow(df),0)) 
  
df$match[combo$Lstart] <- 1

# Accid AAindx match RpepIndx Group

