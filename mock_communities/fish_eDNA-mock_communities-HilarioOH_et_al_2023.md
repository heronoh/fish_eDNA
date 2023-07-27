Evaluation of detection power of three 12S makers for 40 selected
species of the Jequitinhonha & São Francisco Rivers Basins
================
Hilário, OH; Mendes, IS; Sales, NG; Carvalho, DC
10/03/2022

------------------------------------------------------------------------

# Bioinformatics

## Data acquisiton

Download demultiplexed samples from *Base Space* using the *bs*
interface.

``` bash
#navigate to raw-data folder
cd $raw_data_folder/$run_folder;

#authenticate to BaseSpace (only at first log in)
bs auth;

#list datasets from runs on your BaseSpace
bs list datasets;

#create folders to organize fastq files
mkdir ~/runs/run_01mar21/fastq/;        #edna
mkdir ~/runs/run_09fev21/fastq/;        #edna
mkdir ~/runs/run_29jul20/fastq/;        #edna

#download runs from BaseSpace
bs download project -n fish_eDNA -o ~/runs/run_29jul20/fastq/ --extension=fastq.gz;          #primeira corrida LGC
bs download project -n eDNA_2run_B -o ~/runs/run_09fev21/fastq/ --extension=fastq.gz;        #segunda corrida LGC
bs download project -n iSeqRun2_Daniel -o ~/runs/run_01mar21/fastq/ --extension=fastq.gz;    #amostras iSeq ecomol

#organize all fastq files of each run in a single folder
mkdir ~/runs/run_01mar21/fastq/all;
mkdir ~/runs/run_09fev21/fastq/all;
mkdir ~/runs/run_29jul20/fastq/all;

#move all fastqfiles to a single folder
mv ~/runs/run_01mar21/fastq/*/*fastq.gz ~/runs/run_01mar21/fastq/all;
mv ~/runs/run_09fev21/fastq/*/*fastq.gz ~/runs/run_09fev21/fastq/all;
mv ~/runs/run_29jul20/fastq/*/*fastq.gz ~/runs/run_29jul20/fastq/all;
```

<br>

## Load *R libs* and system programs

``` r
# 0 - load libraries ----
{
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggbreak)
  library(ggtree)
  library(phyloseq)
  library(Biostrings)
  library(Matrix)
  library(ShortRead)
  library(dada2)
  library(DECIPHER)
  library(future)
  library(vegan)
  library(ape)
  library(phangorn)
  library(adegenet)
}

#set complete path to cutadapt executable
cutadapt <- "/usr/local/bin/cutadapt"

#important 
prjct_path <- "~/prjcts/fish_eDNA/sfjq"

notes_path <- paste0(prjct_path,"/notes")

results_path <- paste0(prjct_path,"/results")

figs_path <- paste0(results_path,"/figs")

prcj_radical <- "SFJq_fish_metabarcoding"
#path to project data folder were the processed reads will be stored
data_path <- paste0(prjct_path,"/data/reads")
```

<br>

## Quality control

Chech overall quality of sequencing runs for all samples

``` bash
```

<br>

## Demultiplex SFJQ sample (MiFish & NeoFish mixed)

``` bash
#Demultiplex SFJQ sample (MiFish & NeoFish mixed)
#samples MiniSeq LGC
cutadapt -j 79 --no-indels  -g file:/home/heron/prjcts/fish_eDNA/sfjq/data/primers_neo_mi.fasta  -G file:/home/heron/prjcts/fish_eDNA/sfjq/data/primers_neo_mi.fasta  -o /home/heron/runs/run_09fev21/fastq/all/sfjq_dmx/SFJQ-{name1}-{name2}_R1_001.fastq.gz  -p /home/heron/runs/run_09fev21/fastq/all/sfjq_dmx/SFJQ-{name1}-{name2}_R2_001.fastq.gz  /home/heron/runs/run_09fev21/fastq/all/SFJQ-neo-mi_S23_L001_R1_001.fastq  /home/heron/runs/run_09fev21/fastq/all/SFJQ-neo-mi_S23_L001_R2_001.fastq 2> /home/heron/runs/run_09fev21/fastq/all/sfjq_dmx/cut_SFJQ_demux.txt

cp ~/runs/run_09fev21/fastq/all/sfjq_dmx/SFJQ-neo_FWD-neo_REV.* ~/runs/run_01mar21/fastq/all/
cp ~/runs/run_09fev21/fastq/all/sfjq_dmx/SFJQ-mif_FWD-mif_REV.* ~/runs/run_09fev21/fastq/all/

mv ~/runs/run_01mar21/fastq/all/SFJQ-neo-mi_S23_L001_R* ~/runs/run_01mar21/fastq/all/sfjq_dmx/


#samples iSeq Ecomol
cutadapt -j 79 --no-indels  -g file:/home/heron/prjcts/fish_eDNA/sfjq/data/primers_neo_mi.fasta  -G file:/home/heron/prjcts/fish_eDNA/sfjq/data/primers_neo_mi.fasta  -o /home/heron/runs/run_01mar21/fastq/all/sfjq_dmx/Da23-{name1}-{name2}_R1_001.fastq.gz  -p /home/heron/runs/run_01mar21/fastq/all/sfjq_dmx/Da23-{name1}-{name2}_R2_001.fastq.gz  /home/heron/runs/run_01mar21/fastq/all/Da23_S72_L001_R1_001.fastq /home/heron/runs/run_01mar21/fastq/all/Da23_S72_L001_R2_001.fastq 2> /home/heron/runs/run_01mar21/fastq/all/sfjq_dmx/cut_SFJQ_demux.txt

cp ~/runs/run_01mar21/fastq/all/sfjq_dmx/Da23-neo_FWD-neo_REV.* ~/runs/run_01mar21/fastq/all/
cp ~/runs/run_01mar21/fastq/all/sfjq_dmx/Da23-mif_FWD-mif_REV.* ~/runs/run_01mar21/fastq/all/

mv ~/runs/run_01mar21/fastq/all/Da23_S72_L001_R* ~/runs/run_01mar21/fastq/all/sfjq_dmx/
```

## Set path to raw data

``` r
#1 - load runs raw data ----
## All libs are demultiplexed
{
  # PATH to the directory containing raw fastq files after unzipping.
  libs_path1 <- "~/runs/run_29jul20/fastq/all"
  libs_path2 <- "~/runs/run_09fev21/fastq/all" 
  libs_path3 <- "~/runs/run_01mar21/fastq/all" 
}
#check content
list.files(path = libs_path1,pattern = "fastq") 
list.files(path = libs_path2,pattern = "fastq") 
list.files(path = libs_path3,pattern = "fastq") 
```

### Identify sample names radicals

``` r
#2 - get sample names ----

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
{
  all_fnFs1 <- sort(list.files(libs_path1, pattern="_R1_001.fastq", full.names = TRUE))
  all_fnRs1 <- sort(list.files(libs_path1, pattern="_R2_001.fastq", full.names = TRUE))

  all_fnFs2 <- sort(list.files(libs_path2, pattern="_R1_001.fastq", full.names = TRUE))
  all_fnRs2 <- sort(list.files(libs_path2, pattern="_R2_001.fastq", full.names = TRUE))

  all_fnFs3 <- sort(list.files(libs_path3, pattern="_R1_001.fastq", full.names = TRUE))
  all_fnRs3 <- sort(list.files(libs_path3, pattern="_R2_001.fastq", full.names = TRUE))

  all_fnFs <- c(all_fnFs1,all_fnFs2,all_fnFs3)
  all_fnRs <- c(all_fnRs1,all_fnRs2,all_fnRs3)
}
#load csv with primers and respective samples
primers_n_samples <- read.csv(file = "~/prjcts/fish_eDNA/sfjq/data/primers_n_samples_sfjq.csv",
         header = TRUE)


#3 - map sample names to reads files ----
primers_n_samples <- primers_n_samples %>%
  mutate("R1" = "R1",
         "R2" = "R2")

for (sample in 1:nrow(primers_n_samples)) {
  
  primers_n_samples$R1[sample] <-
   all_fnFs[grep(pattern =  paste0("/",primers_n_samples$File_name[sample]),x = all_fnFs)]
  
  primers_n_samples$R2[sample] <-
   all_fnRs[grep(pattern =  paste0("/",primers_n_samples$File_name[sample]),x = all_fnRs)]

}
```

<br>

### Define sample levels

``` r
#4 - set sample levels
primers_n_samples$File_name

sample_levels <- c(
"Da23-mif", "SFJQ-mif",
"Da23-neo", "SFJQ-neo",
"Da20","SFnNorm-mi",
"Da19","SFnNorm-neo",
"Da22","SFNorm-mi",
"Da21","SFNorm-neo",
"pJequei-N-norm-N","pJequei-N-norm-M","pJequei-N-norm-T",
"pJequei-norm-N","pJequei-norm-M","pJequei-norm-T",
"Cassaum","neg-PCR2")

primer_levels <- c("NeoFish", "MiFish", "Teleo","NeoFish/MiFish", "NeoFish/MiFish/Teleo")
```

<br> \#\#\# Remove primers from reads

As the primer-derived sequences are identical, they are not informative
and thus must be removed for the following steps.

<br>

#### Load primer sequences

``` r
#4- identify primers ----

#primers sequences used for each sample
# inosine pairs with A, C, U
#                    T, G, A = IUPAC code:  D
#cutadapt  accepts IUPAC code !!!!!!!!
{
  #neo
  neo_FWD <- "CGCCGTCGCAAGCTTACCCT"
  names(neo_FWD) <- "neo_FWD"
  neo_REV <- "AGTGACGGGCGGTGTGTGC"
  names(neo_REV) <- "neo_REV" 
  
  #mif
  mif_FWD <- "GTCGGTAAAACTCGTGCCAGC"
  names(mif_FWD) <- "mif_FWD"
  mif_REV <- "ACATAGTGGGGTATCTAATCCCAGTTTG"
 # mif_REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" #original
  names(mif_REV) <- "mif_REV" 
  
  #tel
  tel_FWD <- "ACACCGCCCGTCACTCT"
  names(tel_FWD) <- "tel_FWD"
  tel_REV <- "ACTTCCGGTACACTTACCATG"
  names(tel_REV) <- "tel_REV"  
  
  
#creates a list of single row tibbles for each primer
primers <- tibble(Primers = c(neo_FWD,neo_REV,
                              mif_FWD,mif_REV,
                              tel_FWD,tel_REV)) %>% 
  mutate(`Primer names`= names(Primers)) %>% 
  split(1:nrow(.))
}
```

<br>

#### Generate sequences for complement, reverse, and reverse complement for each primer

The function *allOrients* is used to generate all possible orientations
for primers FWD e REV.

``` r
#5 - check primer orientation ----

#function to get all possible primer orientations
allOrients <- function(primers) {
   # Create all orientations of the input sequence
    # Must be a tibble with cols = c(Primers,`Primer names`)
  
   require(Biostrings)
   dna <- Biostrings::DNAString(primers$Primers)  # The Biostrings works w/ DNAString objects rather than character vectors
   orients <- c(Forward = dna, 
                Complement = Biostrings::complement(dna), 
                Reverse = Biostrings::reverse(dna),
                RevComp = Biostrings::reverseComplement(dna))
   names(orients) <- paste0(names(orients))
   
   primer_tbl <- sapply(orients, toString)
   
   primer_tbl <- dplyr::tibble(Sequence = primer_tbl,
                        `Primer orientation` = names(primer_tbl)) %>% 
     dplyr::mutate(Primer = primers$`Primer names`) %>%
     unite(col=`Orientation name`, Primer ,`Primer orientation`,remove = FALSE)
   
   return(primer_tbl)  # Convert back to character vector
}


#apply function 
primers_all_orients <- purrr::map_dfr(primers, allOrients)

names(primers_all_orients$Sequence) <- primers_all_orients$`Orientation name`
```

<br>

#### Remove reads with undetermined bases **(Ns)** and unpaired

Reads with undetermined bases prevent proper primer identification and
ASV determination. These sequences must be removed from the data.

``` r
#6 - pre filter reads with Ns for primer checking ----
# create names for N-cleaned files

primers_n_samples <- primers_n_samples %>%
  mutate("R1 N-cleaned" = "R1 N-cleaned",
         "R2 N-cleaned" = "R2 N-cleaned")

for (sample in 1:nrow(primers_n_samples)) {
  primers_n_samples$`R1 N-cleaned`[sample] <-
   paste0(data_path,"/N-cleaned/",primers_n_samples$File_name[sample],"_R1_N-cleaned.fastq.gz")
  primers_n_samples$`R2 N-cleaned`[sample] <-
   paste0(data_path,"/N-cleaned/",primers_n_samples$File_name[sample],"_R2_N-cleaned.fastq.gz")
}

# remove reads with Ns to make primer filtering more accurate

dada2::filterAndTrim(
  fwd = primers_n_samples$R1, filt = primers_n_samples$`R1 N-cleaned`, 
  rev = primers_n_samples$R2, filt.rev = primers_n_samples$`R2 N-cleaned`,
  maxN = 0, multithread = TRUE, matchIDs = TRUE,
  verbose = TRUE, compress = TRUE)

# pivote table to longer format
primers_n_samples <- primers_n_samples %>% 
  pivot_longer(cols = c(R1,R2,`R1 N-cleaned`,`R2 N-cleaned`),
               names_to = "Stage", values_to = "Read file")
```

<br>

#### Count primer presence on reads

Before primer removal it is possible to count their presence on the
reads. This procedures is carried on independently for each sample.

``` r
#6 - count primer orientation hits ----

#function to count primer on each specific library
primerHits <- function(primer, fn) {
   # Counts number of reads in which the primer is found
   nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
   return(sum(nhits > 0))
}

#function to call primerHits for multiple primers
multi_primerHits <- function(Read_file,primers){
  primer_counts <- purrr::map_df(primers,.f = primerHits, fn = Read_file)
  primer_counts <- primer_counts %>%  mutate(`Read file` = Read_file)
  return(primer_counts)
}
###########

#vector of read files to look on for primers
reads_seqs <- primers_n_samples %>% 
  filter(Stage %in% c("R1 N-cleaned", "R2 N-cleaned")) %>% 
  select(`Read file`) %>% as.list()
 

#named vector of primer sequences
primers_seqs <- primers_all_orients$Sequence


cores_to_be_used <- future::availableCores() - 2 # Usar todos os cores -2 = 78

future::plan(future::multisession(workers = cores_to_be_used))


#count primers
primers_in_Nreads <- furrr::future_map_dfr(reads_seqs$`Read file`, .f = multi_primerHits, primers = primers_seqs, .options = furrr::furrr_options(seed = NULL))

#get sample information into primers_in_Nreads table
primers_in_Nreads <- left_join(primers_in_Nreads,primers_n_samples,by = "Read file")

# 
# primers_in_Nreads_bckp <- primers_in_Nreads
# primers_in_Nreads <- primers_in_Nreads_bckp
```

#### Prepare primer counts for ploting

``` r
#7- prepare primer counts for plots ----

# cat(paste0(colnames(primers_in_Nreads),"\n"))

primers_in_Nreads <-
  primers_in_Nreads %>% 
  select(# `Read file
 File_name, Type, Group, Library, Primer, Run, Stage,
         neo_FWD_Forward, neo_REV_Forward, neo_FWD_Complement, neo_REV_Complement, 
         neo_FWD_Reverse, neo_REV_Reverse, neo_FWD_RevComp, neo_REV_RevComp, 
         mif_FWD_Forward, mif_REV_Forward, mif_FWD_Complement, mif_REV_Complement, 
         mif_FWD_Reverse, mif_REV_Reverse, mif_FWD_RevComp, mif_REV_RevComp, 
         tel_FWD_Forward, tel_REV_Forward, tel_FWD_Complement, tel_REV_Complement, 
         tel_FWD_Reverse, tel_REV_Reverse, tel_FWD_RevComp, tel_REV_RevComp)


#write.csv(x = primer_hits_tbl, file = "/home/heron/prjcts/fish_eDNA/notes/jequiDNApool/csv/primers_hits_in_reads.csv")

str(primers_in_Nreads)
colnames(primers_in_Nreads)
rownames(primers_in_Nreads)

primers_in_Nreads$Primer
primers_in_Nreads$Library


#8- prepare primer counts for plots in ggplot----

#convert primer hits table to long format
primers_in_Nreads_long <- primers_in_Nreads %>% 
  gather(key = Sequences, 
         value = Count,  
         neo_FWD_Forward, neo_FWD_Complement, neo_FWD_Reverse,neo_FWD_RevComp, 
         neo_REV_Forward, neo_REV_Complement, neo_REV_Reverse, neo_REV_RevComp,
         mif_FWD_Forward, mif_FWD_Complement, mif_FWD_Reverse, mif_FWD_RevComp, 
         mif_REV_Forward, mif_REV_Complement, mif_REV_Reverse, mif_REV_RevComp,
         tel_FWD_Forward, tel_FWD_Complement, tel_FWD_Reverse, tel_FWD_RevComp, 
         tel_REV_Forward, tel_REV_Complement, tel_REV_Reverse, tel_REV_RevComp
         ) %>% 
  mutate(Sequences = factor(Sequences,
                            levels = c("neo_FWD_Forward","neo_FWD_RevComp",
                                       "neo_REV_Forward","neo_REV_RevComp",
                                       "neo_FWD_Complement","neo_FWD_Reverse",
                                       "neo_REV_Complement","neo_REV_Reverse",
                                       
                                       
                                       "mif_FWD_Forward","mif_FWD_RevComp",
                                       "mif_REV_Forward","mif_REV_RevComp",
                                       "mif_FWD_Complement","mif_FWD_Reverse",
                                       "mif_REV_Complement","mif_REV_Reverse",
                                       
                                       
                                       "tel_FWD_Forward","tel_FWD_RevComp",
                                       "tel_REV_Forward","tel_REV_RevComp",
                                       "tel_FWD_Complement","tel_FWD_Reverse",
                                       "tel_REV_Complement","tel_REV_Reverse")),
                                       
                                       
         File_name = factor(File_name,levels = sample_levels),
         Run = as.factor(Run),
         Primer = factor(Primer,levels = c("NeoFish",
                                           "MiFish",
                                           "Teleo",
                                           "NeoFish/MiFish",
                                           "NeoFish/MiFish/Teleo"))) 



# PLOT 1: primers counts in reads tile plot - only primers FWD & REV, foward & revcomp ----
primers_tile <- 
  primers_in_Nreads_long %>% 
  # filter(Sequences  %in% c(
  #   "mif_REV_RevComp", "mif_REV_Forward", "mif_FWD_RevComp", "mif_FWD_Forward",
  #   "neo_REV_RevComp", "neo_REV_Forward", "neo_FWD_RevComp", "neo_FWD_Forward",
  #   "tel_REV_RevComp", "tel_REV_Forward", "tel_FWD_RevComp", "tel_FWD_Forward")) %>% 
  mutate(File_name = factor(File_name,levels = sample_levels)) %>% 
  ggplot2::ggplot(aes(y=File_name,x=Sequences,fill=log10(Count)
                      # ,col=Stage
                      )) +
  geom_tile()+
  geom_text(aes(label = Count),size=1)+
  # scale_fill_gradient(low="white", high="darkgreen",trans="log10") +
  scale_fill_gradientn(name = "Primer counts",
                       colours = c("white","darkgreen"),
                       values = c(0,1),
                       na.value ="white") +
  theme_light(base_line_size = 1,base_size = 6) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  geom_hline(yintercept = c(40.5,82.5,86.5,116.5),color = "grey") +
  geom_vline(xintercept = c(4.5,8.5,12.5,16.5),color = "grey") +
  # coord_fixed(ratio = 0.20) +
  xlab("Primers") +
  ylab("Amostra") +
  ggtitle(label = "eDNA 1st, 2nd & 3rd runs",
              subtitle = "Primer presence on sample reads") 
# +
#   facet_wrap(~Run, drop = TRUE)
  # facet_wrap(~Stage, drop = TRUE)
# +
#   facet_wrap(~Primer)
  # annotate(geom = "rect", 
  #          xmin=c(0.5), xmax=c(4.5),  #neo
  #          ymin=0.5, ymax=116.5, 
  #          alpha=0.05,fill="yellow") +
  # annotate(geom = "rect", 
  #          xmin=c(4.5), xmax=c(8.5),  #mif
  #          ymin=c(0,116.5), ymax=c(86.5,158.5), 
  #          alpha=0.05,fill="yellow")+
  # annotate(geom = "rect", 
  #          xmin=c(8.5), xmax=c(12.5),  tel
  #          ymin=c(0,86.5), ymax=c(82.5,158.5), 
  #          alpha=0.05,fill="yellow")

primers_tile

ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/1-primers_in_reads_all_FR.png",
     plot = primers_tile,
     device = "png",
     width = 27,
     height = 40,
     units = "cm",
     dpi = 600)

ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/1-primers_in_reads_all_FR.pdf",
     plot = primers_tile,
     device = "pdf",
     width = 27,
     height = 40,
     units = "cm",
     dpi = 600)


#9- write csv file with primer hits counts per lib ----
write.csv(x = primer_hits_tbl,file = "~/prjcts/fish_eDNA/sfjq/results/primers_hits_tbl.csv",
          row.names = FALSE)
```

### Remove primers from reads

#### Primer removal with ***Cutadapt***

The ***cutadapt*** software
([DOI:10.14806/ej.17.1.200](http://journal.embnet.org/index.php/embnetjournal/article/view/200))
was used for primer removal on read sequences.

``` r
#10 - cutadapt ----

#set or create cutadapt processed reads dir path
path.cut <- file.path(data_path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
```

<br>

#### Generate and execute primer-specific commands

The original DADA2 ITS protocol removes only *FWD* and *REV reverse
complement* sequences. This protocol is adjusted for selecting reads
only of the expected primer and removing the primer.

``` r
# opitional: remove all primers from all reads and samples ----

#10 - map sample names to reads files ----

#name outputs
cutadapt_files <- primers_n_samples %>% 
  filter(Stage %in% c("R1 N-cleaned","R2 N-cleaned")) %>% 
  mutate(`Read file` = str_replace_all(.$`Read file`,pattern = "N-cleaned",replacement = "cutadapt")) %>% 
  mutate(Stage = str_replace_all(.$Stage,pattern = "N-cleaned",replacement = "cutadapt"))


primers_n_samples <- bind_rows(primers_n_samples,cutadapt_files)

#all ----
{ 
          #make reverse complements
        #  the XXX_Complement and XXX_reverse have no hits so were ignored at last plot and from now on
          all_FWD.orients <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("FWD_Forward"))] 
          all_FWD.RC <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("FWD_RevComp"))] 
          all_REV.orients <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("REV_Forward"))] 
          all_REV.RC <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("REV_RevComp"))] 
  
}
  
#remove primers and filter only the reads that contain the expected primer ----
{
  #MiFish ----
  mif_FWD.orients <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("mif_FWD_Forward"))] 
  mif_FWD.RC <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("mif_FWD_RevComp"))] 
  mif_REV.orients <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("mif_REV_Forward"))] 
  mif_REV.RC <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("mif_REV_RevComp"))] 
  
  #NeoFish ----
  neo_FWD.orients <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("neo_FWD_Forward"))] 
  neo_FWD.RC <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("neo_FWD_RevComp"))] 
  neo_REV.orients <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("neo_REV_Forward"))] 
  neo_REV.RC <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("neo_REV_RevComp"))] 
  
  #Teleo ----
  tel_FWD.orients <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("tel_FWD_Forward"))]
  tel_FWD.RC <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("tel_FWD_RevComp"))]
  tel_REV.orients <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("tel_REV_Forward"))]
  tel_REV.RC <-  primers_all_orients$Sequence[primers_all_orients$`Orientation name` %>% grep(pattern = c("tel_REV_RevComp"))]
  }
  
                #creat flags
                # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
                all_R1.flags <- paste("-g", all_FWD.orients, "-a", all_REV.RC)
                # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
                all_R2.flags <- paste("-G", all_REV.orients, "-A", all_FWD.RC)
                
#create primer-specific tags ----
  {
  #MiFish ----
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  mif_R1.flags <- paste("-g", mif_FWD.orients, "-a", mif_REV.RC)
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  mif_R2.flags <- paste("-G", mif_REV.orients, "-A", mif_FWD.RC)
  
  #NeoFish ----
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  neo_R1.flags <- paste("-g", neo_FWD.orients, "-a", neo_REV.RC)
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  neo_R2.flags <- paste("-G", neo_REV.orients, "-A", neo_FWD.RC)
  
  #Teleo ----
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  tel_R1.flags <- paste("-g", tel_FWD.orients, "-a", tel_REV.RC)
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  tel_R2.flags <- paste("-G", tel_REV.orients, "-A", tel_FWD.RC)
  }

#cutadapt files path and names ----
{
        all_fnFs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 cutadapt"]
        names(all_fnFs.cut) <- all_fnFs.cut %>% 
            str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
            str_remove(pattern = "_cutadapt.fastq.gz")
        all_fnRs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 cutadapt"]
        names(all_fnRs.cut) <- all_fnRs.cut %>% 
            str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
            str_remove(pattern = "_cutadapt.fastq.gz")
        
        all_fnFs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 N-cleaned"]
        all_fnRs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 N-cleaned"]

}

{
  #neo ----
  neo_fnFs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 cutadapt" & primers_n_samples$Primer == "NeoFish"]
    names(neo_fnFs.cut) <- neo_fnFs.cut %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
      str_remove(pattern = "_cutadapt.fastq.gz")

  neo_fnRs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 cutadapt" & primers_n_samples$Primer == "NeoFish"]
    names(neo_fnRs.cut) <- neo_fnRs.cut %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
      str_remove(pattern = "_cutadapt.fastq.gz")

  neo_fnFs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 N-cleaned" & primers_n_samples$Primer == "NeoFish"]
    names(neo_fnFs.filtN) <- neo_fnFs.filtN %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/N-cleaned/")) %>% 
      str_remove(pattern = "_N-cleaned.fastq.gz")
  
  neo_fnRs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 N-cleaned" & primers_n_samples$Primer == "NeoFish"]
    names(neo_fnRs.filtN) <- neo_fnRs.filtN %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/N-cleaned/")) %>% 
      str_remove(pattern = "_N-cleaned.fastq.gz")
    
  #mif ----
  mif_fnFs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 cutadapt" & primers_n_samples$Primer == "MiFish"]
    names(mif_fnFs.cut) <- mif_fnFs.cut %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
      str_remove(pattern = "_cutadapt.fastq.gz")

  mif_fnRs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 cutadapt" & primers_n_samples$Primer == "MiFish"]
    names(mif_fnRs.cut) <- mif_fnRs.cut %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
      str_remove(pattern = "_cutadapt.fastq.gz")

  mif_fnFs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 N-cleaned" & primers_n_samples$Primer == "MiFish"]
    names(mif_fnFs.filtN) <- mif_fnFs.filtN %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/N-cleaned/")) %>% 
      str_remove(pattern = "_N-cleaned.fastq.gz")
  
  mif_fnRs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 N-cleaned" & primers_n_samples$Primer == "MiFish"]
    names(mif_fnRs.filtN) <- mif_fnRs.filtN %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/N-cleaned/")) %>% 
      str_remove(pattern = "_N-cleaned.fastq.gz")
    
      
  #tel ----
  tel_fnFs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 cutadapt" & primers_n_samples$Primer == "Teleo"]
    names(tel_fnFs.cut) <- tel_fnFs.cut %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
      str_remove(pattern = "_cutadapt.fastq.gz")

  tel_fnRs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 cutadapt" & primers_n_samples$Primer == "Teleo"]
    names(tel_fnRs.cut) <- tel_fnRs.cut %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
      str_remove(pattern = "_cutadapt.fastq.gz")

  tel_fnFs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 N-cleaned" & primers_n_samples$Primer == "Teleo"]
    names(tel_fnFs.filtN) <- tel_fnFs.filtN %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/N-cleaned/")) %>% 
      str_remove(pattern = "_N-cleaned.fastq.gz")
  
  tel_fnRs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 N-cleaned" & primers_n_samples$Primer == "Teleo"]
    names(tel_fnRs.filtN) <- tel_fnRs.filtN %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/N-cleaned/")) %>% 
      str_remove(pattern = "_N-cleaned.fastq.gz")
    
          
  #neo/mif & neo/mif/tel ----
  nmt_fnFs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 cutadapt" & primers_n_samples$Primer %in% c("NeoFish/MiFish","NeoFish/MiFish/Teleo")]
    names(nmt_fnFs.cut) <- nmt_fnFs.cut %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
      str_remove(pattern = "_cutadapt.fastq.gz")

  nmt_fnRs.cut <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 cutadapt" & primers_n_samples$Primer %in% c("NeoFish/MiFish","NeoFish/MiFish/Teleo")]
    names(nmt_fnRs.cut) <- nmt_fnRs.cut %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/cutadapt/")) %>% 
      str_remove(pattern = "_cutadapt.fastq.gz")

  nmt_fnFs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 N-cleaned" & primers_n_samples$Primer %in% c("NeoFish/MiFish","NeoFish/MiFish/Teleo")]
    names(nmt_fnFs.filtN) <- nmt_fnFs.filtN %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/N-cleaned/")) %>% 
      str_remove(pattern = "_N-cleaned.fastq.gz")
  
  nmt_fnRs.filtN <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 N-cleaned" & primers_n_samples$Primer %in% c("NeoFish/MiFish","NeoFish/MiFish/Teleo")]
    names(nmt_fnRs.filtN) <- nmt_fnRs.filtN %>% 
      str_remove(pattern = paste0("/home/heron/prjcts/fish_eDNA/sfjq/data/reads/N-cleaned/")) %>% 
      str_remove(pattern = "_N-cleaned.fastq.gz")
  
}


# e essa ordem tá perigosa
names(primers_n_samples$`Read file`)<- primers_n_samples$File_name


primers_n_samples$`Read file`[primers_n_samples$Stage %in% c("R1 N-cleaned", "R2 N-cleaned")] %>%  length()

#TODO esse tem q ir pro purrr...


# Run Cutadapt
#   output folder must exist
for(i in 1:40) {
# for(i in seq_along(all_fnFs)) {
system2(cutadapt, args = c(all_R1.flags, all_R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
"-o", all_fnFs.cut[i], "-p", all_fnRs.cut[i], # output files
all_fnFs.filtN[i], all_fnRs.filtN[i],  # input files
"--minimum-length 10")) # guarantee no zerolength reads
}


length(neo_fnFs.cut)
length(mif_fnFs.cut)
length(mif_fnRs.cut)
length(tel_fnFs.cut)
length(nmt_fnFs.cut)


#run cutadapt by primer ----

#MiFish ----
for(i in 1:length(mif_fnFs.cut)) {
# for(i in seq_along(all_fnFs)) {
system2(cutadapt, args = c(mif_R1.flags, mif_R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
"-o", mif_fnFs.cut[i], "-p", mif_fnRs.cut[i], # output files
mif_fnFs.filtN[i], mif_fnRs.filtN[i],  # input files
"--minimum-length 10 --discard-untrimmed")) # guarantee no zerolength reads
}

#NeoFish ----
for(i in 1:length(neo_fnFs.cut)) {
# for(i in seq_along(all_fnFs)) {
system2(cutadapt, args = c(neo_R1.flags, neo_R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
"-o", neo_fnFs.cut[i], "-p", neo_fnRs.cut[i], # output files
neo_fnFs.filtN[i], neo_fnRs.filtN[i],  # input files
"--minimum-length 10 --discard-untrimmed")) # guarantee no zerolength reads
}


#Teleo ----
for(i in 1:length(tel_fnFs.cut)) {
# for(i in seq_along(all_fnFs)) {
system2(cutadapt, args = c(tel_R1.flags, tel_R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
"-o", tel_fnFs.cut[i], "-p", tel_fnRs.cut[i], # output files
tel_fnFs.filtN[i], tel_fnRs.filtN[i],  # input files
"--minimum-length 10 --discard-untrimmed")) # guarantee no zerolength reads
}


#neo/mif & neo/mif/tel ----
  
for(i in 1:length(nmt_fnFs.cut)) {
# for(i in seq_along(all_fnFs)) {
system2(cutadapt, args = c(all_R1.flags, all_R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
"-o", nmt_fnFs.cut[i], "-p", nmt_fnRs.cut[i], # output files
nmt_fnFs.filtN[i], nmt_fnRs.filtN[i],  # input files
"--minimum-length 10 --discard-untrimmed")) # guarantee no zerolength reads
}


# as  Da23-mif & neo & SFJQ-mif & neo were previously demultiplexed, they do not bear the primers, so all reads are filtered out


# have no time to implement a solution, so:
#        cp ../N-cleaned/Da23-mif_R1_N-cleaned.fastq.gz Da23-mif_R1_cutadapt.fastq.gz
#        cp ../N-cleaned/Da23-mif_R2_N-cleaned.fastq.gz Da23-mif_R2_cutadapt.fastq.gz
#        cp ../N-cleaned/Da23-neo_R2_N-cleaned.fastq.gz Da23-neo_R2_cutadapt.fastq.gz
#        cp ../N-cleaned/Da23-neo_R1_N-cleaned.fastq.gz Da23-neo_R1_cutadapt.fastq.gz
#        cp ../N-cleaned/SFJQ-mif_R1_N-cleaned.fastq.gz SFJQ-mif_R1_cutadapt.fastq.gz
#        cp ../N-cleaned/SFJQ-mif_R2_N-cleaned.fastq.gz SFJQ-mif_R2_cutadapt.fastq.gz
#        cp ../N-cleaned/SFJQ-neo_R1_N-cleaned.fastq.gz SFJQ-neo_R1_cutadapt.fastq.gz
#        cp ../N-cleaned/SFJQ-neo_R2_N-cleaned.fastq.gz SFJQ-neo_R2_cutadapt.fastq.gz
```

<br>

#### Check primer removal

Change de library number index in order to check the presence of
remaining primer sequences on each lib data. It is expected that all
removed orientation counts change to zero since the primer sequences are
removed.

``` r
# only ofr didatic purposes (or logical error checking)
# 8 - check for remaining adapters ----
# As a sanity check, we will count the presence of primers in the first cutadapt-ed sample:


#vector of read files to look on for primers
reads_seqs_cut <- primers_n_samples %>% 
  filter(Stage %in% c("R2 cutadapt", "R2 cutadapt")) %>% 
  select(`Read file`) %>% as.list()


#count primers
future::plan(future::multisession(workers = cores_to_be_used))
primers_in_cut_reads <- furrr::future_map_dfr(reads_seqs_cut$`Read file`, .f = multi_primerHits, primers = primers_seqs, .options = furrr::furrr_options(seed = NULL))

 # 
# primers_in_Nreads <- purrr::map_df(reads_seqs,.f = multi_primerHits, primers = primers_seqs)

#get sample information into primers_in_Nreads table
primers_in_cut_reads <- left_join(primers_in_cut_reads,primers_n_samples,by = "Read file")

# primers_in_cut_reads_bckp <- primers_in_cut_reads
# primers_in_cut_reads <- primers_in_cut_reads_bckp

primers_in_cut_reads <-
  primers_in_cut_reads %>% 
  select(# `Read file
 File_name, Type, Group, Library, Primer, Run, Stage,
         neo_FWD_Forward, neo_REV_Forward, neo_FWD_Complement, neo_REV_Complement, 
         neo_FWD_Reverse, neo_REV_Reverse, neo_FWD_RevComp, neo_REV_RevComp, 
         mif_FWD_Forward, mif_REV_Forward, mif_FWD_Complement, mif_REV_Complement, 
         mif_FWD_Reverse, mif_REV_Reverse, mif_FWD_RevComp, mif_REV_RevComp, 
         tel_FWD_Forward, tel_REV_Forward, tel_FWD_Complement, tel_REV_Complement, 
         tel_FWD_Reverse, tel_REV_Reverse, tel_FWD_RevComp, tel_REV_RevComp)

#convert primer hits table to long format
primers_in_cut_reads_long <- primers_in_cut_reads %>% 
  gather(key = Sequences, 
         value = Count,  
         neo_FWD_Forward, neo_FWD_Complement,
         neo_FWD_Reverse,neo_FWD_RevComp, 
         neo_REV_Forward, neo_REV_Complement,
         neo_REV_Reverse, neo_REV_RevComp,
         mif_FWD_Forward, mif_FWD_Complement,
         mif_FWD_Reverse, mif_FWD_RevComp, 
         mif_REV_Forward, mif_REV_Complement,
         mif_REV_Reverse, mif_REV_RevComp,
         tel_FWD_Forward, tel_FWD_Complement,
         tel_FWD_Reverse, tel_FWD_RevComp, 
         tel_REV_Forward, tel_REV_Complement,
         tel_REV_Reverse, tel_REV_RevComp) %>% 
  mutate(Sequences = factor(Sequences,levels = c("neo_FWD_Forward", "neo_REV_Forward",
                                                 "neo_FWD_RevComp", "neo_REV_RevComp",
                                                 "neo_FWD_Complement", "neo_REV_Complement",
                                                 "neo_FWD_Reverse", "neo_REV_Reverse",
                                       
                                                 "mif_FWD_Forward", "mif_REV_Forward",
                                                 "mif_FWD_RevComp", "mif_REV_RevComp",
                                                 "mif_FWD_Complement", "mif_REV_Complement",
                                                 "mif_FWD_Reverse", "mif_REV_Reverse",
                                       
                                                "tel_FWD_Forward", "tel_REV_Forward",
                                                "tel_FWD_RevComp", "tel_REV_RevComp",
                                                "tel_FWD_Complement", "tel_REV_Complement",
                                                "tel_FWD_Reverse", "tel_REV_Reverse")),
         File_name = factor(File_name,levels = sample_levels),
         Run = as.factor(Run),
         Primer = factor(Primer,levels = c("NeoFish","MiFish","Teleo","NeoFish/MiFish","NeoFish/MiFish/Teleo"))) 



# PLOT 1: primers counts in reads tile plot - only primers FWD & REV, foward & revcomp ----
primers_tile_clean <- 
  primers_in_cut_reads_long %>% 
  filter(Sequences  %in% c(
    "mif_REV_RevComp", "mif_REV_Forward", "mif_FWD_RevComp", "mif_FWD_Forward",
    "neo_REV_RevComp", "neo_REV_Forward", "neo_FWD_RevComp", "neo_FWD_Forward",
    "tel_REV_RevComp", "tel_REV_Forward", "tel_FWD_RevComp", "tel_FWD_Forward")) %>% 
  filter(Run %in% c("LGC_MiniSeq_1")) %>% 
  ggplot2::ggplot(aes(y=File_name,x=Sequences,fill=log10(Count))) +
  geom_tile()+
  geom_text(aes(label = Count),size=1)+
  # scale_fill_gradient(low="white", high="darkgreen",trans="log10") +
  scale_fill_gradientn(name = "Primer counts",
                       colours = c("white","darkgreen"),
                       values = c(0,1),
                       na.value ="white") +
  theme_light(base_line_size = 1,base_size = 6) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  geom_hline(yintercept = c(40.5,82.5,86.5,116.5),color = "grey") +
  geom_vline(xintercept = c(4.5,8.5,12.5,16.5),color = "grey") +
  # coord_fixed(ratio = 0.20) +
  xlab("Primers") +
  ylab("Amostra") +
  ggtitle(label = "eDNA 1st, 2nd & 3rd runs",
              subtitle = "Primer presence on sample reads") 
# +
#   facet_wrap(~Run, drop = TRUE)


primers_tile_clean
```

<br><br>

### Quality filtering

Here the **DADA2** pipeline starts.

<br>

#### Set input libs paths

Define the paths to the libraries after *cutadapt* primer removal.

``` r
# 9 - load clean seqs to DADA2 pipe ----

all_fnFs.cut <- c(mif_fnFs.cut,neo_fnFs.cut,tel_fnFs.cut,nmt_fnFs.cut)
all_fnRs.cut <- c(mif_fnRs.cut,neo_fnRs.cut,tel_fnRs.cut,nmt_fnRs.cut)


all_fnFs.cut
all_fnRs.cut
```

<br>

#### View pre-filtering quality profiles

<br>

#### Set quality filtering output files names

``` r
# 11 - quality filter preparation ----


#name outputs
Qfilter_files <- primers_n_samples %>% 
  filter(Stage %in% c("R1 N-cleaned","R2 N-cleaned")) %>% 
  mutate(`Read file` = str_replace_all(.$`Read file`,pattern = "N-cleaned",replacement = "Qfiltered")) %>% 
  mutate(Stage = str_replace_all(.$Stage,pattern = "N-cleaned",replacement = "Qfiltered"))

primers_n_samples <- bind_rows(primers_n_samples,Qfilter_files)

#rename files so all can be traceble
names(primers_n_samples$`Read file`) <- primers_n_samples$File_name


# Qfiltered files path and names

all_filtFs <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R1 Qfiltered"]
names(all_filtFs) <- all_filtFs %>% 
    str_remove(pattern = paste0(data_path,"/Qfiltered/")) %>% 
    str_remove(pattern = "_Qfiltered.fastq.gz")
all_filtRs <- primers_n_samples$`Read file`[primers_n_samples$Stage == "R2 Qfiltered"]
names(all_filtRs) <- all_filtRs %>% 
    str_remove(pattern = paste0(data_path,"/Qfiltered/")) %>% 
    str_remove(pattern = "_Qfiltered.fastq.gz")
```

<br>

#### Quality filtering

On this step it is possible to filter by size but, as we have already
removed primers from the beginning/end of the reads, it is expected that
the remaining sequences are already trimmed to lengths compatible with
their respective amplicons. Thus, no length trimming was conducted.

``` r
# 12 - dada filtering ----

# We’ll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.

length(all_fnFs.cut)
length(all_filtFs)
length(all_fnRs.cut)
length(all_filtRs)

#sorting for next step not to mix files
all_fnFs.cut <- base::sort(all_fnFs.cut)
all_filtFs <- base::sort(all_filtFs)
all_fnRs.cut <- base::sort(all_fnRs.cut)
all_filtRs <- base::sort(all_filtRs)
names(all_fnFs.cut)
names(all_filtFs)
names(all_fnRs.cut)
names(all_filtRs)

#all

all_filtered_out <- dada2::filterAndTrim(fwd = all_fnFs.cut,
                                         filt = all_filtFs, 
                                         rev = all_fnRs.cut,
                                         filt.rev = all_filtRs,
                                         # truncLen=c(240,160),
                                         maxN=0,
                                         maxEE=c(2,2),
                                         # truncQ=2,
                                         rm.phix=TRUE,
                                         compress=TRUE, multithread=TRUE,verbose = TRUE,matchIDs = TRUE) # On Windows set multithread=FALSE
head(all_filtered_out)
```

<br>

#### View post-filtering quality profiles

``` r
#check quality profile after filtering and trimming
plotQualityProfile(all_filtFs[2])
plotQualityProfile(all_filtRs[2:8])

plotQualityProfile(mif_filtFs[])
plotQualityProfile(mif_filtRs[])
```

<br>

### Identify error rates intrinsic to sequencing

``` r
# 13 - learn error rates ----

#Learn the Error Rates
primers_n_samples$Run %>%  unique()

#run LGC_MiniSeq_1 ----
run1_errF <- learnErrors(primers_n_samples$`Read file`[primers_n_samples$Run == "LGC_MiniSeq_1" & 
                                                         primers_n_samples$Stage == "R1 Qfiltered"], 
                        multithread=TRUE,randomize = TRUE)
run1_errR <- learnErrors(primers_n_samples$`Read file`[primers_n_samples$Run == "LGC_MiniSeq_1" & 
                                                         primers_n_samples$Stage == "R2 Qfiltered"], 
                        multithread=TRUE,randomize = TRUE)

#run LGC_MiniSeq_2 ----
run2_errF <- learnErrors(primers_n_samples$`Read file`[primers_n_samples$Run == "LGC_MiniSeq_2" & 
                                                         primers_n_samples$Stage == "R1 Qfiltered"], 
                        multithread=TRUE,randomize = TRUE)
run2_errR <- learnErrors(primers_n_samples$`Read file`[primers_n_samples$Run == "LGC_MiniSeq_2" & 
                                                         primers_n_samples$Stage == "R2 Qfiltered"], 
                        multithread=TRUE,randomize = TRUE)

#run ecomol_iSeq ----
run3_errF <- learnErrors(primers_n_samples$`Read file`[primers_n_samples$Run == "ecomol_iSeq" & 
                                                         primers_n_samples$Stage == "R1 Qfiltered"], 
                        multithread=TRUE,randomize = TRUE)
run3_errR <- learnErrors(primers_n_samples$`Read file`[primers_n_samples$Run == "ecomol_iSeq" & 
                                                         primers_n_samples$Stage == "R2 Qfiltered"], 
                        multithread=TRUE,randomize = TRUE)


#   # plotErrors(run3_errF, nominalQ=TRUE)
#   # plotErrors(run3_errR, nominalQ=TRUE)
```

<br>

### Dereplication: grouping into ASVs

On this step each library is reduced to its unique composing sequences
and their counts.

``` r
# 14 - dada dereplication ----


names(primers_n_samples$`Read file`) <- primers_n_samples$File_name ###!!!!!!!!!!!!! Everything must have unique names from here

unique(primers_n_samples$`Read file`)
names(primers_n_samples$`Read file`)


#run LGC_MiniSeq_1 ----
LGC_1_derep_forward <- derepFastq(primers_n_samples$`Read file`[
  primers_n_samples$Run == "LGC_MiniSeq_1" & 
    primers_n_samples$Stage == "R1 Qfiltered"], verbose=TRUE)
# names(all_derep_forward) <- all_sample.names

LGC_1_derep_reverse <- derepFastq(primers_n_samples$`Read file`[
  primers_n_samples$Run == "LGC_MiniSeq_1" & 
    primers_n_samples$Stage == "R2 Qfiltered"], verbose=TRUE)
# names(all_derep_reverse) <- all_sample.names

LGC_1_dadaFs <- dada(LGC_1_derep_forward, err=run1_errF, multithread=TRUE)
LGC_1_dadaRs <- dada(LGC_1_derep_reverse, err=run1_errR, multithread=TRUE)

#run LGC_MiniSeq_2 ----
LGC_2_derep_forward <- derepFastq(primers_n_samples$`Read file`[
  primers_n_samples$Run == "LGC_MiniSeq_2" & 
    primers_n_samples$Stage == "R1 Qfiltered"], verbose=TRUE)
# names(all_derep_forward) <- all_sample.names

LGC_2_derep_reverse <- derepFastq(primers_n_samples$`Read file`[
  primers_n_samples$Run == "LGC_MiniSeq_2" & 
    primers_n_samples$Stage == "R2 Qfiltered"], verbose=TRUE)
# names(all_derep_reverse) <- all_sample.names

LGC_2_dadaFs <- dada(LGC_2_derep_forward, err=run2_errF, multithread=TRUE)
LGC_2_dadaRs <- dada(LGC_2_derep_reverse, err=run2_errR, multithread=TRUE)

#run ecomol_iSeq ----
ecomol_derep_forward <- derepFastq(primers_n_samples$`Read file`[
  primers_n_samples$Run == "ecomol_iSeq" & 
    primers_n_samples$Stage == "R1 Qfiltered"], verbose=TRUE)
# names(all_derep_forward) <- all_sample.names

ecomol_derep_reverse <- derepFastq(primers_n_samples$`Read file`[
  primers_n_samples$Run == "ecomol_iSeq" & 
    primers_n_samples$Stage == "R2 Qfiltered"], verbose=TRUE)
# names(all_derep_reverse) <- all_sample.names

ecomol_dadaFs <- dada(ecomol_derep_forward, err=run3_errF, multithread=TRUE)
ecomol_dadaRs <- dada(ecomol_derep_reverse, err=run3_errR, multithread=TRUE)


all_dadaFs <- c(ecomol_dadaFs,LGC_2_dadaFs,LGC_1_dadaFs)
all_dadaRs <- c(ecomol_dadaRs,LGC_2_dadaRs,LGC_1_dadaRs)

names(all_dadaFs)
```

<br>

### Merge read pairs

On this step the forward an reverse reads are merged, by overlap, in
order to reconstruct the insert full sequence. As we have samples from
three runs, they are all worked independently.

``` r
# 15 - merge read pairs ----

#run1 ----
run1_mergers <- mergePairs(dadaF = LGC_1_dadaFs,
                          derepF = LGC_1_derep_forward,
                          dadaR = LGC_1_dadaRs,
                          derepR = LGC_1_derep_reverse,
                          minOverlap = 20,
                          maxMismatch = 0,   # changed from 0 to 1 since a lot was being left out for single mismatch
                          returnRejects = TRUE,
                          verbose=TRUE)

#run2 ----
run2_mergers <- mergePairs(dadaF = LGC_2_dadaFs,
                          derepF = LGC_2_derep_forward,
                          dadaR = LGC_2_dadaRs,
                          derepR = LGC_2_derep_reverse,
                          minOverlap = 20,
                          maxMismatch = 0,   # changed from 0 to 1 since a lot was being left out for single mismatch
                          returnRejects = TRUE,
                          verbose=TRUE)

#run1 ----
run3_mergers <- mergePairs(dadaF = ecomol_dadaFs,
                          derepF = ecomol_derep_forward,
                          dadaR = ecomol_dadaRs,
                          derepR = ecomol_derep_reverse,
                          minOverlap = 20,
                          maxMismatch = 0,   # changed from 0 to 1 since a lot was being left out for single mismatch
                          returnRejects = TRUE,
                          verbose=TRUE)

all_mergers <- c(run3_mergers,run2_mergers,run1_mergers)

names(all_mergers)

length(all_dadaFs)
length(all_dadaRs)
head(all_mergers[[12]])
length(all_dadaFs)
names(all_mergers)
str(all_mergers)
class(all_mergers)


# all_seqtab <- makeSequenceTable(samples = c(run3_mergers,run2_mergers,run1_mergers))   #talvez essa função aceite varios mergers
all_seqtab <- makeSequenceTable(samples = all_mergers)
dim(all_seqtab)
View(all_seqtab)
str(all_seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(all_seqtab)))
table(nchar(getSequences(all_seqtab))) %>% plot()


names(all_dadaFs)
names(all_derep_forward)
names(all_dadaRs)
names(all_derep_reverse)
```

<br>

### Remove *chimeras*

*Chimeras* are artificial read pairs that might have been generated
erroneously on sequencing. The **DADA2** package estimates the
probability of a sequence to be chimeric given the abundancy of its
parental sequnces. After chimeric sequences removal, the remaining ASVs
length distribution is assessed. On further steps it will be used to
restrict analisys to ASVs compatible with each primer amplicons’ length
interval, in order to keep of unexpected ASVs.

``` r
# 16 - remove chimeras ----


# any(colnames(C1conc_seqtab) %in% colnames(all_seqtab))

all_seqtab.nochim <- removeBimeraDenovo(all_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(all_seqtab.nochim)
sum(all_seqtab.nochim)/sum(all_seqtab) # =  0.8404743 , perda de 16% na abundancia
#count proportion of ASVs of a given length
table(nchar(getSequences(all_seqtab.nochim)))
table(nchar(getSequences(all_seqtab.nochim))) %>% plot()
rownames(all_seqtab.nochim)



View(all_seqtab.nochim)
str(all_seqtab.nochim)
```

<br>

### Count reads and remaining ASVs

``` r
# 17 - count reads proportion throughout the pipeline ----

getN <- function(x) sum(getUniques(x))

#preparing subtables with named rows to combine latter
#raw files

names(primers_n_samples$`Read file`) <- primers_n_samples$Library 

raw_reads <- primers_n_samples %>% filter(Stage %in% c("R1","R2")) 

raw_reads_counts <- ShortRead::countFastq(dirPath = raw_reads$`Read file`) %>% as_tibble(rownames = "Read file")
raw_reads_counts <- left_join(x = raw_reads_counts, y = (raw_reads %>%  mutate(`Read file` = basename(`Read file`)) 
                                                         ),by = "Read file")

tbl_raw_FWD <- raw_reads_counts[raw_reads_counts$Stage %in% c("R1"),] %>% select(File_name, records) %>% `colnames<-`(c("File_name", "Raw FWD"))
tbl_raw_REV <- raw_reads_counts[raw_reads_counts$Stage %in% c("R2"),] %>% select(File_name, records) %>% `colnames<-`(c("File_name", "Raw REV"))



tbl_Denoised_FWD <- (sapply(all_dadaFs, getN) %>% as_tibble(rownames = "File_name")) %>% `colnames<-`(c("File_name", "Denoised FWD"))
tbl_Denoised_REV <- (sapply(all_dadaRs, getN) %>% as_tibble(rownames = "File_name")) %>% `colnames<-`(c("File_name", "Denoised REV"))
tbl_Merged <- (rowSums(all_seqtab) %>% as_tibble(rownames = "File_name")) %>% `colnames<-`(c("File_name", "Merged"))
tbl_Non_chimeric <- (rowSums(all_seqtab.nochim) %>% as_tibble(rownames = "File_name")) %>% `colnames<-`(c("File_name", "Non-chimeric"))

# combine all counts by sample to plot

all_track <- all_filtered_out %>%  as_tibble(rownames = "File_name") %>% 
  mutate(`File_name` = str_remove(string = `File_name`, pattern = "_R1_cutadapt.fastq.gz")) %>% 
  left_join(tbl_raw_FWD,by = "File_name") %>% 
  left_join(tbl_raw_REV,by = "File_name") %>% 
  left_join(tbl_Denoised_FWD,by = "File_name") %>% 
  left_join(tbl_Denoised_REV,by = "File_name") %>% 
  left_join(tbl_Merged,by = "File_name") %>% 
  left_join(tbl_Non_chimeric,by = "File_name") %>% 
  left_join(primers_n_samples[primers_n_samples$Stage == "R1",],by = "File_name") %>% 
  select(!c("Stage", "Read file"))



colnames(all_track) <- c("File_name","N-cleaned", "Filtered","Raw FWD", "Raw REV", "Denoised FWD", "Denoised REV", "Merged", "Non-Chimeric", "Type", "Group", "Library", "Primer", "Run")




# Combine tables together (if there is more than one)
track_tbl <- bind_rows(all_track)


{
all_track$File_name[all_track$File_name == "Cassaum"] <- "Positive Control\n(P.glauca)"
all_track$File_name[all_track$File_name == "Da19"] <- "Da19"
all_track$File_name[all_track$File_name == "Da20"] <- "Da20"
all_track$File_name[all_track$File_name == "Da21"] <- "Da21"
all_track$File_name[all_track$File_name == "Da22"] <- "Da22"
all_track$File_name[all_track$File_name == "Da23-mif"] <- "Da23-mif"
all_track$File_name[all_track$File_name == "Da23-neo"] <- "Da23-neo"
all_track$File_name[all_track$File_name == "neg-PCR2"] <- "neg-PCR2"
all_track$File_name[all_track$File_name == "pJequei-N-norm-M"] <- "Non-normalized JQmc\nMiFish"
all_track$File_name[all_track$File_name == "pJequei-N-norm-N"] <- "Non-normalized JQmc\nNeoFish"
all_track$File_name[all_track$File_name == "pJequei-N-norm-T"] <- "Non-normalized JQmc\nTeleo"
all_track$File_name[all_track$File_name == "pJequei-norm-M"] <- "Normalized JQmc\nMiFish"
all_track$File_name[all_track$File_name == "pJequei-norm-N"] <- "Normalized JQmc\nNeoFish"
all_track$File_name[all_track$File_name == "pJequei-norm-T"] <- "Normalized JQmc\nTeleo"
all_track$File_name[all_track$File_name == "SFJQ-mif"] <- "Normalized SFJQmc\nMiFish"
all_track$File_name[all_track$File_name == "SFJQ-neo"] <- "Normalized SFJQmc\nNeoFish"
all_track$File_name[all_track$File_name == "SFnNorm-mi"] <- "Non-normalized SFmc\nMiFish"
all_track$File_name[all_track$File_name == "SFnNorm-neo"] <- "Non-normalized SFmc\nNeoFish"
all_track$File_name[all_track$File_name == "SFNorm-mi"] <- "Normalized SFmc\nMiFish"
all_track$File_name[all_track$File_name == "SFNorm-neo"] <- "Normalized SFmc\nNeoFish"
}


# save reads counts table



writexl::write_xlsx(x = all_track,
                    path = "/home/heron/prjcts/fish_eDNA/sfjq/results/sfjq_read_counts_along_quality_control.xlsx",
                    col_names = TRUE,format_headers = TRUE)






# plot reads proportion troughout the pipeline ----


track_tbl$Primer %>% unique()
track_tbl$File_name %>% unique()

#TODO
# https://bhaskarvk.github.io/colormap/
#https://www.thinkingondata.com/something-about-viridis-library/
#set colors here ss




#perda por filtrar N
all_track %>% mutate(perda = `N-cleaned`/`Raw FWD`)


#18 - set colors for downstream plots ----

# colors 
scales::show_col(colors5)



colors5 <- c("#017504","#000791","#820000","#780058","#ff5500") #neo,mi,tel,all
colors_norm <- c("#017504","#4fc952",
                 "#000791","#3862eb",
                 "#820000","#bf4b4b")
scales::show_col(colors_norm)
scales::show_col(colors5)

#PLOT2 - sample track plot ----

# # track_tbl$File_name %>% paste0('"\n"',collapse = "") %>% cat()
# track_tbl$File_name %>% unique() %>% base::sort() %>% paste0(collapse = '\n') %>%  cat()
# sample_levels


 {
track_tbl$File_name[track_tbl$File_name == "Cassaum"] <- "Positive Control\n(P.glauca)"
track_tbl$File_name[track_tbl$File_name == "Da19"] <- "Da19"
track_tbl$File_name[track_tbl$File_name == "Da20"] <- "Da20"
track_tbl$File_name[track_tbl$File_name == "Da21"] <- "Da21"
track_tbl$File_name[track_tbl$File_name == "Da22"] <- "Da22"
track_tbl$File_name[track_tbl$File_name == "Da23-mif"] <- "Da23-mif"
track_tbl$File_name[track_tbl$File_name == "Da23-neo"] <- "Da23-neo"
#track_tbl$File_name[track_tbl$File_name == "neg-PCR2"] <-
track_tbl$File_name[track_tbl$File_name == "pJequei-N-norm-M"] <- "Non-normalized JQmc\nMiFish"
track_tbl$File_name[track_tbl$File_name == "pJequei-N-norm-N"] <- "Non-normalized JQmc\nNeoFish"
track_tbl$File_name[track_tbl$File_name == "pJequei-N-norm-T"] <- "Non-normalized JQmc\nTeleo"
track_tbl$File_name[track_tbl$File_name == "pJequei-norm-M"] <- "Normalized JQmc\nMiFish"
track_tbl$File_name[track_tbl$File_name == "pJequei-norm-N"] <- "Normalized JQmc\nNeoFish"
track_tbl$File_name[track_tbl$File_name == "pJequei-norm-T"] <- "Normalized JQmc\nTeleo"
track_tbl$File_name[track_tbl$File_name == "SFJQ-mif"] <- "Normalized SFJQmc\nMiFish"
track_tbl$File_name[track_tbl$File_name == "SFJQ-neo"] <- "Normalized SFJQmc\nNeoFish"
track_tbl$File_name[track_tbl$File_name == "SFnNorm-mi"] <- "Non-normalized SFmc\nMiFish"
track_tbl$File_name[track_tbl$File_name == "SFnNorm-neo"] <- "Non-normalized SFmc\nNeoFish"
track_tbl$File_name[track_tbl$File_name == "SFNorm-mi"] <- "Normalized SFmc\nMiFish"
track_tbl$File_name[track_tbl$File_name == "SFNorm-neo"] <- "Normalized SFmc\nNeoFish"
}




sample_levels <- c(
"Da23-mif", "Normalized SFJQmc\nMiFish",
"Da23-neo", "Normalized SFJQmc\nNeoFish",
"Da20","Non-normalized SFmc\nMiFish",
"Da19","Non-normalized SFmc\nNeoFish",
"Da22","Normalized SFmc\nMiFish",
"Da21","Normalized SFmc\nNeoFish",
"Non-normalized JQmc\nNeoFish","Non-normalized JQmc\nMiFish","Non-normalized JQmc\nTeleo",
"Normalized JQmc\nNeoFish","Normalized JQmc\nMiFish","Normalized JQmc\nTeleo",
"Positive Control\n(P.glauca)","neg-PCR2")


# Prepare counts for ploting ----

  track_tbl <- track_tbl %>%
  gather(key = "Stage",
        value = "Read Number",
        "Raw REV","Raw FWD",
        "N-cleaned", "Filtered", "Denoised FWD",
        "Denoised REV", "Merged", "Non-Chimeric") %>%
  mutate(Stage = factor(Stage, levels = c("Non-Chimeric", "Merged", "Denoised REV", "Denoised FWD", "Filtered","N-cleaned", "Raw REV","Raw FWD"))) %>%
  mutate(
    Primer = factor(Primer, levels = c("NeoFish", "MiFish", "Teleo","NeoFish/MiFish",  "NeoFish/MiFish/Teleo")),
    File_name = factor(File_name,levels = sample_levels))

  options(scipen = 22)
  
  # track_tbl %>% base::sort(track_tbl$Sample) 

# ploting ----
    
  track_plot <- track_tbl %>% 
    # filter(Run %in% c("LGC_MiniSeq_1", "LGC_MiniSeq_2")) %>% 
    # filter(Group %in% c("DNA pool")) %>% 
    ggplot(aes(y = Stage,x = `Read Number`, fill = Primer, group = Stage)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept = 300000, col = 1, linetype = 2) +
    scale_fill_manual(labels = c("NeoFish", "MiFish", "Teleo", "NeoFish/MiFish","NeoFish/MiFish/Teleo"),
                      values = alpha(colour = colors5,
                                     alpha =  0.75)) +
    labs(title = "LGC eDNA 1st & 2nd runs",
         subtitle = "Read counts per library and filtering step",
         x = "Read counts",
         y = "Data filtering step")+
    facet_wrap(~File_name,ncol = 6) +
    coord_fixed(ratio = 60000) +
    theme_bw(base_size = 7) +
    theme(legend.position = "bottom") +
    theme(axis.title = ggtext::element_markdown())

track_plot 







# save plot
ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_sample_track_plot.png",
     plot = track_plot,
     device = "png",
     width = 12,
     height = 10,
     units = "cm",
     dpi = 600)

# save plot
ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_sample_track_plot.svg",
     plot = track_plot,
     device = "svg",
     width = 12,
     height = 10,
     units = "cm",
     dpi = 600)


ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_sample_track_plot.pdf",
     plot = track_plot,
     device = "pdf",
     width = 20,
     height = 16,
     units = "cm",
     dpi = 600)
```

Reads proportions are displayed below. The experimental design intended
the same read yield for all libs, between 200K a 250K reads. The
deviations of this reange are probably due to dosage/pipetting errors.

<br>

<br>

## Classify taxonomy

On this step the ASVs identified by the **DADA2** pipeline, jointly for
all libraries of each primer, are associated (or not) to any of the
sequences on the Reference 12S Sequences Database. DADA2 has two
strategies to identify taxa. The first, *assignSpecies*, identify
perfect matches of the ASVs in the Reference Database. The second,
*assignTaxonomy*, use a RDP Naive Bayesian Classifier algorithm (Wang,
2007) with kmer size 8 and 100 bootstrap replicates to associate ASVs to
the Reference Database Sequences. In the latter, the taxonomy ranks
classification is proportional to the sequence similarity, although this
relation is not yet clear to us.

``` r
#19 - classify taxonomy exactly ----

all_sps <- dada2::assignSpecies(seqs = all_seqtab.nochim,allowMultiple = 10,
                         refFasta =  "/home/heron/prjcts/fish_eDNA/data/refs/db/LGC/jul21/dada_tax_fullDB_order_SPs_jul21.fasta",
                         tryRC=TRUE,
                         n = 20000,
                         verbose = TRUE)


#check how many ASVs were exactly identified as species

View(all_sps)

      all_csv_sp <- all_sps %>% as_tibble() %>% mutate(ASV = rownames(all_sps))
      colnames(all_csv_sp) <- c("exact Genus", "exact Species", "ASV")
```

``` r
#20 - classify taxonomy ----

all_taxa <- dada2::assignTaxonomy(seqs = all_seqtab.nochim,
                     refFasta =  "/home/heron/prjcts/fish_eDNA/data/refs/db/LGC/jul21/dada_tax_fullDB_order_jul21.fasta",
                           multithread=TRUE, tryRC=TRUE,taxLevels = c("Kingdom","Phylum","Class","Order","Family", "Genus", "Species","Specimen","Basin"),outputBootstraps = TRUE)


all_taxa$boot
```

<br><br>

Here the **DADA2** pipeline ends.

<br><br><br><br>

## Phyloseq

On this step the ASVs associated to taxonomic ranks by **DADA2** and
their respective counts by library, are combined using the **Phyloseq**
package.

<br>

### Generate sample metadata table

Here the experiment metadata is associated to each sample.

``` r
# 22 - create sample table ----

#create a sample table for each primer

# primers_n_samples$File_name

all_samdf <- unique(primers_n_samples[1:6])

samdf<- all_samdf

#rownames must me assigned in order to the next step to work

samdf <- samdf %>% as.data.frame()
rownames(samdf) <- samdf$File_name
```

<br>

This sample metadata table was created with the information available
for the samples analyzed on this first run. This table must be
customized for each experiment.

<br>

<br><br>

### **Phyloseq** data interpretation

``` r
#23 - interpret dada on phyloseq ----

all_ps <- phyloseq(otu_table(all_seqtab.nochim, taxa_are_rows = FALSE),
                   sample_data(samdf),
                   tax_table(all_taxa$tax))
                   # tax_table(all_taxa))

rownames(all_seqtab.nochim)
```

<br>

### Merge and Flex Phyloseq results

Many different graphics can be generated, together or in isolation, for
all primers/libraries and taxonomic ranks.

``` r
#24 - merge ps analisys ----

#melt phyloseq object into tbl
all_ps_tbl <- psmelt(all_ps) %>% as_tibble() %>% filter(Abundance >= 1)

colnames(all_ps_tbl)[colnames(all_ps_tbl) == "OTU"] <- "ASV"

unique(all_ps_tbl$ASV)
# unique(neo_ps_tbl$ASV)
# unique(mif_ps_tbl$ASV)

all_ps_tbl$Sample %>%  unique()
all_ps_tbl$Primer %>%  unique()

#concatenate exact species table

all_ps_tbl <- left_join(by = "ASV",x=all_ps_tbl,y= all_csv_sp)

# backup table
# all_ps_tbl_bckp <- all_ps_tbl
# all_ps_tbl <- all_ps_tbl_bckp
```

## Identify reads using NCBI BLASTn in house

``` r
# blastn ----
# Annotate all ASVs by blastN

asvs_blast <- all_ps_tbl$ASV %>% unique() %>% as.character() 

class(asvs_blast)

####### Lucio RQ - Execute shell commands #######
shell_exec <- function(cmd, .envir = parent.frame()) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    rlang::abort(message = "Package `processx` package is not installed.")
  }
  cmd_res <- processx::run(
    command = "bash",
    args = c("-c", glue::glue(cmd, .envir = .envir)), echo_cmd = FALSE
  )
  return(cmd_res)
}
########################

######### function to get fasta names from db based on subjectIDs #############
get_fasta_header <- function(id, db_path = "/data/databases/nt/nt") {
  # id <- "JQ365494.1"
  command <- "blastdbcmd -db {db_path} -entry {id} -outfmt %t"
  result <- shell_exec(cmd = command)
  return(result$stdout)
}
################################################################################

#################### function to run blast for each ASV/ASV ####################
run_blastn <- function(asv, db_path = "/data/databases/nt/nt", num_alignments = 3, num_thread = 40) {
  # blast_cmd <- "echo -e '>seq1\n{asv}' | blastn -db {db_path} -outfmt 6 -perc_identity 95 -qcov_hsp_perc 95 -num_threads {as.character(num_thread)} -num_alignments {as.character(num_alignments)}"
  blast_cmd <- "echo -e '>seq1\n{asv}' | blastn -db {db_path} -outfmt 6 -max_hsps 1 -perc_identity 95 -qcov_hsp_perc 95 -num_threads {as.character(num_thread)} -num_alignments {as.character(num_alignments)}"
  blast_res <- shell_exec(cmd = blast_cmd)
  return(blast_res)
}

## 
get_blastn_results <- function(asv, num_thread = 40) {
  blast_res <- run_blastn(asv, num_thread = num_thread)
  if (blast_res$status != 0) {
    rlang::abort(message = "Blast has not run correctly.")
  }
  `%>%` <- dplyr::`%>%`
  
  if (blast_res$stdout == "") {
    # rlang::inform(glue::glue("Sequence {asv} not found."))
    df_to_return <- tibble::tibble(`ASV` = asv)
    return(df_to_return)
  }
  blast_table <- blast_res$stdout %>%
    readr::read_delim(delim = "\t",
      col_names = c("query","subject","indentity","length","mismatches","gaps",
                  "query start","query end","subject start","subject end",
                  "e-value","bitscore"),
    trim_ws = TRUE, comment = "#"
  )
  
  blast_table$`subject header` <- purrr::map_chr(blast_table$subject, get_fasta_header)
  blast_table <- dplyr::relocate(blast_table, `subject header`)
  blast_table <- tibble::rowid_to_column(blast_table, var = "res")
  
  blast_table <- tidyr::pivot_wider(blast_table, id_cols = subject,  names_from = res, values_from = seq_len(ncol(blast_table)), names_glue = "{res}_{.value}")
  
  blast_table <- blast_table %>%
    dplyr::mutate(`ASV` = asv) %>%
    dplyr::relocate(starts_with("3_")) %>%
    dplyr::relocate(starts_with("2_")) %>%
    dplyr::relocate(starts_with("1_")) %>%
    dplyr::relocate(`ASV`)
  return(blast_table)
}
###############################################################################

# Versões paralelas
cores_to_be_used <- future::availableCores() - 2 # Usar todos os cores -2 = 78
future::plan(future::multisession(workers = cores_to_be_used))


blast_res <- furrr::future_map_dfr(asvs_blast, get_blastn_results, num_thread = 1, .options = furrr::furrr_options(seed = NULL))

# readr::write_csv(blast_res, "otu_blast_res_3_hit.csv")

# blast_res_bckp <- blast_res
# blast_res <- blast_res_bckp

nrow(blast_res)
dim(blast_res)

blast_res <- blast_res %>%  filter(`1_res` == 1 ) #remover o que não deu nada            

str(blast_res)
```

### Rename samples for plots

``` r
primers_n_samples
sample_levels

{
all_ps_tbl$File_name[all_ps_tbl$File_name == "Cassaum"] <- "Positive Control\n(P.glauca)"
all_ps_tbl$File_name[all_ps_tbl$File_name == "Da19"] <- "Non-normalized SFmc\nNeoFish B"
all_ps_tbl$File_name[all_ps_tbl$File_name == "Da20"] <- "Non-normalized SFmc\nMiFish B"
all_ps_tbl$File_name[all_ps_tbl$File_name == "Da21"] <- "Normalized SFmc\nNeoFish B"
all_ps_tbl$File_name[all_ps_tbl$File_name == "Da22"] <- "Normalized SFmc\nMiFish B"
all_ps_tbl$File_name[all_ps_tbl$File_name == "Da23-mif"] <- "SFJQ-mif B"
all_ps_tbl$File_name[all_ps_tbl$File_name == "Da23-neo"] <- "SFJQ-neo B"
all_ps_tbl$File_name[all_ps_tbl$File_name == "neg-PCR2"] <-"neg-PCR2"
all_ps_tbl$File_name[all_ps_tbl$File_name == "pJequei-N-norm-M"] <- "Non-normalized JQmc\nMiFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "pJequei-N-norm-N"] <- "Non-normalized JQmc\nNeoFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "pJequei-N-norm-T"] <- "Non-normalized JQmc\nTeleo"
all_ps_tbl$File_name[all_ps_tbl$File_name == "pJequei-norm-M"] <- "Normalized JQmc\nMiFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "pJequei-norm-N"] <- "Normalized JQmc\nNeoFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "pJequei-norm-T"] <- "Normalized JQmc\nTeleo"
all_ps_tbl$File_name[all_ps_tbl$File_name == "SFJQ-mif"] <- "Normalized SFJQmc\nMiFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "SFJQ-neo"] <- "Normalized SFJQmc\nNeoFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "SFnNorm-mi"] <- "Non-normalized SFmc\nMiFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "SFnNorm-neo"] <- "Non-normalized SFmc\nNeoFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "SFNorm-mi"] <- "Normalized SFmc\nMiFish"
all_ps_tbl$File_name[all_ps_tbl$File_name == "SFNorm-neo"] <- "Normalized SFmc\nNeoFish"
}

{
all_ps_tbl$Sample[all_ps_tbl$Sample == "Cassaum"] <- "Positive Control\n(P.glauca)"
all_ps_tbl$Sample[all_ps_tbl$Sample == "Da19"] <- "Non-normalized SFmc\nNeoFish B"
all_ps_tbl$Sample[all_ps_tbl$Sample == "Da20"] <- "Non-normalized SFmc\nMiFish B"
all_ps_tbl$Sample[all_ps_tbl$Sample == "Da21"] <- "Normalized SFmc\nNeoFish B"
all_ps_tbl$Sample[all_ps_tbl$Sample == "Da22"] <- "Normalized SFmc\nMiFish B"
all_ps_tbl$Sample[all_ps_tbl$Sample == "Da23-mif"] <- "Normalized SFJQmc\nMiFish B"
all_ps_tbl$Sample[all_ps_tbl$Sample == "Da23-neo"] <- "Normalized SFJQmc\nNeoFish B"
all_ps_tbl$Sample[all_ps_tbl$Sample == "neg-PCR2"] <- "neg-PCR2"
all_ps_tbl$Sample[all_ps_tbl$Sample == "pJequei-N-norm-M"] <- "Non-normalized JQmc\nMiFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "pJequei-N-norm-N"] <- "Non-normalized JQmc\nNeoFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "pJequei-N-norm-T"] <- "Non-normalized JQmc\nTeleo"
all_ps_tbl$Sample[all_ps_tbl$Sample == "pJequei-norm-M"] <- "Normalized JQmc\nMiFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "pJequei-norm-N"] <- "Normalized JQmc\nNeoFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "pJequei-norm-T"] <- "Normalized JQmc\nTeleo"
all_ps_tbl$Sample[all_ps_tbl$Sample == "SFJQ-mif"] <- "Normalized SFJQmc\nMiFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "SFJQ-neo"] <- "Normalized SFJQmc\nNeoFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "SFnNorm-mi"] <- "Non-normalized SFmc\nMiFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "SFnNorm-neo"] <- "Non-normalized SFmc\nNeoFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "SFNorm-mi"] <- "Normalized SFmc\nMiFish"
all_ps_tbl$Sample[all_ps_tbl$Sample == "SFNorm-neo"] <- "Normalized SFmc\nNeoFish"
}




all_ps_tbl$Group %>%  unique() %>% base::sort()


# all_ps_tbl$Group <- all_ps_tbl$Group %>% unfactor()
{
all_ps_tbl$Group[all_ps_tbl$Group %in% c("JqNnorm")] <- "Non-normalized JQmc"
all_ps_tbl$Group[all_ps_tbl$Group %in% c("JqNorm")] <- "Normalized JQmc"
all_ps_tbl$Group[all_ps_tbl$Group %in% c("SFnNorm")] <- "Non-normalized SFmc"
all_ps_tbl$Group[all_ps_tbl$Group %in% c("SFNorm")] <- "Normalized SFmc"
all_ps_tbl$Group[all_ps_tbl$Group %in% c("SFJQ")] <- "Normalized SFJQmc"
all_ps_tbl$Group[all_ps_tbl$Group %in% c("Positive control")] <- "Positive control\n(P. glauca)"
}
all_ps_tbl$Group %>% unique()

all_ps_tbl$Group <- all_ps_tbl$Group %>%  factor(levels = c("Non-normalized JQmc",
                                                            "Normalized JQmc",
                                                            "Non-normalized SFmc",
                                                            "Normalized SFmc",
                                                            "Normalized SFJQmc",
                                                            "Positive control\n(P. glauca)"))


sample_levels <- c(
"Normalized SFJQmc\nMiFish B", "Normalized SFJQmc\nMiFish",
"Normalized SFJQmc\nNeoFish B", "Normalized SFJQmc\nNeoFish",

"Normalized SFmc\nMiFish B","Normalized SFmc\nMiFish",
"Non-normalized SFmc\nMiFish B","Non-normalized SFmc\nMiFish",
"Normalized SFmc\nNeoFish B","Normalized SFmc\nNeoFish",
"Non-normalized SFmc\nNeoFish B","Non-normalized SFmc\nNeoFish",

"Non-normalized JQmc\nMiFish","Normalized JQmc\nMiFish",
"Non-normalized JQmc\nNeoFish","Normalized JQmc\nNeoFish",
"Non-normalized JQmc\nTeleo","Normalized JQmc\nTeleo",

"Positive Control\n(P.glauca)","neg-PCR2")

# all_ps_tbl_blast$Sample[all_ps_tbl_blast$Sample %in% c("pool não-normalizado\nMiFish")] 

all_ps_tbl$Sample <- all_ps_tbl$Sample %>%  factor(levels = sample_levels)


class(asvs_blast)

  
all_ps_tbl_blast <- left_join(x = all_ps_tbl,y = blast_res,by = "ASV")

colnames(all_ps_tbl_blast)

#all_ps_tbl_blast_bckp <- all_ps_tbl_blast
#all_ps_tbl_blast <- all_ps_tbl_blast_bckp
```

## ASVs seqs

``` r
#25 - recover all ASVs sequences to prepare fasta ----


#all ----
# giving our seq headers more manageable names (ASV_1, ASV_2...)
# all_asv_seqs <- tibble("ASV" = colnames(seqtab.nochim))
all_asv_seqs <- tibble("ASV" = asvs_blast)

all_asv_seqs <- all_asv_seqs %>% 
  mutate("ASV length" = nchar(ASV),
         "ASV header" = as.character(""))

all_asv_seqs <- all_asv_seqs[base::order(all_asv_seqs$`ASV length`),]
  for (i in 1:nrow(all_asv_seqs)) {

    all_asv_seqs$`ASV header`[i] <- paste0(">ASV_", i, "_", all_asv_seqs$`ASV length`[i], "bp")

  }


#combine ASV headers and all_ps_tbl
all_ps_tbl_blast <- dplyr::left_join(x = all_ps_tbl_blast,    
                               y = all_asv_seqs,
                               by = "ASV" )


# making and writing out a fasta of our final ASV seqs with tax
for (asv in 1:nrow(all_asv_seqs)) {
  
  tax <- all_ps_tbl_blast %>% 
    filter(ASV == all_asv_seqs$ASV[asv]) %>% 
    select("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Specimen") %>% 
    unique() %>% 
    paste0(collapse = "|")
  
  all_asv_seqs$`ASV header`[asv] <- paste0(all_asv_seqs$`ASV header`[asv],"_",tax)
  
  # if (condition) {
  # fazer algum teste pra ver ser ta certo
  # }
}

#write fasta file with ASVs and Taxonomy
all_asv_fasta <- c(rbind(all_asv_seqs$`ASV header`, all_asv_seqs$ASV))

write(all_asv_fasta, "~/prjcts/fish_eDNA/sfjq/results/sfjq_all_ASVs_all_primers.fasta")
```

### SWARM - ASVs to OTUs

``` r
# #swarm

asvs_abd <- all_ps_tbl_blast %>%
  group_by(`ASV`,`ASV header`) %>%
  mutate("ASV total abundance" = sum(Abundance)) %>%
  select(c(`ASV`,`ASV header`,`ASV total abundance`)) %>%
  unique() %>%
  mutate(`ASV header` = paste0(`ASV header`,"_",`ASV total abundance`))


asvs_abd$ASV %>% unique()
asvs_abd$`ASV header` %>% unique()


#write fasta file with ASVs and  abundance
# all_asv_fasta_abd <- c(rbind(asvs_abd$`ASV header primer`, asvs_abd$`ASV`))
all_asv_fasta_abd <- c(rbind(asvs_abd$`ASV header`, asvs_abd$`ASV`))

# write(all_asv_fasta_abd, "~/prjcts/fish_eDNA/sfjq/results/sfjq_ASVs_abd_primer.fasta")
write(all_asv_fasta_abd, paste0(results_path,"/sfjq_ASVs_abd.fasta"))

# ~/prjcts/fish_eDNA/sfjq/swarm$ swarm -t 50 ~/prjcts/fish_eDNA/sfjq/results/sfjq_ASVs_abd.fasta -s sfjq_swarm.stats -o sfjq_swarm.out -w sfjq_representative_OTUs.fasta -i sfjq_swarm.structure -f


swarm_clust <- readr::read_lines("~/prjcts/fish_eDNA/sfjq/swarm/sfjq_swarm.out")





asvs_abd <- asvs_abd %>% mutate("OTU"= 0)



for (asv in 1:nrow(asvs_abd)){
  for (line in 1:length(swarm_clust)) {
    if (str_detect(string =  swarm_clust[line],
                   pattern = str_remove(asvs_abd$`ASV header`[asv],
                                        pattern = ">"))) {
  asvs_abd$OTU[asv] <- line
    }
    }
}



all_ps_tbl_blast <- left_join(x = all_ps_tbl_blast,y = asvs_abd[,c(1,3,4)],by="ASV" ) 

# all_ps_tbl_blast %>% select(`final ID`,OTU) %>% View() 
# all_ps_tbl_blast %>% select(`final ID`,OTU) %>% select(OTU) %>% unique() 
# all_ps_tbl_blast %>% select(ASV,`final ID`,OTU) %>% select(ASV) %>% unique() 


all_ps_tbl_blast$OTU %>% unique()
all_ps_tbl_blast$Group %>% unique()
```

## Calculate sample abundances —-

``` r
#add ASV legth to table
# all_ps_tbl_blast_bckp2 <- all_ps_tbl_blast
# all_ps_tbl_blast <- all_ps_tbl_blast_bckp2


all_ps_tbl_blast <- all_ps_tbl_blast %>% 
  mutate("Relative abundance to all samples" = 0,
         "Relative abundance on sample" = 0,
         "Sample total abundance" = 0)

abd_total <- sum(all_ps_tbl_blast$Abundance)




all_ps_tbl_blast <- all_ps_tbl_blast %>%
  group_by(Sample) %>%
  mutate("Sample total abundance" = sum(Abundance),
         "Relative abundance to all samples" = Abundance/abd_total,
         "Relative abundance on sample" = Abundance/`Sample total abundance`) %>%
  ungroup()
```

### Set final identification from all possibilities

``` r
all_ps_tbl_blast <- all_ps_tbl_blast %>% 
  mutate(`exact GenSp` = paste(`exact Genus`,`exact Species`,sep=" "))



all_ps_tbl_blast <- all_ps_tbl_blast %>% 
  mutate("final ID" = if_else((`exact Species` %in% c(NA,"NA", "NA NA")),
                              if_else((Species %in% c(NA,"NA")),
                                      if_else(Genus %in% c(NA,"NA"),
                                              substr(as.character(`1_subject header`),1,30),
                                              Genus),
                                      Species),
                              as.character(`exact GenSp`)))
```

\#Group/correct species for ploting

``` r
# all_ps_tbl_blast_bckp3 <- all_ps_tbl_blast

#Species detected
all_ps_tbl_blast$`final ID` %>% unique() %>% base::sort()
{
all_ps_tbl_blast$`final ID`[(all_ps_tbl_blast$`final ID` %in% c("Astyanax aff_fasciatus","Astyanax cf_fasciatus"))] <- "Astyanax fasciatus"
all_ps_tbl_blast$`final ID`[(all_ps_tbl_blast$`final ID` %in% c("Astyanax cf_lacustris"))] <- "Astyanax lacustris"
all_ps_tbl_blast$`final ID`[(all_ps_tbl_blast$`final ID` %in% c("Characidium sp"))] <- "Characidium"
all_ps_tbl_blast$`final ID`[(all_ps_tbl_blast$`final ID` %in% c("Hypostomus sp"))] <- "Hypostomus"
all_ps_tbl_blast$`final ID`[(all_ps_tbl_blast$`final ID` %in% c("Hoplias malabaricus/sp"))] <- "Hoplias malabaricus"
all_ps_tbl_blast$`final ID`[(all_ps_tbl_blast$`final ID` %in% c("Rhamdia aff_quelen"))] <- "Rhamdia quelen"
all_ps_tbl_blast$`final ID`[(all_ps_tbl_blast$`final ID` %in% c("Coptodon zillii KAUM:I:90126 m"))] <- "Coptodon zillii"
all_ps_tbl_blast$`final ID`[(all_ps_tbl_blast$`final ID` %in% c("Trachelyopterus cf_galeatus/galeatus"))] <- "Trachelyopterus galeatus"
}
```

### Identify primers expected ASV legth range

``` r
# create ranges column
all_ps_tbl_blast <- all_ps_tbl_blast %>%
   mutate("Expected length" = "FALSE")


# fill ranges column with expected primer insert ranges

for (asv in 1:nrow(all_ps_tbl_blast)) {
   if (all_ps_tbl_blast$Primer[asv] == "NeoFish") {
      if (all_ps_tbl_blast$`ASV length`[asv] >= 185 && all_ps_tbl_blast$`ASV length`[asv] <= 200) {
         all_ps_tbl_blast$`Expected length`[asv] <- "in range"
      }else{
         all_ps_tbl_blast$`Expected length`[asv] <- "out of range"
         }

   }
   if (all_ps_tbl_blast$Primer[asv] == "MiFish") {
      if (all_ps_tbl_blast$`ASV length`[asv] >= 165 && all_ps_tbl_blast$`ASV length`[asv] <= 180) {
         all_ps_tbl_blast$`Expected length`[asv] <- "in range"
      }else{
         all_ps_tbl_blast$`Expected length`[asv] <- "out of range"
      }

   }
  if (all_ps_tbl_blast$Primer[asv] == "Teleo") {
      if (all_ps_tbl_blast$`ASV length`[asv] >= 60 && all_ps_tbl_blast$`ASV length`[asv] <= 75) {
         all_ps_tbl_blast$`Expected length`[asv] <- "in range"
      }else{
         all_ps_tbl_blast$`Expected length`[asv] <- "out of range"
      }

   }
  if (all_ps_tbl_blast$Primer[asv] == "NeoFish/MiFish") {
      if (all_ps_tbl_blast$`ASV length`[asv] >= 165 && all_ps_tbl_blast$`ASV length`[asv] <= 200) {
         all_ps_tbl_blast$`Expected length`[asv] <- "in range"
      }else{
         all_ps_tbl_blast$`Expected length`[asv] <- "out of range"
      }

   }
  if (all_ps_tbl_blast$Primer[asv] == "NeoFish/MiFish/Teleo") {
      if (all_ps_tbl_blast$`ASV length`[asv] %in% c(60:75,165:180,185:200)) {
         all_ps_tbl_blast$`Expected length`[asv] <- "in range"
      }else{
         all_ps_tbl_blast$`Expected length`[asv] <- "out of range"
      }

   }

}


#factorize comlumn
all_ps_tbl_blast$`Expected length` <- as.factor(all_ps_tbl_blast$`Expected length`)
```

\#Reorder table

``` r
paste0(colnames(all_ps_tbl_blast),"\n") %>%  cat()


# all_ps_tbl_blast_bckp4 <- all_ps_tbl_blast
# all_ps_tbl_blast <- all_ps_tbl_blast_bckp4


all_ps_tbl_blast <- 
  all_ps_tbl_blast %>% 
  select(c("Sample","Group","Type","Primer","File_name","Library","Run",
           "final ID",
           "Abundance",
           "Relative abundance to all samples",
           "Relative abundance on sample",
           "Sample total abundance",
           "Kingdom","Phylum","Class","Order","Family",
           "Genus","Species","Specimen","Basin",
           "exact Genus","exact Species",
           "exact GenSp",
           "1_subject header","1_subject",
           "1_indentity","1_length",
           # "1_mismatches","1_gaps",
           # "1_query start","1_query end","1_subject start",
           # "1_subject end","1_e-value","1_bitscore",
           "2_subject header","2_subject",
           "2_indentity","2_length",
           # "2_mismatches","2_gaps",
           # "2_query start","2_query end","2_subject start",
           # "2_subject end","2_e-value","2_bitscore",
           "3_subject header","3_subject",
           "3_indentity","3_length",
           # "3_mismatches","3_gaps",
           # "3_query start","3_query end","3_subject start",
           # "3_subject end","3_e-value","3_bitscore",
           "ASV","ASV length","ASV header","Expected length","OTU"
           ))

# paste0(colnames(all_ps_tbl_blast),"\n") %>%  cat()
names(all_ps_tbl_blast)[which(names(all_ps_tbl_blast)=="ASV")] <- "ASV (Sequence)"
names(all_ps_tbl_blast)[which(names(all_ps_tbl_blast)== "ASV length")] <- "ASV size (pb)"
```

\#\#\#save complete table

``` r
#order by abundance

smp_abd_ID <- all_ps_tbl_blast[rev(base::order(all_ps_tbl_blast$Abundance)),] %>% 
  filter(`Abundance` > 0) 

dim(smp_abd_ID)

writexl::write_xlsx(x = smp_abd_ID,
                    path = "/home/heron/prjcts/fish_eDNA/sfjq/results/sfjq_all_analysis_info_06-03-22.xlsx",
                    col_names = TRUE,format_headers = TRUE)














ASVs_per_sample <- all_ps_tbl_blast %>% 
  # filter(Run %in% c("LGC_MiniSeq_1", "LGC_MiniSeq_2")) %>%
  filter(!`final ID` %in% c(NA,"NA")) %>% 
  mutate("Sample" = factor(Sample,levels = sample_levels)) %>%
  group_by(Sample) %>%
   summarize("Library" = unique(`Library`),
     "Group" = unique(Group),
             "Primer" = unique(Primer),
             "Total ASV" = length(unique(`ASV (Sequence)`[Abundance != 0])),
             "ASVs out of range" = length(unique(`ASV (Sequence)`[Abundance != 0 & `Expected length` == "out of range"])),
             "ASVs in range" = length(unique(`ASV (Sequence)`[Abundance != 0 & `Expected length` == "in range"]))
             ,
             "Identified Species" = length(unique(`final ID`[Abundance != 0 & `Expected length` == "in range"]))
             ) 


writexl::write_xlsx(x = ASVs_per_sample,
                    path = "/home/heron/prjcts/fish_eDNA/sfjq/results/sfjq_ASVs_per_sample_06-07-21.xlsx",
                    col_names = TRUE,format_headers = TRUE)
```

## curing: manual checking of the species assignment

``` r
all_ps_tbl_bl_cur <- smp_abd_ID



all_ps_tbl_bl_cur <- all_ps_tbl_bl_cur %>% 
  mutate("revised final ID" = `final ID`)


all_ps_tbl_bl_cur$`revised final ID` %>% unique() %>%  sort()




#correct misidentifications one by one
{
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("Prochilodus")] <- "Prochilodus argenteus/hartii"
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("Prochilodus argenteus")] <- "Prochilodus argenteus/hartii"
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("Pimelodus")] <- "Pimelodus pohli"
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("Astyanax")] <- "Astyanax lacustris"
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("Hypostomus")] <- "Hypostomus alatus"
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("NA elegans/gilbert")] <- "Cyphocharax gilbert"
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("NA lepidura/xenodon")] <- "Curimatella lepidura"
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("Roeboides xenodon")] <- "Curimatella lepidura"
all_ps_tbl_bl_cur$`revised final ID`[all_ps_tbl_bl_cur$`revised final ID` %in% c("Curimatella lepidura")] <- "Roeboides xenodon"





all_ps_tbl_bl_cur$`revised final ID`[(all_ps_tbl_bl_cur$`revised final ID` %in% c("Hoplias brasiliensis/intermedius")) & 
                                       all_ps_tbl_bl_cur$Group %in% c("Non-normalized SFmc","Normalized SFmc")] <- "Hoplias intermedius"


all_ps_tbl_bl_cur$`revised final ID`[(all_ps_tbl_bl_cur$`revised final ID` %in% c("Hoplias brasiliensis/intermedius")) & 
                                       all_ps_tbl_bl_cur$Group %in% c("Non-normalized JQmc","Normalized JQmc")] <- "Hoplias brasiliensis"



all_ps_tbl_bl_cur$`revised final ID`[(all_ps_tbl_bl_cur$`revised final ID` %in% c("Hoplias intermedius")) & 
                                       all_ps_tbl_bl_cur$Group %in% c("Normalized SFJQmc")] <- "Hoplias brasiliensis/intermedius"




all_ps_tbl_bl_cur$`revised final ID`[(all_ps_tbl_bl_cur$`revised final ID` %in% c("Hoplias intermedius")) & 
                                       all_ps_tbl_bl_cur$Group %in% c("Non-normalized JQmc","Normalized JQmc")] <- "Hoplias brasiliensis"







}






# all_ps_tbl_bl_cur[all_ps_tbl_bl_cur$`revised final ID` %in% c("Pimelodus"),] %>% select(`ASV (Sequence)`) 
# <- "Prochilodus argenteus/hartii"
```

### Identify expected species

As this was a controlled experiment, we must identify the species that
were expected and those who were not.

``` r
#duplicate final ID so we can identify partially identified species

# all_ps_tbl_blast_bckp5 <- all_ps_tbl_blast 
# all_ps_tbl_blast <- all_ps_tbl_blast_bckp5
# 
#                           all_ps_tbl_blast <- all_ps_tbl_bl_cur %>%
#                             mutate(`revised final ID` = `final ID`)



# expctd_sps_tbl <- read.csv(file = "~/prjcts/fish_eDNA/sfjq/data/sfjq_species.csv",
# expctd_sps_tbl <- read.csv(file = "~/prjcts/fish_eDNA/sfjq/data/sfjq_species_06mar22.csv",
# expctd_sps_tbl <- read.csv(file = "/home/heron/prjcts/fish_eDNA/sfjq/data/sfjq_species_10mar22.csv",
expctd_sps_tbl <- read.csv(file = "/home/heron/prjcts/fish_eDNA/sfjq/data/sfjq_species_02jun22.csv",
                           header = TRUE, check.names = FALSE) %>% as_tibble() 

#mudei o T. chalceus de 0.0073 para 0.0299 no pool SF

colnames(expctd_sps_tbl)[colnames(expctd_sps_tbl) == "Species"] <- "revised final ID"

# expctd_sps_tbl <- expctd_sps_tbl %>% 
#   mutate(`Percentage on respective Non-norm pool` = as.numeric(`Percentage on respective Non-norm pool`),
#          `Percentage on respective Norm pool` = as.numeric(`Percentage on respective Norm pool`))

# all_ps_tbl_blast$`final ID` %>% unique() %>% base::sort()%>% paste0(collapse = '", \n"') %>% cat()
# all_ps_tbl_blast$`revised final ID` %>% unique() %>% base::sort()%>% paste0(collapse = '", \n"') %>% cat()

# all_ps_tbl_blast$`revised final ID` %>% unique()
# all_ps_tbl_blast$`revised final ID` %>% base::sort() %>% unique()
```

### Proportions table

As this was a controlled experiment, we must identify the species that
were expected and those who were not.

``` r
# tabelas de proporções ----
# JQmc ----

# criando a tabela de proporções do JQ


# if used left_join, p hartii wont show up
all_ps_tbl_jq <- full_join(
  all_ps_tbl_bl_cur[all_ps_tbl_bl_cur$Sample %in% c(
  "Normalized JQmc\nNeoFish", "Non-normalized JQmc\nNeoFish",
  "Normalized JQmc\nMiFish", "Non-normalized JQmc\nMiFish",
  "Normalized JQmc\nTeleo", "Non-normalized JQmc\nTeleo"),],
  expctd_sps_tbl[expctd_sps_tbl$Pool=="JQ",c(2:9)],
          by = "revised final ID") 

# all_ps_tbl_jq %>% select(`revised final ID`,Abundance,`Relative abundance on sample`, 64:73) %>% View()

all_ps_tbl_jq$`Expected Species` <- "not expected"

all_ps_tbl_jq$`revised final ID` %>% unique()

all_ps_tbl_jq <- all_ps_tbl_jq %>% mutate( 
  # `Expected Species`=if_else((!is.na(`Full name`) & (.$`revised final ID` %in% expected_sps)),"expected","not expected")
  `Expected Species`=if_else((!is.na(`Full name`)),"expected","not expected")
  )

colnames(all_ps_tbl_jq)


jq_proportions <- all_ps_tbl_jq %>%
  select(c(1:24,37:50)) %>%
  pivot_longer(cols = c("Percentage on respective Norm pool", "Percentage on respective Non-norm pool" ),
               names_to = "Pool",values_to = "Proportion") %>% 
  # mutate(Proportion = Proportion*100) %>% 
  mutate(Sample = factor(Sample,levels = c(
    "Normalized JQmc\nNeoFish", "Non-normalized JQmc\nNeoFish",
    "Normalized JQmc\nMiFish", "Non-normalized JQmc\nMiFish",
    "Normalized JQmc\nTeleo", "Non-normalized JQmc\nTeleo"
    # ,
    # "Non-normalized SFmc\nMiFish","Non-normalized SFmc\nNeoFish",
    # "Normalized SFmc\nMiFish","Normalized SFmc\nNeoFish",
    # "Normalized SFJQmc\nMiFish", "Normalized SFJQmc\nNeoFish"
    ))) %>%
  filter((Sample %in% c("Non-normalized JQmc\nNeoFish", "Non-normalized JQmc\nMiFish", "Non-normalized JQmc\nTeleo"
                        # , "Non-normalized SFmc\nMiFish","Non-normalized SFmc\nNeoFish"
                        ) & Pool %in% c("Percentage on respective Non-norm pool"
                                        ) )|(Sample %in% c("Normalized JQmc\nNeoFish", "Normalized JQmc\nMiFish", "Normalized JQmc\nTeleo"
                                                           # , "Normalized SFmc\nMiFish", "Normalized SFmc\nNeoFish", "Normalized SFJQmc\nMiFish", "Normalized SFJQmc\nNeoFish"
                                                           ) & Pool %in% c("Percentage on respective Norm pool") )) %>% 
  group_by(Sample,`revised final ID`) %>%
  summarize(Proportion = unique(Proportion),
            `Relative abundance on sample` = sum(`Relative abundance on sample`),
            Sample = unique(Sample),
            # `revised final ID` = unique(`revised final ID`),
            `revised final ID` = unique(`revised final ID`),
            Primer = unique(Primer),
            `Expected Species` = unique(`Expected Species`),
            Pool = unique(Pool),
            `Num ASVs` = length(`ASV (Sequence)`),
            `Num OTUs` = length(unique(OTU))) %>% 
  ungroup()


#SFmc ----

# criando a tabela de proporções do 

all_ps_tbl_sf <- full_join(
  all_ps_tbl_bl_cur[all_ps_tbl_bl_cur$Sample %in% c(
  "Normalized SFmc\nNeoFish", "Non-normalized SFmc\nNeoFish",
  "Normalized SFmc\nMiFish", "Non-normalized SFmc\nMiFish"),],
  expctd_sps_tbl[expctd_sps_tbl$Pool=="SF",c(2:9)],
          by = "revised final ID")

# all_ps_tbl_jq %>% select(`revised final ID`,Abundance,`Relative abundance on sample`, 64:73) %>% View()

all_ps_tbl_sf$`Expected Species` <- "not expected"

all_ps_tbl_sf$`revised final ID` %>% unique()

all_ps_tbl_sf <- all_ps_tbl_sf %>% mutate( 
  # `Expected Species`=if_else((!is.na(`Full name`) & (.$`revised final ID` %in% expected_sps)),"expected","not expected")
  `Expected Species`=if_else((!is.na(`Full name`)),"expected","not expected")
  )

colnames(all_ps_tbl_sf)


sf_proportions <- all_ps_tbl_sf %>%
  # select(c(1:15,63:68,73:74)) %>% 
  select(c(1:24,37:50)) %>%
  pivot_longer(cols = c("Percentage on respective Norm pool", "Percentage on respective Non-norm pool" ),
               names_to = "Pool",values_to = "Proportion") %>% 
  # mutate(Proportion = Proportion*100) %>% 
  mutate(Sample = factor(Sample,levels = c(
    "Normalized SFmc\nNeoFish", "Non-normalized SFmc\nNeoFish",
    "Normalized SFmc\nMiFish", "Non-normalized SFmc\nMiFish"
    # ,
    # "Non-normalized SFmc\nMiFish","Non-normalized SFmc\nNeoFish",
    # "Normalized SFmc\nMiFish","Normalized SFmc\nNeoFish",
    # "Normalized SFJQmc\nMiFish", "Normalized SFJQmc\nNeoFish"
    ))) %>%
  filter((Sample %in% c("Non-normalized SFmc\nNeoFish", "Non-normalized SFmc\nMiFish"
                        # , "Non-normalized SFmc\nMiFish","Non-normalized SFmc\nNeoFish"
                        ) & Pool %in% c("Percentage on respective Non-norm pool"
                                        ) )|(Sample %in% c("Normalized SFmc\nNeoFish", "Normalized SFmc\nMiFish"
                                                           # , "Normalized SFmc\nMiFish", "Normalized SFmc\nNeoFish", "Normalized SFJQmc\nMiFish", "Normalized SFJQmc\nNeoFish"
                                                           ) & Pool %in% c("Percentage on respective Norm pool") )) %>% 
  group_by(Sample,`revised final ID`) %>% 
  summarize(Proportion = unique(Proportion),
            `Relative abundance on sample` = sum(`Relative abundance on sample`),
            Sample = unique(Sample),
            # `revised final ID` = unique(`revised final ID`),
            `revised final ID` = unique(`revised final ID`),
            Primer = unique(Primer),
            `Expected Species` = unique(`Expected Species`),
            Pool = unique(Pool),
            `Num ASVs` = length(`ASV (Sequence)`),
            `Num OTUs` = length(unique(OTU))) %>% 
  ungroup()


#SFJQmc ----
# View(expctd_sps_tbl[expctd_sps_tbl$Pool=="SFJQ",c(2:9)])
# criando a tabela de proporções do SFJQ

all_ps_tbl_sfjq <- full_join(
  all_ps_tbl_bl_cur[all_ps_tbl_bl_cur$Sample %in% c(
  "Normalized SFJQmc\nNeoFish", "Normalized SFJQmc\nMiFish"),],
  expctd_sps_tbl[expctd_sps_tbl$Pool=="SFJQ",c(2:9)],
          by = "revised final ID")


# all_ps_tbl_jq %>% select(`revised final ID`,Abundance,`Relative abundance on sample`, 64:73) %>% View()

all_ps_tbl_sfjq$`Expected Species` <- "not expected"

all_ps_tbl_sfjq$`revised final ID` %>% unique()

all_ps_tbl_sfjq <- all_ps_tbl_sfjq %>% mutate( 
  # `Expected Species`=if_else((!is.na(`Full name`) & (.$`revised final ID` %in% expected_sps)),"expected","not expected")
  `Expected Species`=if_else((!is.na(`Full name`)),"expected","not expected")
  )

colnames(all_ps_tbl_sfjq)


sfjq_proportions <- all_ps_tbl_sfjq %>%
  select(c(1:24,37:50)) %>%
  # select(c(1:15,63:68,73:74)) %>% 
  pivot_longer(cols = c("Percentage on respective Norm pool", "Percentage on respective Non-norm pool" ),
               names_to = "Pool",values_to = "Proportion") %>% 
  # mutate(Proportion = Proportion*100) %>% 
  mutate(Sample = factor(Sample,levels = c(
    "Normalized SFJQmc\nNeoFish", "Normalized SFJQmc\nMiFish"
    ))) %>%
  filter((Sample %in% c("Normalized SFJQmc\nNeoFish", "Normalized SFJQmc\nMiFish"
                        ) & Pool %in% c("Percentage on respective Norm pool"
                                        ) )) %>% 
  group_by(Sample,`revised final ID`) %>% 
  summarize(Proportion = unique(Proportion),
            `Relative abundance on sample` = sum(`Relative abundance on sample`),
            Sample = unique(Sample),
            # `revised final ID` = unique(`revised final ID`),
            `revised final ID` = unique(`revised final ID`),
            Primer = unique(Primer),
            `Expected Species` = unique(`Expected Species`),
            Pool = unique(Pool),
            `Num ASVs` = length(`ASV (Sequence)`),
            `Num OTUs` = length(unique(OTU))) %>% 
  ungroup()



#this table will be used
all_ps_tbl_sfjq_full <-  dplyr::bind_rows(all_ps_tbl_sfjq,all_ps_tbl_sf,all_ps_tbl_jq)



all_ps_tbl_sfjq_full %>% colnames() %>% unique() %>% paste0(collapse = '",\n"') %>% cat()

# all_ps_tbl_sfjq_full <- all_ps_tbl_sfjq_full %>% 
#   select(c(
#     "Sample",
#     "Group",
#     "Type",
#     "Primer",
#     "File_name",
#     "Library",
#     "Run",
#     "final ID",
#     "Abundance",
#     "Relative abundance to all samples",
#     "Relative abundance on sample",
#     "Sample total abundance",
#     "Kingdom",
#     "Phylum",
#     "Class",
#     "Order",
#     "Family",
#     "Genus",
#     "Species",
#     "Specimen",
#     "Basin",
#     "exact Genus",
#     "exact Species",
#     "exact GenSp",
#     # "1_subject header",
#     # "1_subject",
#     # "1_indentity",
#     # "1_length",
#     # "2_subject header",
#     # "2_subject",
#     # "2_indentity",
#     # "2_length",
#     # "3_subject header",
#     # "3_subject",
#     # "3_indentity",
#     # "3_length",
#     "ASV (Sequence)",
#     "ASV size (pb)",
#     "ASV header",
#     "Expected length",
#     "OTU",
#     "revised final ID",
#     "Full name",
#     "LGCdb ID",
#     "DNA concentration QUBIT (ng/ul)",
#     "Input DNA",
#     "Total DNA in NN pool",
#     "Percentage on respective Non-norm pool",
#     "Percentage on respective Norm pool",
#     "Expected Species") )
```

### Pearson correlations between DNA input and sequence yield

``` r
#JQmc ----
# {
# jq_df_neo <- jq_proportions %>% 
#   # filter(`Expected Species` %in% c("expected")) %>% 
#   pivot_wider(names_from = Pool,values_from = `Relative abundance on sample`,id_cols = c(
#   `final ID`,Primer),) %>% filter(Primer %in% c("NeoFish")) %>% 
#   select(c(1,3,4)) %>% as.data.frame() %>% na.exclude()
#   
# 
# jq_df_mif <- jq_proportions %>% 
#   filter(`Expected Species` %in% c("expected")) %>% 
#   pivot_wider(names_from = Pool,values_from = `Relative abundance on sample`,id_cols = c(
#   `final ID`,Primer),) %>% filter(Primer %in% c("MiFish")) %>% 
#   select(c(1,3,4)) %>% as.data.frame() %>% na.exclude()
#   
# 
# jq_df_tel <- jq_proportions %>% 
#   filter(`Expected Species` %in% c("expected")) %>% pivot_wider(names_from = Pool,values_from = `Relative abundance on sample`,id_cols = c(
#   `final ID`,Primer),) %>% filter(Primer %in% c("Teleo")) %>% 
#   select(c(1,3,4)) %>% as.data.frame() %>% na.exclude()
# }
# 
# 
# 
# rownames(jq_df_neo) <-jq_df_neo$`final ID`
# rownames(jq_df_mif) <-jq_df_mif$`final ID`
# rownames(jq_df_tel) <-jq_df_tel$`final ID`
# 
# #neo
# cor.test(x = jq_df_neo$`Percentage on respective Norm pool`,y=jq_df_neo$`Percentage on respective Non-norm pool` ,method = "pearson")
# #mif
# cor.test(x = jq_df_mif$`Percentage on respective Norm pool`,y=jq_df_mif$`Percentage on respective Non-norm pool` ,method = "pearson")
# #teleo
# cor.test(x = jq_df_tel$`Percentage on respective Norm pool`,y=jq_df_tel$`Percentage on respective Non-norm pool` ,method = "pearson")
#   
#   
#   
#   
# 
# 
# #SFmc ----
# 
# sf_df_neo <- sf_proportions %>% 
#   filter(`Expected Species` %in% c("expected")) %>% 
#   pivot_wider(names_from = Pool,values_from = `Relative abundance on sample`,id_cols = c(
#   `final ID`,Primer),) %>% filter(Primer %in% c("NeoFish")) %>% 
#   select(c(1,3,4)) %>% as.data.frame() %>% na.exclude()
#   
# 
# sf_df_mif <- sf_proportions %>% 
#   filter(`Expected Species` %in% c("expected")) %>% 
#   pivot_wider(names_from = Pool,values_from = `Relative abundance on sample`,id_cols = c(
#   `final ID`,Primer),) %>% filter(Primer %in% c("MiFish")) %>% 
#   select(c(1,3,4)) %>% as.data.frame() %>% na.exclude()
#   
# 
# 
# 
# rownames(sf_df_neo) <-sf_df_neo$`final ID`
# rownames(sf_df_mif) <-sf_df_mif$`final ID`
# 
# #neo
# cor.test(x = sf_df_neo$`Percentage on respective Norm pool`,y=sf_df_neo$`Percentage on respective Non-norm pool` ,method = "pearson")
# #mif
# cor.test(x = sf_df_mif$`Percentage on respective Norm pool`,y=sf_df_mif$`Percentage on respective Non-norm pool` ,method = "pearson")
# 
# 
# # ------------------------




#correlação entre o DNA input e reads ABD -----

# JQmc
{
jq_df_neo_norm <- jq_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Normalized JQmc\nNeoFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

jq_df_mif_norm <- jq_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Normalized JQmc\nMiFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

jq_df_tel_norm <- jq_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Normalized JQmc\nTeleo")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

jq_df_neo_skew <- jq_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Non-normalized JQmc\nNeoFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

jq_df_mif_skew <- jq_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Non-normalized JQmc\nMiFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

jq_df_tel_skew <- jq_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Non-normalized JQmc\nTeleo")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

#SF 
sf_df_neo_norm <- sf_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Normalized SFmc\nNeoFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

sf_df_mif_norm <- sf_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Normalized SFmc\nMiFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

sf_df_neo_skew <- sf_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Non-normalized SFmc\nNeoFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

sf_df_mif_skew <- sf_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Non-normalized SFmc\nMiFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

# SFJQ
sfjq_df_neo_norm <- sfjq_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Normalized SFJQmc\nNeoFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()

sfjq_df_mif_norm <- sfjq_proportions %>% 
  ungroup() %>% 
  filter(Sample %in% c("Normalized SFJQmc\nMiFish")) %>% 
  select(c(2,3,4)) %>% as.data.frame() %>% na.exclude()
}




#correlações 
#JQmc
cor.test(x = jq_df_neo_norm$Proportion ,
         y = jq_df_neo_norm$`Relative abundance on sample` ,
         method = "pearson")


cor.test(x = jq_df_neo_skew$Proportion ,
         y = jq_df_neo_skew$`Relative abundance on sample` ,
         method = "pearson")



cor.test(x = jq_df_mif_norm$Proportion ,
         y = jq_df_mif_norm$`Relative abundance on sample` ,
         method = "pearson")


cor.test(x = jq_df_mif_skew$Proportion ,
         y = jq_df_mif_skew$`Relative abundance on sample` ,
         method = "pearson")


cor.test(x = jq_df_tel_norm$Proportion ,
         y = jq_df_tel_norm$`Relative abundance on sample` ,
         method = "pearson")


cor.test(x = jq_df_tel_skew$Proportion ,
         y = jq_df_tel_skew$`Relative abundance on sample` ,
         method = "pearson")

#SF
cor.test(x = sf_df_neo_norm$Proportion ,
         y = sf_df_neo_norm$`Relative abundance on sample` ,
         method = "pearson",)


cor.test(x = sf_df_neo_skew$Proportion ,
         y = sf_df_neo_skew$`Relative abundance on sample` ,
         method = "pearson")



cor.test(x = sf_df_mif_norm$Proportion ,
         y = sf_df_mif_norm$`Relative abundance on sample` ,
         method = "pearson")


cor.test(x = sf_df_mif_skew$Proportion ,
         y = sf_df_mif_skew$`Relative abundance on sample` ,
         method = "pearson")
```

## Richness analysis on Vegan

``` r
library(vegan)
# data(dune)
# decorana(dune)

# class(dune)
#1- prepare data for entry in vegan ----

all_ps_blst_vegan <- all_ps_tbl_bl_cur %>% 
  filter(`Expected length` %in% c("in range")) %>% 
  mutate("Normalization" = str_split(.$Group,pattern = " ",2,simplify = TRUE)[,1],
         "Mock Community" = str_split(.$Group,pattern = " ",2,simplify = TRUE)[,2]) %>% 
  filter(!(Sample %in% c("Positive Control\n(P.glauca)","neg-PCR2"))) %>%   #remove control samples
  select(c(Sample,Group,Normalization,`Mock Community`,Type,Primer,File_name,Library,Run,`final ID`,`Relative abundance on sample`)) %>% 
  group_by(Sample,`final ID`,Group,Type,Primer,File_name,Library,Run,Normalization,`Mock Community`) %>% 
  summarise(`Relative abundance on sample` = sum(`Relative abundance on sample`)) %>% 
  pivot_wider(c(Sample,Group,Type,Primer,File_name,Library,Run,Normalization,`Mock Community`),names_from = `final ID` ,values_from = `Relative abundance on sample`) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  # mutate(Library = unfactor(Library)) %>% 
  mutate("Sample number" = 0) %>% 
  ungroup()  %>% 
  select(`Sample number`, 1:(ncol(.)-1)) %>% 
  mutate(Normalization = factor(Normalization))

#2- associate sample numbers to sample names ----
for (sample in 1:nrow(all_ps_blst_vegan)) {
  all_ps_blst_vegan$`Sample number`[sample] <- sample 
  
}

#tirando as amostras da ecomol pra facilitar

all_ps_blst_vegan <- all_ps_blst_vegan[all_ps_blst_vegan$Run %in% c("LGC_MiniSeq_1", "LGC_MiniSeq_2"),] 




colnames(all_ps_blst_vegan)
hist(colSums(all_ps_blst_vegan[,-c(1:10)]))
hist(rowSums(all_ps_blst_vegan[,-c(1:10)]))
all_ps_blst_vegan[,-c(1:10)]

all_ps_blst_vegan %>% select(Sample, `Sample number`)
# all_ps_blst_vegan %>% select(`Sample number`, 1:(ncol(.)-1))

#3- create data.frame of species counts: rownames are Sample numbers ----

all_ps_blst_vegan_df <- all_ps_blst_vegan %>% 
  # filter(Run %in% c("LGC_MiniSeq_1", "LGC_MiniSeq_2")) %>% 
  select(-c("Sample", "Group", "Type", "Primer", "File_name", "Library", "Run","Normalization","Mock Community")) %>% 
  select(base::sort(colnames(.))) %>% 
  as.data.frame() 

#4- name rows as Sample numbers and remove column ----
row.names(all_ps_blst_vegan_df) <- all_ps_blst_vegan_df$`Sample number`
all_ps_blst_vegan_df <- all_ps_blst_vegan_df %>% 
  select(-c(`Sample number`))

library(vegan)
#5- 

all_ps_ord <- decorana(veg = all_ps_blst_vegan_df)

all_ps_ord %>% summary()

all_ps_ord %>% str()

all_ps_ord$cproj


plot(all_ps_ord)
plot(all_ps_ord,type = "p")
plot(all_ps_ord,type = "c") 

points(all_ps_ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(all_ps_ord, display = "sites", cex=0.7, col="blue")
text(all_ps_ord, display = "spec", cex=0.7, col="blue")



#6- NMDS analisys ----



# library(vegan)
# data(varespec)
#6a- Calculate distances ----
all_ps_vg_dist <- vegdist(all_ps_blst_vegan_df, method="bray")

vegan::scores(all_ps_vg_dist)

# all_ps_vg_dist_metaMDS <- metaMDS(comm = all_ps_vg_dist, autotransform = FALSE) 
# actually autotransform = FALSE doesn't seem to change the results

# plot(all_ps_vg_dist_metaMDS)

# all_ps_vg_dist_metaMDS_2 <- metaMDS(comm = all_ps_vg_dist, distance = "bray", k =2)

# plot(all_ps_vg_dist_metaMDS_2)

#selecionar apenas espécies esperadas?

all_ps_blst_vegan_df %>% ncol()
all_ps_blst_vegan_df <- all_ps_blst_vegan_df[,(colnames(all_ps_blst_vegan_df) %in% expected_sps)]


all_ps_blst_vegan_df %>% ncol()
all_ps_vg_dist <- vegdist(all_ps_blst_vegan_df, method="bray")

all_ps_ord <- decorana(veg = all_ps_blst_vegan_df)

all_ps_ord %>% summary()

all_ps_ord %>% str()

all_ps_ord$cproj
all_ps_ord


plot(all_ps_ord)
plot(all_ps_ord,type = "p")
plot(all_ps_ord,type = "c") 
vegan::scores(all_ps_vg_dist)




# all_ps_blst_vegan_df[,(colnames(all_ps_blst_vegan_df) %in% expected_sps)] %>% colnames()
# all_ps_blst_vegan_df%>% colnames()



all_ps_vegan_ord_meta <- metaMDS(veg = all_ps_blst_vegan_df, comm = all_ps_vg_dist)
# actually autotransform = FALSE doesn't seem to change the results
plot(all_ps_vegan_ord_meta, type = "t")


all_ps_vegan_ord_meta %>% str()
all_ps_vegan_ord_meta$stress


  
#6b- extract NMDS scores from results
  
all_vegan_meta <- (vegan::scores(all_ps_vegan_ord_meta) %>% tidyr::as_tibble(rownames = "Sample number")) %>% mutate(`Sample number` = as.numeric(`Sample number`))
            # all_vegan_meta <- as.data.frame(vegan::scores(all_ps_vegan_ord_meta))
            
            #Using the scores function from vegan to extract the site scores and convert to a data.frame
            
            # all_vegan_meta$`Sample number` <- rownames(all_vegan_meta) %>% as.numeric()  
            
            # all_vegan_meta %>% left_join()# create a column of site names, from the rownames of data.scores
            
            # all_vegan_meta <- all_vegan_meta  %>% as_tibble() # create a column of site names, from the rownames of data.scores

#7- bring NMDS scores to complete table

all_vegan_meta_tbl <- left_join(x = unique(all_ps_blst_vegan[,c(1:10)]),y = all_vegan_meta, by = "Sample number") %>% 
  mutate(Primer=factor(Primer,levels = c("NeoFish", "MiFish", "Teleo")),
         `Mock Community`=factor(`Mock Community`))





library(factoextra)
library(ggforce)



nmds_PLOT <- all_vegan_meta_tbl %>% 
  # filter(Run %in% c("LGC_MiniSeq_1", "LGC_MiniSeq_2")) %>% 
  ggplot(aes(x = NMDS1,y = NMDS2, col = Primer,shape = Normalization,label = Sample,Group = `Mock Community`))+
    # stat_ellipse()+ 
  geom_point(size = 11)+
  theme_light(base_size = 18) +
  theme(legend.position="bottom") +
  coord_fixed(ratio = 1) +
  # ggrepel::geom_label_repel(label.size = 0.8,size = 3,min.segment.length = 2) +
  # ggrepel::geom_text_repel(col="black",size = 3,min.segment.length = 2) +
  # scale_shape_manual() %>% 
  scale_color_manual(
    # labels = c("NeoFish", "MiFish", "Teleo", "NeoFish/MiFish", "NeoFish/MiFish/Teleo"),
    labels = c("NeoFish", "MiFish", "Teleo"),
                     values = alpha(colour = colors_norm[c(1,3,5)] ))+
  annotate(geom = "text",
           x=c(0.275),y=c(-0.275),label=paste0("Stress: ",round(all_ps_vegan_ord_meta$stress,digits = 4)),size=5) +

    # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(inherit.aes = FALSE,
                             aes(x = NMDS1,y = NMDS2,
                                 group=`Mock Community`,
                                 label=`Mock Community`),
                             n = 100,
                             expand = 0.03,
                             label.fontsize = 20,con.cap = 0.1) 
  
    # facet_wrap(~`Mock Community`,ncol = 2)
  
#   
nmds_PLOT
  # 
# ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_NMDS.pdf",
ggsave(file = "~/outros/sfjq_temp/SFJQ_NMDS.pdf",
     plot = nmds_PLOT,
     device = "pdf",
     width = 40,
     height =25,
     units = "cm",
     dpi = 300)
# 
# ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_NMDS.svg",
#      plot = nmds_PLOT,
#      device = "svg",
#      width = 14,
#      height =10,
#      dpi = 600)


ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_NMDS.png",
     plot = nmds_PLOT,
     device = "png",
     width = 31,
     height =20,
     units = "cm",
     dpi = 300)




# anosim----

################################# anosin function###############################
anosin_auto<- function(tbl,cols_out, Coluna){
  df <- tbl[,c(1,(1+cols_out):ncol(tbl))] %>% 
    as.data.frame() %>% 
    `rownames<-`(.$`Sample number`) %>% 
    select(-c("Sample number"))
  
  ano <- anosim(df, grouping = tbl[[Coluna]],
       permutations = 9999, distance = "bray", strata = NULL)
  return(ano)
  
}
################################################################################


#Primers ----
#Todas MC juntas
all_ps_blst_vegan %>% 
  # filter(`Mock Community` %in% c("SFmc")) %>% 
anosin_auto(cols_out = 10,Coluna = "Primer")

#Apenas JQ
all_ps_blst_vegan %>% 
  filter(`Mock Community` %in% c("JQmc")) %>% 
anosin_auto(cols_out = 10,Coluna = "Primer")

#Apenas SF
all_ps_blst_vegan %>% 
  filter(`Mock Community` %in% c("SFmc")) %>% 
anosin_auto(cols_out = 10,Coluna = "Primer")

#Apenas SFJQ
all_ps_blst_vegan %>% 
  filter(`Mock Community` %in% c("SFJQmc")) %>% 
anosin_auto(cols_out = 10,Coluna = "Primer")


#Normalization ----
#Todas MC juntas
all_ps_blst_vegan %>% 
anosin_auto(cols_out = 10,Coluna = "Normalization")

#Apenas JQ
all_ps_blst_vegan %>% 
  filter(`Mock Community` %in% c("JQmc")) %>% 
anosin_auto(cols_out = 10,Coluna = "Normalization")

#Apenas SF
all_ps_blst_vegan %>% 
  filter(`Mock Community` %in% c("SFmc")) %>% 
anosin_auto(cols_out = 10,Coluna = "Normalization")

#Apenas SFJQ
all_ps_blst_vegan %>% 
  filter(`Mock Community` %in% c("SFJQmc")) %>% 
anosin_auto(cols_out = 10,Coluna = "Normalization")


all_vegan_meta_tbl %>% 
  filter(Run %in% c("LGC_MiniSeq_1","LGC_MiniSeq_2")) %>% 
  select(Sample,`Sample number`,Primer)
```

``` r
# all_ps_tbl_blast_bckp6 <- all_ps_tbl_blast

all_ps_tbl_blast <- all_ps_tbl_bl_cur 




#28- ASVs plots by sample and species ----

options(set.seed(seed = 13))

#28a- ASV size distribution - alphabetical----
all_ps_tbl_blast$Run %>% unique()
  
ASV_size_by_Sample <- all_ps_tbl_blast %>%
  filter(!Sample %in% c("Positive Control\n(P.glauca)")) %>% 
  filter(!Run %in% c("ecomol_iSeq")) %>% 
  mutate(Sample = factor(Sample,levels = sample_levels)) %>% 
  mutate(Primer = factor(Primer,levels = c("NeoFish", "MiFish", 
                                           "Teleo", "NeoFish/MiFish/Teleo"))) %>% 
  ggplot(aes(y=Sample,
             x=`ASV size (pb)`,
             colour = Primer,
             size=`Relative abundance on sample`,
             shape=`Expected length`
             )) +
  geom_jitter(height = 0.2,
              width = 0) +
  ggplot2::scale_colour_manual(
                     values = ggplot2::alpha(colour = colors5[1:4] ,alpha =  0.3)) +
  coord_fixed(ratio = 8) +
  scale_x_continuous(breaks = c(20,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340),expand = c(0.02,0.02)) +
  xlab("ASV length (bp)") +
  ylab("Sample") +
  ggtitle(label = "SFJQ mock communities ",
          subtitle = "All ASVs found in samples, by length and abundance") +
  theme_bw(base_size = 15) +
  theme(legend.position = "right") 

ASV_size_by_Sample


ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/3-ASV_size_by_sample.png",
     plot = ASV_size_by_Sample,
     device = "png",
     width = 18,
     height = 10,
     dpi = 600)

ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/3-ASV_size_by_sample.svg",
     plot = ASV_size_by_Sample,
     device = "svg",
     width = 18,
     height = 10,
     dpi = 600)

ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/3-ASV_size_by_sample.pdf",
     plot = ASV_size_by_Sample,
     device = "pdf",
     width = 18,
     height = 10,
     dpi = 600)
dev.off()



##################### paper ################################# ----

#29 - ASVS - sample X species (12S) - distribution: blast---- 
all_ps_tbl_blast$Group %>%unique()

#Jq & SF DNA pools - ASVs size and species by sample----

all_ps_tbl_blast$Sample %>% unfactor() %>% unique() %>% base::sort()
all_ps_tbl_blast$`revised final ID` %>% unique() %>% base::sort() %>% paste0(collapse = '", \n"') %>% cat()
all_ps_tbl_blast$`final ID` %>% unique() %>% base::sort() %>% paste0(collapse = '", \n"') %>% cat()


pool_labels <- c(
  "Normalized JQmc" = "Jequitinhonha\nNormalized DNA pool",
  "Non-normalized JQmc" = "Jequitinhonha\nNon-Normalized DNA pool",
  "Normalized SFJQmc" = "São Francisco & Jequitinhonha\nNormalized DNA pool",
  "Non-normalized SFmc" = "São Francisco\nNon-Normalized DNA pool",
  "Normalized SFmc" = "São Francisco\nNormalized DNA pool"
)

pools_levels <- c(
"Normalized JQmc",
"Non-normalized JQmc",
"Normalized SFJQmc",
"Normalized SFmc",
"Non-normalized SFmc"
)

#final ID levels
{
finalID_levels<- c(
#jq
"Astyanax lacustris", 
"Australoheros sp", 
"Delturus brevis", 
"Eugerres brasilianus", 
"Hypomasticus steindachneri", 

"Hypostomus nigrolineatus", 
"Megaleporinus elongatus", 
"Megaleporinus garmani", 
"Moenkhausia costae", 
"Rhamdia quelen", 
"Steindachneridion amblyurum", 
"Wertheimeria maculata",
#ambos
"Gymnotus carapo", 
"Hoplias", 
"Hoplias brasiliensis", 

"Hoplias intermedius", 
"Hoplias malabaricus", 
"Prochilodus", 
"Prochilodus argenteus", 
 
"Prochilodus costatus", 
"Steindachnerina elegans", 
"Trachelyopterus galeatus",
#sf

"Astyanax fasciatus", 
"Brycon orthotaenia", 
"Characidium lagosantense", 
"Crenicichla lepidota", 
# "Curimatella lepidura", esse é na vdd roeboides
"Roeboides xenodon",
"Eigenmannia virescens", 
"Franciscodoras marmoratus", 
"Hypostomus alatus", 
"Imparfinis minutus", 
"Leporinus reinhardti", 
"Microglanis leptostriatus", 
"Moenkhausia sanctaefilomenae", 
"Myleus micans", 
"Pamphorichthys hollandi", 
"Phalloceros uai", 

"Pimelodus maculatus", 
"Pimelodus pohli", 
"Pseudoplatystoma corruscans", 
"Pterygoplichthys etentaculatus", 
"Serrasalmus brandtii", 
"Tetragonopterus chalceus", 

#partial
"Astyanax", 
"Hoplias brasiliensis/intermedius",
"Hypostomus", 
"Pimelodus", 
"Prochilodus argenteus/hartii",
#trash


"NA elegans/gilbert", 
"NA lepidura/xenodon", 
"Acestrorhynchus lacustris",
"Coptodon zillii",
"Cyphocharax gilbert",
"Geophagus brasiliensis",
"Planaltina myersi", 
"Poecilia reticulata mitochondr"
)
}
finalID_levels[finalID_levels %>% duplicated()]



# final ID levels 2
{
finalID_levels <- c(
"Acestrorhynchus lacustris", 
"Acinocheirodon melanogramma", 
"Astyanax fasciatus", 
"Astyanax lacustris", 
"Australoheros sp", 
"Bos taurus", 
"Brycon orthotaenia", 
"Characidium lagosantense", 
"Coptodon zillii", 
"Crenicichla lepidota", 
# "Curimatella lepidura", agora é roeboides
"Roeboides xenodon",
"Cyphocharax gilbert", 
"Delturus brevis", 
"Eigenmannia virescens", 
"Eugerres brasilianus", 
"Franciscodoras marmoratus", 
"Geophagus brasiliensis", 
"Gymnotus carapo", 
"Hoplias brasiliensis", 
"Hoplias intermedius", 
"Hoplias malabaricus", 
"Hypomasticus steindachneri", 
"Hypostomus alatus", 
"Hypostomus nigrolineatus", 
"Imparfinis minutus", 
"Leporinus reinhardti", 
"Megaleporinus elongatus", 
"Megaleporinus garmani", 
"Microglanis leptostriatus", 
"Moenkhausia costae", 
"Moenkhausia sanctaefilomenae", 
"Myleus micans", 
"Pamphorichthys hollandi", 
"Phalloceros uai", 
"Pimelodus maculatus", 
"Pimelodus pohli", 
"Planaltina myersi", 
"Prionace glauca", 
"Prochilodus argenteus", 
"Prochilodus costatus", 
"Prochilodus hartii",
"Pseudoplatystoma corruscans", 
"Pterygoplichthys etentaculatus", 
"Rhamdia quelen", 
"Roeboides xenodon", 
"Serrasalmus brandtii", 
"Steindachneridion amblyurum", 
"Tetragonopterus chalceus", 
"Trachelyopterus galeatus", 
"Wertheimeria maculata",

#partial
"Astyanax", 
"Hoplias brasiliensis/intermedius", 
"Hypostomus", 
"Pimelodus", 
"Prochilodus", 
"Prochilodus argenteus/hartii")
}

sps_remove <- c(
NA,"NA",
"Acestrorhynchus lacustris", 
"Acinocheirodon melanogramma", 
"Bos taurus", 
"Coptodon zillii", 
"Curimatella lepidura", 
"Eugerres brasilianus", 
"Geophagus brasiliensis", 
"Leporinus reinhardti", 
"Moenkhausia costae", 
"Planaltina myersi", 
"Prionace glauca", 
"Pseudoplatystoma corruscans"
)
```

### Count ASVs

Now we will count total ASVs and which were considered *in range* or not
for each sample

# Images arcticle

## Fold-change bar plots

\#\#Ven diagram

\#\#Upset plot

\#\#gráficos de correlação DNA/RRA

``` r
colnames(all_ps_tbl_sfjq_full_uniq)
all_ps_tbl_sfjq_full_uniq$Group
all_ps_tbl_sfjq_full_uniq$MC
all_ps_tbl_sfjq_full_uniq$Normalization




######função pra plotar lm

# lm_eqn <- function(df){
#     m <- lm(y ~ x, df);
#     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#          list(a = format(unname(coef(m)[1]), digits = 2),
#               b = format(unname(coef(m)[2]), digits = 2),
#              r2 = format(summary(m)$r.squared, digits = 3)))
#     as.character(as.expression(eq));
# }





#inclinação das retas ----

#LM
 models_corr <- all_ps_tbl_sfjq_full_uniq %>% 
   filter(Normalization %in% c("Non-normalized")) %>% 
   # filter(Primer %in% c("NeoFish")) %>% 
  group_by(Primer,MC) %>% 
   select(c("Percentage on respective Non-norm pool", "Relative abundance on sample")) %>% 
  do(model = lm(formula = `Percentage on respective Non-norm pool` ~ `Relative abundance on sample`, data = .))
   # lm(formula = `Percentage on respective Non-norm pool` ~ `Relative abundance on sample`)
   # aov(formula = `Percentage on respective Non-norm pool` ~ `Relative abundance on sample`)

models_corr$model
models_corr$Primer
models_corr$MC



models_corr$model[[1]] %>% broom::glance(model) #NeoFish JQmc
models_corr$model[[2]] %>% broom::glance(model) #NeoFish SFmc
models_corr$model[[3]] %>% broom::glance(model) #MiFish JQmc
models_corr$model[[4]] %>% broom::glance(model) #MiFish SFmc
models_corr$model[[5]] %>% broom::glance(model) #Teleo JQmc

#AOV
models_corr <- all_ps_tbl_sfjq_full_uniq %>% 
   filter(Normalization %in% c("Non-normalized")) %>% 
   # filter(Primer %in% c("NeoFish")) %>% 
  group_by(Primer) %>% 
   select(c("Percentage on respective Non-norm pool", "Relative abundance on sample")) %>% 
  do(model = aov(formula = `Percentage on respective Non-norm pool` ~ `Relative abundance on sample`, data = .))
   # lm(formula = `Percentage on respective Non-norm pool` ~ `Relative abundance on sample`)
   # aov(formula = `Percentage on respective Non-norm pool` ~ `Relative abundance on sample`)

models_corr$model
models_corr$Primer

models_corr$model[[1]] %>% broom::tidy()
models_corr$model[[2]] %>% broom::tidy()
models_corr$model[[3]] %>% broom::tidy()



#JQmc ----
# jq_coef <- lm(`Percentage on respective Non-norm pool` ~ `Relative abundance on sample`,
#    (all_ps_tbl_sfjq_full_uniq %>% 
#       filter(MC %in% c("JQmc")) %>% 
#       filter(Normalization %in% c("Non-normalized"))))




############################ tentando com a tib com norm e n norm na mesma coluna
tab_curated_SFJQ_all_pools %>% colnames()


sfjq_sp_corr <-
  tab_curated_SFJQ_all_pools %>% 
  ggplot(aes(x=`input DNA (%)`,
             y=`RRA (%)`,
             col=Primer,
             shape=Pool))+
  geom_point() + 
  # geom_smooth(method=lm) +
  # coord_fixed()+
  scale_colour_manual(values = alpha(colour = colors_norm[c(1,3,5)] ,alpha =  0.8)) +
  # scale_shape_manual(drop=TRUE) +
  xlab("Input DNA (%)") +
  ylab("Relative Read Abundance (%)") +
  ggtitle("SFmc & JQmc: Correlation between\nInput DNA and RRA") +
  scale_x_sqrt(breaks=c(0,0.0001,0.001,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1)*100) +
  scale_y_sqrt(breaks=c(0,0.0001,0.001,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1)*100) +
  coord_fixed(ratio = 1)+
  geom_smooth(method=lm) +
  theme_bw(base_size = 10) +
  # facet_wrap(MC~Primer,ncol = 3) 
  facet_wrap(Primer~Group,ncol = 5) 
# +
#   scale_x_log10() +
#   scale_y_log10()
# minor_breaks = mb
# minor_breaks = mb

ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_RRA_DNA_bySp2.png", plot = sfjq_sp_corr, device = "png", width = 30, height = 15, units = "cm", dpi = 600)









###########################




sfjq_sp_corr <- all_ps_tbl_sfjq_full_uniq %>% 
  # filter(MC %in% c("JQmc")) %>%
  filter(Normalization %in% c("Non-normalized")) %>%
  # mutate(`Percentage on respective Non-norm pool`= if_else(`Percentage on respective Non-norm pool` %in% c(NA,"NA"),0,`Percentage on respective Non-norm pool`)) %>% 
  # mutate(`Relative abundance on sample`= if_else(`Relative abundance on sample` %in% c(NA,"NA"),0,`Relative abundance on sample`)) %>% View()
  
  ggplot(aes(x=`Percentage on respective Non-norm pool`*100,
             y=`Relative abundance on sample`*100,
             col=Primer,
             shape=MC
             ))+
  geom_point() + 
  # geom_smooth(method=lm) +
  # coord_fixed()+
  scale_colour_manual(values = alpha(colour = colors_norm[c(1,3,5)] ,alpha =  0.8)) +
  # scale_shape_manual(drop=TRUE) +
  xlab("Input DNA (%)") +
  ylab("Relative Read Abundance (%)") +
  ggtitle("SFmc & JQmc: Correlation between\nInput DNA and RRA") +
  scale_x_sqrt(breaks=c(0,0.0001,0.001,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1)*100) +
  scale_y_sqrt(breaks=c(0,0.0001,0.001,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1)*100) +
  coord_fixed(ratio = 1)+
  geom_smooth(method=lm) +
  theme_bw(base_size = 10) +
  # facet_wrap(MC~Primer,ncol = 3) 
  facet_wrap(Primer~Normalization,ncol = 3) 
# +
#   scale_x_log10() +
#   scale_y_log10()
# minor_breaks = mb
# minor_breaks = mb

ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_RRA_DNA_bySp.png", plot = sfjq_sp_corr, device = "png", width = 30, height = 15, units = "cm", dpi = 600)




#SFmc ----


library("plotly")

sf_sp_corr <- all_ps_tbl_sfjq_full_uniq %>% 
  filter(MC %in% c("SFmc")) %>% 
  filter(Normalization %in% c("Non-normalized")) %>% 
  # mutate(`Percentage on respective Non-norm pool`= if_else(`Percentage on respective Non-norm pool` %in% c(NA,"NA"),0,`Percentage on respective Non-norm pool`)) %>% 
  # mutate(`Relative abundance on sample`= if_else(`Relative abundance on sample` %in% c(NA,"NA"),0,`Relative abundance on sample`)) %>% View()
  
  ggplot(aes(x=`Percentage on respective Non-norm pool`*100,
             y=`Relative abundance on sample`*100,
             col=Primer,
             shape=MC))+
  geom_point() + 
  # geom_smooth(method=lm) +
  # coord_fixed()+
  scale_colour_manual(values = alpha(colour = colors_norm[c(1,3,5)] ,alpha =  0.8)) +
  xlab("input DNA (%)") +
  ylab("Relative Read Abundance (%)") +
  ggtitle("SFmc: Correlation between\nInput DNA and RRA") +
  scale_x_sqrt(breaks=c(0,0.0001,0.001,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1)*100) +
  scale_y_sqrt(breaks=c(0,0.0001,0.001,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1)*100) +
  coord_fixed(ratio = 1)+
  geom_smooth(method=lm) +
  theme_bw(base_size = 10) +
  facet_wrap(~Primer,ncol = 3) 
    
    
    
    ggplotly(sf_sp_corr)


ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SF_RRA_DNA_bySp.png", plot = sf_sp_corr, device = "png", width = 20, height = 15, units = "cm", dpi = 600)
```

## Gráficos Daniel

``` r
# set colors ----

colors6 <- c("#19821c","#8dbc8f", #neo norm & non norm
             "#171d9a","#8d90c7", #mif norm & non norm
             "#8c1717","#c18d8d"  #tel norm & non norm
             )

scales::show_col(colors6)






# carregando tabela 
# https://docs.google.com/spreadsheets/d/1sTsOaI999E9Py_f4ll5j838r7CVc991a8mPmd8LvykU/edit#gid=0

tab_curated_SFJQ_all_pools <- read.csv("/home/heron/prjcts/fish_eDNA/sfjq/data/SFJQ_all_pools-tabelas_curadas_para_o_R-SFJQmc-copy.csv",check.names = F) %>% 
  as_tibble()

tab_curated_SFJQ_all_pools <- tab_curated_SFJQ_all_pools %>% 
    filter(`input DNA (%)` != 0 & `RRA (%)` !=0) %>%
  group_by(Pool,Normalization,Primer,Species) %>% 
  summarize(Pool = unique(Pool),
            Normalization = unique(Normalization),
            Status = unique(Status),
            Species = unique(Species),
            Primer = unique(Primer),
            `Num ASVs` = `Num ASVs`,
            `Num OTUs` = `Num OTUs`,
            `input DNA (%)` = `input DNA (%)`,
            `RRA (%)` = `RRA (%)`
            # `Expected species` = length(unique(`revised final ID`[`Expected Species` %in% c("expected")])),
            # `Expected species list` = list(unique(base::sort(`revised final ID`[`Expected Species` %in% c("expected")]))),
            # `revised final ID`= unique(`revised final ID`),
            # `RRA (%)` = sum(`Relative abundance on sample`),
            # `Percentage on respective Norm pool` = unique(`Percentage on respective Norm pool`),
            # `Percentage on respective Non-norm pool` = unique(`Percentage on respective Non-norm pool`),
            # `Relative abundance on sample` = unique(`Relative abundance on sample`)
            ) %>%
  ungroup() %>% 
  # mutate(`revised final ID`=factor(`revised final ID`
  #                                  # , levels = rev(finalID_levels)
  #                                  )) %>% 
  mutate(Primer = factor(Primer, levels = c("NeoFish", "MiFish", "Teleo", "NeoFish/MiFish","NeoFish/MiFish/Teleo"))) %>% 
  mutate(`Fold change (RRA/DNA input)` = gtools::foldchange(denom = `input DNA (%)`,num = `RRA (%)`)) %>% 
  unique()
  # mutate(Group = factor(Group,levels = c("Normalized SFJQmc", "Non-normalized JQmc", "Normalized JQmc", "Non-normalized SFmc", "Normalized SFmc" ))) %>% 
  # mutate(MC = str_remove(Group,"Normalized |Non-normalized ")) %>% 
  # mutate(Normalization = str_remove(Group," SFmc| JQmc| SFJQmc")) 


tab_curated_SFJQ_all_pools$`Fold change (RRA/DNA input)` %>% sort() %>% duplicated()

tab_curated_SFJQ_all_pools[c(which(tab_curated_SFJQ_all_pools$`Fold change (RRA/DNA input)` %>% as.character() %>% duplicated()),
                             which(tab_curated_SFJQ_all_pools$`Fold change (RRA/DNA input)` %>% as.character() %>% duplicated())-1),] %>% View()




# SFJQmc ----
## build tree ----
### read 12S db seqs for the species present in pools
SFJQ_Sps_seqs <- Biostrings::readDNAStringSet(filepath = "~/prjcts/fish_eDNA/sfjq/data/trees/SFJQmc-38_SPs_fused-unique.fas") %>%
# SFJQ_Sps_seqs <- Biostrings::readDNAStringSet(filepath = "~/outros/sfjq_temp/trees/SFJQmc-38_SPs_fused-unique.fas") %>% 
  DECIPHER::RemoveGaps()

### align seqs
SFJQ_Sps_algn <- DECIPHER::AlignSeqs(myXStringSet = SFJQ_Sps_seqs, 
                                      refinements = 100,
                                      iterations = 100,
                                      verbose = TRUE)

### generate distance matrix
SFJQ_Sps_dist <- DECIPHER::DistanceMatrix(myXStringSet = SFJQ_Sps_algn,
                                            includeTerminalGaps = FALSE,
                                            correction = "Jukes-Cantor",
                                            processors = 20,
                                            verbose = TRUE)

### generate dendrogram/tree from alignment and distance matrix
SFJQmc_tree <- ape::nj(SFJQ_Sps_dist)
        # tree <- phangorn::NJ(SFJQ_Sps_dist)
class(SFJQmc_tree)

### save tree as newick
ape::write.tree(phy = SFJQmc_tree,file = "~/prjcts/fish_eDNA/sfjq/data/trees/SFJQmc-38_SPs_fused-unique_APEtree.nwk")
# ape::write.tree(phy = SFJQmc_tree,file = "~/outros/sfjq_temp/trees/SFJQmc-38_SPs_fused-unique_APEtree.nwk")

### read tree from file (or stay with the same object)
# SFJQmc_tree <- read.tree("~/outros/sfjq_temp/trees/SFJQmc-38_SPs_fused-unique_APEtree.nwk")
SFJQmc_tree <- read.tree("~/prjcts/fish_eDNA/sfjq/data/trees/SFJQmc-38_SPs_fused-unique_APEtree.nwk")


## species metadata ----
### read table with pools species, input DNA and RRA
# tab_curated_SFJQ_all_pools <- read.csv("~/outros/sfjq_temp/SFJQ_all_pools-tabelas_curadas_para_o_R-SFJQmc.csv",check.names = F)
# tab_curated_SFJQ_all_pools <- read.csv("/home/heron/prjcts/fish_eDNA/sfjq/data/SFJQ_all_pools-tabelas_curadas_para_o_R-SFJQmc.csv",check.names = F) %>% 
#   as_tibble()
# 
# tab_curated_SFJQ_all_pools




tab_curated_SFJQ_all_pools %>% colnames()
tab_curated_SFJQ_all_pools %>% unique()
all_ps_tbl_sfjq_full_uniq %>% colnames()










# tab_curated_SFJQ_all_pools  <- all_ps_tbl_sfjq_full_uniq



#converting the corrected table to the format required for next steps


tab_curated_SFJQ_all_pools <-
  all_ps_tbl_sfjq_full_uniq %>%
  pivot_longer(c("Percentage on respective Norm pool", "Percentage on respective Non-norm pool"),
               names_to = "DNA Originary Pool", values_to = "input DNA (%)") %>% 
  pivot_longer(c("Recovered proportion norm", "Recovered proportion non norm"),
               names_to = "Fold Change Originary Pool", values_to = "Fold Change") %>% 
  select(-c(
    # "DNA Originary Pool", "RRA Originary Pool",
            # "Sample,"
            # "Group",
            "Run"
            )) %>% 
  # colnames()
  # View()
  rename(
   Species = `revised final ID`,
  `Num OTUs` = `OTUs`,
  `Num ASVs` = `ASVs`,
  `Num IDs` = `IDs`,
  # `RRA (%)` = `Relative abundance on sample`,
  # `Expected species`, 
  # `Expected species list`,
  # `Fold Change` = `Ratio`,
  `Pool` = `MC`) %>%
  filter((`DNA Originary Pool` %in% c("Percentage on respective Norm pool") & `Normalization` %in% c("Normalized")) |
           (`DNA Originary Pool` %in% c("Percentage on respective Non-norm pool") & `Normalization` %in% c("Non-normalized")) ) %>% 
  filter((`Fold Change Originary Pool` %in% c("Recovered proportion norm") & `Normalization` %in% c("Normalized")) |
           (`Fold Change Originary Pool` %in% c("Recovered proportion non norm") & `Normalization` %in% c("Non-normalized")) ) %>% 
  select(-c("DNA Originary Pool", "Fold Change Originary Pool")) %>% 
  mutate(Status = if_else(`Expected species` == 0, "Contamination","Expected")) %>% 
  mutate(`RRA (%)` = `RRA (%)` * 100,
         `input DNA (%)` = `input DNA (%)` * 100) %>% unite(Normalization, Primer, col = "Primer_norm",remove = F,sep = " ") %>% 
  mutate(Primer_norm = factor(Primer_norm, levels = c("Normalized NeoFish","Normalized MiFish","Normalized Teleo","Non-normalized NeoFish","Non-normalized MiFish","Non-normalized Teleo"))) %>% View()

# all_ps_tbl_sfjq_full_uniq %>% colnames()


 

# tab_curated_SFJQ_all_pools[sort(c(which(tab_curated_SFJQ_all_pools$`Fold Change` %>% as.character() %>% duplicated()),which(tab_curated_SFJQ_all_pools$`Fold Change`%>% as.character() %>% duplicated())-1)),] %>% View()
# tab_curated_SFJQ_all_pools[tab_curated_SFJQ_all_pools$`Fold Change` %>% duplicated(),] %>% View()

### tidy table
tab_curated_SFJQ <- tab_curated_SFJQ_all_pools %>%
  as_tibble() %>%
  filter(Status == "Expected") %>%
  filter(Pool == "SFJQmc") %>%
  # filter(MC == "SFJQmc") %>% 
  # filter(Pool == "SFJQmc") %>% 
  mutate(Primer = factor(Primer, levels = c("NeoFish","MiFish","Teleo")),
         Species = factor(Species),
         # `revised final ID` = factor(`revised final ID`),
         Normalization = factor(Normalization),
         Pool = factor(Pool),
         Status = factor(Status),
         # `RRA (%)` = as.numeric(`Relative abundance on sample`)*100,
         `RRA (%)` = as.numeric(`RRA (%)`),
         `input DNA (%)` = as.numeric(`input DNA (%)`))


## correct tips labels ----
### check tip names in tree
SFJQmc_tree$tip.label


### rename tree tips to mach table
SFJQmc_tree$tip.label
tab_curated_SFJQ$Species %>% unique() %>% sort()

{
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1339_Prochilodus_costatus_JQ_2860"] <- "Prochilodus costatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_5612_Prochilodus_argenteus_hartii"] <- "Prochilodus argenteus/hartii"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_6586_Steindachnerina_elegans-Cyphocharax_gilbert"] <- "Cyphocharax gilbert"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1230_Serrasalmus_brandtii"] <- "Serrasalmus brandtii"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1162_Myleus_micans"] <- "Myleus micans"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_2984_Hypomasticus_steindachneri"] <- "Hypomasticus steindachneri"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_4290_Megaleporinus_garmani"] <- "Megaleporinus garmani"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_1762_Megaleporinus_elongatus"] <- "Megaleporinus elongatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1037_Brycon_orthotaenia"] <- "Brycon orthotaenia"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1381_Characidium_lagosantense"] <- "Characidium lagosantense"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "1136_Acestrorhynchus_lacustris"] <- "Acestrorhynchus lacustris"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_7822_Hoplias_malabaricus"] <- "Hoplias malabaricus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1377_Hoplias_intermedius_brasiliensis"] <- "Hoplias brasiliensis/intermedius"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_0916_Roeboides_xenodon"] <- "Roeboides xenodon"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_0901_Tetragonopterus_chalceus"] <- "Tetragonopterus chalceus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1113_Moenkhausia_sanctaefilomenae"] <- "Moenkhausia sanctaefilomenae"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SFJQ_7812_Moenkhausia_costae"] <- "Moenkhausia costae"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1153_Astyanax_cf_fasciatus"] <- "Astyanax fasciatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_7880_Astyanax_lacustris"] <- "Astyanax lacustris"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1141_Pterygoplichthys_etentaculatus"] <- "Pterygoplichthys etentaculatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1041_Hypostomus_alatus"] <- "Hypostomus alatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_7893_Hypostomus_nigrolineatus"] <- "Hypostomus nigrolineatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1248_Gymnotus_carapo_JQ_1631"] <- "Gymnotus carapo"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1117_Eigenmannia_virescens"] <- "Eigenmannia virescens"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_1936_Delturus_brevis"] <- "Delturus brevis"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_0708_Franciscodoras_marmoratus"] <- "Franciscodoras marmoratus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_7817_Wertheimeria_maculata"] <- "Wertheimeria maculata"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1243_Trachelyopterus_galeatus_JQ_5675"] <- "Trachelyopterus galeatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1403_Pimelodus_pohli"] <- "Pimelodus pohli"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1280_Pimelodus_maculatus"] <- "Pimelodus maculatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "1135_Pseudoplatystoma_corruscans"] <- "Pseudoplatystoma corruscans"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_2996_Steindachneridion_amblyurum"] <- "Steindachneridion amblyurum"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1104_Microglanis_leptostriatus"] <- "Microglanis leptostriatus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1161_Imparfinis_minutus"] <- "Imparfinis minutus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_5645_Rhamdia_aff_quelen"] <- "Rhamdia quelen"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SFJQ_2943_Eugerres_brasilianus"] <- "Eugerres brasilianus"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1368_Phalloceros_uai"] <- "Phalloceros uai"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1099_Pamphorichthys_hollandi"] <- "Pamphorichthys hollandi"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "SF_1264_Crenicichla_lepidota"] <- "Crenicichla lepidota"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "JQ_6584_Australoheros_sp"] <- "Australoheros sp"
SFJQmc_tree$tip.label[SFJQmc_tree$tip.label == "304_Geophagus_brasiliensis"] <- "Geophagus brasiliensis"
}

### check if all names have correspondences
SFJQmc_tree$tip.label %in% tab_curated_SFJQ$Species
tab_curated_SFJQ$Species %in% SFJQmc_tree$tip.label
tab_curated_SFJQ$Species[!tab_curated_SFJQ$Species %in% SFJQmc_tree$tip.label]





### convert ape tree to prettier ggtree object
SFJQmc_tree_plot <- ggtree(SFJQmc_tree) + 
  theme_tree2() +
  geom_tiplab(offset = 0,align = T)+ 
  xlim(0, 0.42)


### save tree plot
# ggsave(file = "~/prjcts/fish_eDNA/sfjq/data/trees/SFJQmc-38_SPs_fused-unique-tree_plot.pdf", plot = SFJQmc_tree_plot, device = "pdf", width = 24, height = 24, units = "cm", dpi = 600)
ggsave(file = "~/outros/sfjq_temp/trees/SFJQmc-38_SPs_fused-unique-tree_plot.pdf", plot = SFJQmc_tree_plot, device = "pdf", width = 24, height = 24, units = "cm", dpi = 600)



### extract tips order from tree to reproduce on plots
SPs_order_in_SFJQ_tree <- ggtree::get_taxa_name(SFJQmc_tree_plot) %>% rev()

            # SPs_order_in_SFJQ_tree <- extract_tree_data(tree_plot) %>% 
            #     dplyr::filter(isTip) %>% 
            #     dplyr::pull(label)



## plots ----
### RRA histogram ----

# library(ggdist)



SFJQmc_RRA_plot <-
  tab_curated_SFJQ %>%
  filter(Pool %in% c("SFJQmc")) %>% 
  mutate(Species = factor(Species,levels = SPs_order_in_SFJQ_tree)) %>% 
  ggplot(aes(y = Species,
             x = `RRA (%)`,
             fill = Primer, 
             col = Normalization
             ))+
  geom_bar(stat = "identity", size = 0.3,width = .75,
           # alpha = 0.7,
           position = position_dodge(preserve = "single" ,width = 1.1)) +
           # position = position_dodgejust(preserve = "single" ,width = 1.2)) +
  scale_color_manual(values = c("#000000","#747474")) +
  scale_fill_manual(values = colors6[c(1,3)]) +
  geom_point(aes(y = Species,
                 x = `input DNA (%)`),
             shape = "|",
             size = 3,
             colour = "#000000") +
  scale_x_break(c(13, 30),
                scales = "fixed") +
  scale_x_continuous(breaks=c(0,5,10,13,30,32)) +
  
  xlab("Relative read abundance (%)")+ 
  # opts(axis.title.y = theme_text(vjust=-0.5))
 theme(axis.text = element_text(vjust = -0.5))

## save plot 

# ggsave(file = "~/prjcts/fish_eDNA/sfjq/data/trees/SFJQmc-38_SPs_fused-unique-RRA_barplot.pdf", plot = SFJQmc_RRA_plot, device = "pdf", width = 18, height = 24, units = "cm", dpi = 600)
dev.off()
ggsave(file = "~/outros/sfjq_temp/trees/SFJQmc-38_SPs_fused-unique-RRA_barplot.pdf", plot = SFJQmc_RRA_plot, device = "pdf", width = 30, height = 24, units = "cm", dpi = 600)




# SFmc ----
## build tree ----
### read 12S db seqs for the species present in pools and select only pool species
names(SFJQ_Sps_seqs) %>% sort() %>% paste0(collapse = '",\n"') %>% cat()

SF_Sps_seqs <- SFJQ_Sps_seqs[c("SF_0708_Franciscodoras_marmoratus", "SF_0901_Tetragonopterus_chalceus",
                               "SF_0916_Roeboides_xenodon", "SF_1037_Brycon_orthotaenia",
                               "SF_1041_Hypostomus_alatus", "SF_1099_Pamphorichthys_hollandi",
                               "SF_1104_Microglanis_leptostriatus", "SF_1113_Moenkhausia_sanctaefilomenae",
                               "SF_1117_Eigenmannia_virescens", "SF_1141_Pterygoplichthys_etentaculatus",
                               "SF_1153_Astyanax_cf_fasciatus", "SF_1161_Imparfinis_minutus",
                               "SF_1162_Myleus_micans", "SF_1230_Serrasalmus_brandtii",
                               "SF_1243_Trachelyopterus_galeatus_JQ_5675", "SF_1248_Gymnotus_carapo_JQ_1631",
                               "SF_1264_Crenicichla_lepidota", "SF_1280_Pimelodus_maculatus",
                               "SF_1339_Prochilodus_costatus_JQ_2860", "SF_1368_Phalloceros_uai",
                               "SF_1377_Hoplias_intermedius_brasiliensis", "SF_1381_Characidium_lagosantense",
                               "SF_1403_Pimelodus_pohli")]
### align seqs
SF_Sps_algn <- DECIPHER::AlignSeqs(myXStringSet = SF_Sps_seqs, 
                                      refinements = 100,
                                      iterations = 100,
                                      verbose = TRUE)

### generate distance matrix
SF_Sps_dist <- DECIPHER::DistanceMatrix(myXStringSet = SF_Sps_algn,
                                            includeTerminalGaps = FALSE,
                                            correction = "Jukes-Cantor",
                                            processors = 20,
                                            verbose = TRUE)

### generate dendrogram/tree from alignment and distance matrix
SFmc_tree <- ape::nj(SF_Sps_dist)
        # tree <- phangorn::NJ(SFJQ_Sps_dist)
class(SFmc_tree)

### save tree as newick
# ape::write.tree(phy = SFmc_tree,file = "~/prjcts/fish_eDNA/sfjq/data/trees/SFmc-23_SPs_fused-unique_APEtree.nwk")
ape::write.tree(phy = SFmc_tree,file = "~/outros/sfjq_temp/trees/SFmc-23_SPs_fused-unique_APEtree.nwk")

### read tree from file (or stay with the same object)
SFmc_tree <- read.tree("~/prjcts/fish_eDNA/sfjq/data/trees/SFmc-23_SPs_fused-unique_APEtree.nwk")
SFmc_tree <- read.tree("~/outros/sfjq_temp/trees/SFmc-23_SPs_fused-unique_APEtree.nwk")

## species metadata ----
### read table with pools species, input DNA and RRA

### tidy table
tab_curated_SF <- tab_curated_SFJQ_all_pools %>% 
  filter(Pool %in% c("SFmc")) %>% 
  filter(Status %in% c("Expected")) %>%
  # filter(MC == "SFJQmc") %>% 
  # filter(Pool == "SFJQmc") %>% 
  mutate(Primer = factor(Primer, levels = c("NeoFish","MiFish","Teleo")),
         Species = factor(Species),
         # `revised final ID` = factor(`revised final ID`),
         Normalization = factor(Normalization),
         Pool = factor(Pool),
         Status = factor(Status),
         # `RRA (%)` = as.numeric(`Relative abundance on sample`)*100,
         `RRA (%)` = as.numeric(`RRA (%)`),
         `input DNA (%)` = as.numeric(`input DNA (%)`))


## correct tips labels ----
### check tip names in tree
SFmc_tree$tip.label


### rename tree tips to mach table
SFmc_tree$tip.label
tab_curated_SF$Species %>% unique() %>% sort()

{
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1339_Prochilodus_costatus_JQ_2860"] <- "Prochilodus costatus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1230_Serrasalmus_brandtii"] <- "Serrasalmus brandtii"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1162_Myleus_micans"] <- "Myleus micans"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1037_Brycon_orthotaenia"] <- "Brycon orthotaenia"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1381_Characidium_lagosantense"] <- "Characidium lagosantense"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1377_Hoplias_intermedius_brasiliensis"] <- "Hoplias intermedius"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_0916_Roeboides_xenodon"] <- "Roeboides xenodon"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_0901_Tetragonopterus_chalceus"] <- "Tetragonopterus chalceus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1113_Moenkhausia_sanctaefilomenae"] <- "Moenkhausia sanctaefilomenae"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1153_Astyanax_cf_fasciatus"] <- "Astyanax fasciatus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1141_Pterygoplichthys_etentaculatus"] <- "Pterygoplichthys etentaculatus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1041_Hypostomus_alatus"] <- "Hypostomus alatus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1248_Gymnotus_carapo_JQ_1631"] <- "Gymnotus carapo"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1117_Eigenmannia_virescens"] <- "Eigenmannia virescens"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_0708_Franciscodoras_marmoratus"] <- "Franciscodoras marmoratus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1243_Trachelyopterus_galeatus_JQ_5675"] <- "Trachelyopterus galeatus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1403_Pimelodus_pohli"] <- "Pimelodus pohli"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1280_Pimelodus_maculatus"] <- "Pimelodus maculatus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1104_Microglanis_leptostriatus"] <- "Microglanis leptostriatus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1161_Imparfinis_minutus"] <- "Imparfinis minutus"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1368_Phalloceros_uai"] <- "Phalloceros uai"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1099_Pamphorichthys_hollandi"] <- "Pamphorichthys hollandi"
SFmc_tree$tip.label[SFmc_tree$tip.label == "SF_1264_Crenicichla_lepidota"] <- "Crenicichla lepidota"
}

### check if all names have correspondences
SFmc_tree$tip.label %in% tab_curated_SF$Species
tab_curated_SF$Species %in% SFmc_tree$tip.label
tab_curated_SF$Species[!tab_curated_SF$Species %in% SFmc_tree$tip.label]





### convert ape tree to prettier ggtree object
SFmc_tree_plot <- ggtree(SFmc_tree) + 
  theme_tree2() +
  geom_tiplab(offset = 0,align = T)+ 
  xlim(0, 0.42)


### save tree plot
# ggsave(file = "~/prjcts/fish_eDNA/sfjq/data/trees/SFmc-23_SPs_fused-unique-tree_plot.pdf", 
ggsave(file = "~/outros/sfjq_temp/trees/SFmc-23_SPs_fused-unique-tree_plot.pdf", 
       plot = SFmc_tree_plot, 
       device = "pdf", 
       width = 24, height = 20, 
       units = "cm", dpi = 600)



### extract tips order from tree to reproduce on plots
SPs_order_in_SF_tree <- get_taxa_name(SFmc_tree_plot) %>% rev()

            # SPs_order_in_SF_tree <- extract_tree_data(tree_plot) %>% 
            #     dplyr::filter(isTip) %>% 
            #     dplyr::pull(label)




## plots ----
### RRA histogram ----

SFmc_RRA_plot <-
  tab_curated_SF %>%
  # filter(Normalization %in% c("Normalized")) %>%
  unite(Normalization, Primer, col = "Primer_norm",remove = F,sep = " ") %>%
  mutate(Primer_norm = factor(Primer_norm, levels = c("Normalized NeoFish","Normalized MiFish","Non-normalized NeoFish","Non-normalized MiFish"))) %>%
  mutate(Species = factor(Species,levels = SPs_order_in_SF_tree)) %>% 
  ggplot(aes(y = Species,
             x = `RRA (%)`,
             # fill = Primer, 
             fill = Primer_norm, 
             col = Normalization
             ))+
  geom_bar(stat = "identity", size = 0.3,width = 1,
           # alpha = 0.7,
           position = position_dodge(preserve = "single" ,width = 1.5)) +
           # position = position_dodgejust(preserve = "single" ,width = 1.2)) +
    scale_fill_manual(values = colors6[c(1,3,2,4)], name = "") +
  geom_point(aes(y = Species,
                 x = `input DNA (%)`,
             colour = Normalization),
             shape = "|",
             size = 3) +
    scale_color_manual(values = c("#747474","#000000")) +
  xlab("Relative read abundance (%)") 
# +
#   scale_color_manual(values = c("#000000","#848484"))

## save plot 

# ggsave(file = "~/prjcts/fish_eDNA/sfjq/data/trees/SFmc-23_SPs_fused-unique-RRA_barplot.pdf", 
#      plot = SFmc_RRA_plot, device = "pdf", width = 24, height = 20, units = "cm", dpi = 600)
ggsave(file = "~/outros/sfjq_temp/trees/SFmc-23_SPs_fused-unique-RRA_barplot.pdf", 
       plot = SFmc_RRA_plot, device = "pdf", width = 30, height = 30, units = "cm", dpi = 600)

















# JQmc ----

## build tree ----
### read 12S db seqs for the species present in pools and select only pool species
names(SFJQ_Sps_seqs) %>% sort() %>% paste0(collapse = '",\n"') %>% cat()

JQ_Sps_seqs <- SFJQ_Sps_seqs[c("JQ_1762_Megaleporinus_elongatus",
                               "JQ_1936_Delturus_brevis",
                               "JQ_2984_Hypomasticus_steindachneri",
                               "JQ_2996_Steindachneridion_amblyurum",
                               "JQ_4290_Megaleporinus_garmani",
                               "JQ_5612_Prochilodus_argenteus_hartii",
                               "JQ_5645_Rhamdia_aff_quelen",
                               "JQ_6584_Australoheros_sp",
                               "JQ_6586_Steindachnerina_elegans-Cyphocharax_gilbert",
                               "JQ_7817_Wertheimeria_maculata",
                               "JQ_7822_Hoplias_malabaricus",
                               "JQ_7880_Astyanax_lacustris",
                               "JQ_7893_Hypostomus_nigrolineatus",
                               "SF_1243_Trachelyopterus_galeatus_JQ_5675",
                               "SF_1248_Gymnotus_carapo_JQ_1631",
                               "SF_1339_Prochilodus_costatus_JQ_2860",
                               "SF_1377_Hoplias_intermedius_brasiliensis")]
### align seqs
JQ_Sps_algn <- DECIPHER::AlignSeqs(myXStringSet = JQ_Sps_seqs, 
                                      refinements = 100,
                                      iterations = 100,
                                      verbose = TRUE)

### generate distance matrix
JQ_Sps_dist <- DECIPHER::DistanceMatrix(myXStringSet = JQ_Sps_algn,
                                            includeTerminalGaps = FALSE,
                                            correction = "Jukes-Cantor",
                                            processors = 20,
                                            verbose = TRUE)

### generate dendrogram/tree from alignment and distance matrix
JQmc_tree <- ape::nj(JQ_Sps_dist)
        # tree <- phangorn::NJ(SFJQ_Sps_dist)
class(JQmc_tree)

### save tree as newick
# ape::write.tree(phy = JQmc_tree,file = "~/prjcts/fish_eDNA/sfjq/data/trees/JQmc-23_SPs_fused-unique_APEtree.nwk")
ape::write.tree(phy = JQmc_tree,file = "~/outros/sfjq_temp/trees/JQmc-23_SPs_fused-unique_APEtree.nwk")

### read tree from file (or stay with the same object)
JQmc_tree <- read.tree("~/prjcts/fish_eDNA/sfjq/data/trees/JQmc-23_SPs_fused-unique_APEtree.nw")
# JQmc_tree <- read.tree("~/outros/sfjq_temp/trees/JQmc-23_SPs_fused-unique_APEtree.nwk")

## species metadata ----
### read table with pools species, input DNA and RRA

### tidy table
tab_curated_JQ <- tab_curated_SFJQ_all_pools %>% 
  filter(Pool %in% c("JQmc")) %>% 
  filter(Status %in% c("Expected")) %>%
  # filter(MC == "SFJQmc") %>% 
  # filter(Pool == "SFJQmc") %>% 
  mutate(Primer = factor(Primer, levels = c("NeoFish","MiFish","Teleo")),
         Species = factor(Species),
         # `revised final ID` = factor(`revised final ID`),
         Normalization = factor(Normalization),
         Pool = factor(Pool),
         Status = factor(Status),
         # `RRA (%)` = as.numeric(`Relative abundance on sample`)*100,
         `RRA (%)` = as.numeric(`RRA (%)`),
         `input DNA (%)` = as.numeric(`input DNA (%)`))


## correct tips labels ----
### check tip names in tree
JQmc_tree$tip.label


### rename tree tips to mach table
JQmc_tree$tip.label
tab_curated_JQ$Species %>% unique() %>% sort()

{
JQmc_tree$tip.label[JQmc_tree$tip.label == "SF_1339_Prochilodus_costatus_JQ_2860"] <- "Prochilodus costatus"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_5612_Prochilodus_argenteus_hartii"] <- "Prochilodus argenteus/hartii"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_6586_Steindachnerina_elegans-Cyphocharax_gilbert"] <- "Cyphocharax gilbert"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_2984_Hypomasticus_steindachneri"] <- "Hypomasticus steindachneri"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_4290_Megaleporinus_garmani"] <- "Megaleporinus garmani"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_1762_Megaleporinus_elongatus"] <- "Megaleporinus elongatus"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_7822_Hoplias_malabaricus"] <- "Hoplias malabaricus"
JQmc_tree$tip.label[JQmc_tree$tip.label == "SF_1377_Hoplias_intermedius_brasiliensis"] <- "Hoplias brasiliensis"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_7880_Astyanax_lacustris"] <- "Astyanax lacustris"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_7893_Hypostomus_nigrolineatus"] <- "Hypostomus nigrolineatus"
JQmc_tree$tip.label[JQmc_tree$tip.label == "SF_1248_Gymnotus_carapo_JQ_1631"] <- "Gymnotus carapo"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_1936_Delturus_brevis"] <- "Delturus brevis"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_7817_Wertheimeria_maculata"] <- "Wertheimeria maculata"
JQmc_tree$tip.label[JQmc_tree$tip.label == "SF_1243_Trachelyopterus_galeatus_JQ_5675"] <- "Trachelyopterus galeatus"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_2996_Steindachneridion_amblyurum"] <- "Steindachneridion amblyurum"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_5645_Rhamdia_aff_quelen"] <- "Rhamdia quelen"
JQmc_tree$tip.label[JQmc_tree$tip.label == "JQ_6584_Australoheros_sp"] <- "Australoheros sp"
}

### check if all names have correspondences
JQmc_tree$tip.label %in% tab_curated_JQ$Species
tab_curated_JQ$Species %in% JQmc_tree$tip.label
tab_curated_JQ$Species[!tab_curated_JQ$Species %in% JQmc_tree$tip.label]




#invert branches to match Siluriformes position on the top


JQmc_tree$edge
JQmc_tree %>% as_tibble() %>% View()


JQmc_tree  %>% 
  ggtree() + 
  # geom_text(aes(label=node)) + 
  theme_tree2() +
  geom_tiplab(offset = 0,align = T)+ 
  xlim(0, 0.42)%>% 
  ggtree::flip(node1 = 21,node2 = 20)
 

rotateNodes(tree = JQmc_tree, "all") %>%
  ggtree() + 
  theme_tree2() +
  geom_tiplab(offset = 0,align = T)+ 
  xlim(0, 0.42)














### convert ape tree to prettier ggtree object
JQmc_tree_plot <- ggtree(JQmc_tree) + 
  theme_tree2() +
  geom_tiplab(offset = 0,align = T)+ 
  xlim(0, 0.42)


### save tree plot
# ggsave(file = "~/prjcts/fish_eDNA/sfjq/data/trees/JQmc-17_SPs_fused-unique-tree_plot.pdf", 
ggsave(file = "~/outros/sfjq_temp/trees/JQmc-17_SPs_fused-unique-tree_plot.pdf", 
       plot = JQmc_tree_plot, 
       device = "pdf", 
       width = 24, height = 24, 
       units = "cm", dpi = 600)



### extract tips order from tree to reproduce on plots
SPs_order_in_JQ_tree <- get_taxa_name(JQmc_tree_plot) %>% rev()

            # SPs_order_in_JQ_tree <- extract_tree_data(tree_plot) %>% 
            #     dplyr::filter(isTip) %>% 
            #     dplyr::pull(label)




## plots ----
### RRA histogram ----

JQmc_RRA_plot <-
  tab_curated_JQ %>%
    unite(Normalization, Primer, col = "Primer_norm",remove = F,sep = " ") %>%
  mutate(Primer_norm = factor(Primer_norm, levels = c("Normalized NeoFish","Normalized MiFish","Normalized Teleo","Non-normalized NeoFish","Non-normalized MiFish","Non-normalized Teleo"))) %>%
  mutate(Species = factor(Species,levels = SPs_order_in_JQ_tree)) %>% 
  ggplot(aes(y = Species,
             x = `RRA (%)`,
             fill = Primer_norm, 
             col = Normalization
             ))+
  geom_bar(stat = "identity", size = 0.3,width = 1,
           # alpha = 0.7,
           position = position_dodge(preserve = "single" ,width = 1.5)) +
           # position = position_dodgejust(preserve = "single" ,width = 1.2)) +
  scale_fill_manual(values = colors6[c(1,3,5,2,4,6)], name = "") +
  geom_point(aes(y = Species,
                 x = `input DNA (%)`,
             colour = Normalization),
             shape = "|",
             size = 3) +
    scale_color_manual(values = c("#747474","#000000")) +
  xlab("Relative read abundance (%)") 


## save plot 

# ggsave(file = "~/prjcts/fish_eDNA/sfjq/data/trees/JQmc-23_SPs_fused-unique-RRA_barplot.pdf", plot = JQmc_RRA_plot, device = "pdf", width = 24, height = 20, units = "cm", dpi = 600)
ggsave(file = "~/outros/sfjq_temp/trees/JQmc-17_SPs_fused-unique-RRA_barplot.pdf", plot = JQmc_RRA_plot, device = "pdf", width = 30, height = 30, units = "cm", dpi = 600)














#arvore ----


# SFJQmc_tree <- read.tree("~/prjcts/fish_eDNA/sfjq/data/12S_full_SFJQmc_fused_SPs_e_3_contams.nwk")




ggplot(SFJQmc_tree) + geom_tree() + theme_tree()

# This is convenient shorthand

# # tree_plot <- 
#   ggtree(SFJQmc_tree) + 
#   theme_tree2() +
#   geom_tiplab(offset = 0,align = T)+ 
#   xlim(0, 0.42) 

  
  tree_plot
  
  
  
    
#TODO create function to build double graph of tree and bars

  
  
# newick tree to plot along the graph
  tree4plot <- SFJQmc_tree

  
# table for ggplot to go alongside the tree (long format, can have factors)
  tbl4plot <- tab_curated_SFJQ
  
# plot to put by side (y axis must be the species in tree (column name == Species))
  plot4tree <- SFJQmc_RRA_plot
  
  
  
### generate plot from plylo object (.nwk read by ape::read.tree)
  

class(SFJQmc_tree)

tree4plot$edge.length %>% sort() %>% sum() 
tree4plot %>% str()


tree_plot$data

# tree_plot <-
  ggtree(tr = tree4plot,layout = "rectangular") +
  theme_tree2() +
  # geom_tiplab(offset = 0,align = T)+
  geom_tiplab(align = T)+
  xlim(0, 0.42)
  
  
### extract tips order from tree to reproduce on plots
SPs_order_in_tree <- ggtree::get_taxa_name(SFJQmc_tree_plot) %>% rev()
  
  
  
  plot4tree$data %>% 
    dplyr::mutate(Species = factor(Species, levels = SPs_order_in_tree))
    
  
  
    
    
    ## check tip names in tree
JQmc_tree$tip.label







### convert ape tree to prettier ggtree object
JQmc_tree_plot <- ggtree(JQmc_tree) + 
  theme_tree2() +
  geom_tiplab(offset = 0,align = T)+ 
  xlim(0, 0.42)
    
    
treeNbar_plot  <- function(){
    
  }  
  
    
  
  
  
  
  
  
  

# https://www.r-bloggers.com/2016/12/add-layer-to-specific-panel-of-facet_plot-output-2/
# facet_plot(tree_plot, panel = 'Stacked Barplot', 
#            data = tab_curated_SFJQ, geom = geom_histogram,
#            mapping = aes(x = Species,y =`RRA (%)`,  fill = as.factor(Primer)),
#            stat='identity',position = 'dodge' )
# 
# # 
# p3 <- facet_plot(tree_plot, panel='bar', data=tab_curated_SFJQ, geom=geom_bar, 
#                  aes(x=`RRA (%)`, y=Species,fill = Primer),
#                  stat = "identity",
#                  position = "dodge")


#https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mbe/35/12/10.1093_molbev_msy194/3/msy194_supp.pdf?Expires=1648204036&Signature=gZC6A1vfyaUWOyL~fKn8wxgY~3fbTBI3jPOGbVtwZSzv3jlXISjCahA37gwR3QTr6oN0SK-bdwAlHQyaPpkdj2~yi5scNQXAQrUi0EQNOqkOo3HUvFvCr-Dir2y7N03vIo5urr1n2idrPclTXtTRtiu7avn255T5eg~cXv0NBNUgiVFcwHHnZ81qQUrSdiA54wIvEs~RF18DYkp-Gla1CJT0eUGuYF8LfFXG5Dq1CgcZV0qGs0fKgfIKRlAT~AP25Xxkdh20RzAkqgBFvxp0JazrVOz5uvdok3uSu3023etErTxhaW7rm67VkCUBVxRgtG8GdFT3fOFJAsPg26Wagw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA

# 
# tree_plot %<+% tab_curated_SFJQ + 
#   geom_histogram(aes(y = Species,
#                      x = `RRA (%)`,
#                      fill = Primer),
#                  stat = "identity",
#                  position = "dodge")

colnames(tab_curated_SFJQ)[colnames(tab_curated_SFJQ) == "Species"] <- "tip.label"



p2 <- 
  facet_plot(p = tree_plot, panel = "SNP", data = tab_curated_SFJQ, geom = geom_histogram,
                 mapping=aes(y = tip.label, x = `RRA (%)`, fill = Primer),
                 stat = "identity",
                 position = "dodge") +
# %>%
#   facet_plot("Trait", bar_data, ggstance::geom_barh,
#              aes(x = dummy_bar_value, color = location, fill = location),
#              stat = "identity", width = .6) +
  theme_tree2(legend.position=c(.05, .85))
print(p2)



p2 <- tree_plot + geom_facet(panel = "RRA (%)",
                       data = tab_curated_SFJQ,
                       geom = geom_bar,
                       mapping = aes(x=`RRA (%)`, fill = Primer),orientation = 'y',stat="identity")





extract_tree_data <- function(tree_disp, displayorder=TRUE) {
    td_out <- tree_disp$data
    if (displayorder) {
       td_out <- dplyr::arrange(td_out,y)
    }
    return(td_out)
}

SPs_order_in_tree <- extract_tree_data(tree_plot) %>% 
    dplyr::filter(isTip) %>% 
    dplyr::pull(label)









## add resuts to tree ----

# import tree ----
```

## dendrograms of ASVs & db species

## fold change standard deviation

## NeoFish redesign

<br>

**This is a partial report, intended to show the current state of
analyses. Many procedures and conclusions might change as the pipeline
evolves. If you notice errors/mistakes/typos, or have any suggestions,
we would be glad to know. *<heronoh@gmail.com>***
