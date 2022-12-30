# ChIPing: A pipeline for ChIP-seq data analysis 

Authors: 
* González Isa, Jacob (jacobgisa17@gmail.com) [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/jacobgisa.svg?style=social&label=Follow%20%40JacobGIsa)](https://twitter.com/jacobgisa)

* López Asensio, Pilar (pilarlopezasensio@gmail.com)
* Pérez López, José Ignacio (naperlopez@gmail.com) [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/nach_per.svg?style=social&label=Follow%20%40nach_per)](https://twitter.com/nach_per)

#### Course: Bioinformatics and Genomic Analysis
### BSc in Biochemistry
#### University of Seville
![](https://www.us.es/sites/default/files/logoPNG_3.png) 

 

## 1. Introduction 

This repository consists of three bash scripts (**chipipe.sh**, **sample_proc.sh** and **peak_call.sh**) and one R script (**chip_bash.R**). It runs in Unix environment.

ChIPing aims at analysing any given ChIP-seq samples of Transcription Factors (TFs) in *Arabidopsis thaliana* model organism. Each **sample** has a ChIP and control/input sequencing data. The number of samples is users' choice! 


## 2. Dependencies and installation

This pipeline requires the installation of:

  * [`Slurm cluster management`](https://slurm.schedmd.com/)
  * [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [`SAMTOOLS`](https://sourceforge.net/projects/samtools/files/samtools/)
  * [`MACS2`](https://github.com/macs3-project/MACS)
  * [`R`](https://www.r-project.org), with the following packages:
    - [`BiocManager`](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
    - [`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
    - [`TxDb.Athaliana.BioMart.plantsmart28`](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html)
    - [`DO.db`](http://bioconductor.org/packages/release/data/annotation/html/DO.db.html)
    - [`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)
    - [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
    - [`pathview`](https://bioconductor.org/packages/release/bioc/html/pathview.html)
  * [`Homer`](http://homer.ucsd.edu/homer/download.html)

Once dependencies have been installed, you might download this repository from GitHub, pasting in your command line:
```bash
git clone https://github.com/jacgonisa/ChIPing.git
```
## 3. Summary

1. **chipipe.sh**
   * Load and read parameters.
   * Generation of workspace.
   * Creation of a genome index ([`bowtie2-build`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)).
   * Processing individual samples with next script:
 2. **sample_proc.sh**
    * Quality control ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).
    * Mapping to the reference genome ([`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)).
    * Generation of sorted .bam file ([`samtools`](http://www.htslib.org)).
    * Peak calling ([`masc2 callpeak`](https://github.com/macs3-project/MACS)).
     When finished for each replicate, next script:
 3. **peak_call.sh**
    * Intersection of replicates' results.
    * Finding motifs ([`findMotifsGenome.pl`](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)).
    * Visualization and statistics with next script:
 4. **chip_bash.R**
    * Load parameters.
    * Reading of peak file ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)).
    * Definition of promoter region ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)).
    * Peak distribution along the genome.
    * Peak annotation ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html),[`DO.db`](http://bioconductor.org/packages/release/data/annotation/html/DO.db.html)).
    * Selection of peaks within promoters and **regulome** identification.
    * GO terms enrichment ([`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html),[`TxDb.Athaliana.BioMart.plantsmart28`](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html),[`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)).
    * KEGG terms enrichment ([`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html),[`TxDb.Athaliana.BioMart.plantsmart28`](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html),[`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html),[`pathview`](https://bioconductor.org/packages/release/bioc/html/pathview.html)).
    
 
## 4. Usage

```sh
bash chipipe.sh <params_file.txt>
```

The main script is **chipipe.sh**. The input is a file containing some parameters to be specified. An example of the model file **params_file.txt** is found ChIPing/my_test/params_file.txt
* General parameters:
> - **experiment_name.** The name the folder that contains everything.
> - **number_replicas.** The number of replicates.
* Directories:
> - **installation_directory.** Where you install the package.
> - **working_directory.** Where your analysis are saved.

* Paths:
> - **path_genome.** The path to access the genome file of the organism; e.g. /home/user/user_genomes/organism_genome.fa
> - **path_annotation.** The path to access the genome's annotation file of the organism; e.g. /home/user/user_annotations/organism_anno.gtf
> - **path_sample_chip_**<**i**> (i= 1,2,3...) The path to access the ChIP-seq data of each sample; e.g. /home/user/user_experiment/sample_chip_i.fq.gz. If you have paired end files, you must write both paths in the same row, separated by space.
> - **path_sample_input_**<**i**> (i= 1,2,3...) The path to access the input data of each sample; e.g. /home/user/user_experiment/sample_input_i.fq.gz. If you have paired end files, you must write both paths in the same row, separated by space.

* Customisable parameters:
> - **universe_chromosomes.**. The ID(s) of the chromosome(s) of your organism you want to use as your genetic universe for GO and KEGG terms enrichment, separated by commas without spaces; e.g. 2,3. In order to use all the available chromosomes, write "all".
> - **p_value_cutoff_go.** The p-value threshold for GO terms enrichment statistical analysis.
> - **p_value_cutoff_kegg.** The p-value threshold for KEGG pathways enrichment statistical analysis.
> - **type_of_peak.** The shape of the peaks. Two posible values: 1 (narrow peaks; for TF binding) or 2 (broad peaks; for histone modifications).
> - **single_or_paired**. The type of reads of the files. The value of this parameter must be either 1 (single end reads) or 2 (paired end reads). For this pipeline, only single-end reads can be run.
> - **tss_upstream.** The upstream number of base pairs that defines the TSS region. Only positive values.
> - **tss_downstream.** The downstream number of bases pairs that defines the TSS region. Only positive values.


## 5.Output
### Working Directory 
 * **`experiment_name/genome`:**
   - Reference Genome (`/genome.fa`) 
   - Index of Reference Genome
 * **`experiment_name/annotation`:**
   - Reference Annotation (`/annotation.gtf`) 
 * **`experiment_name/samples/replica_<i>`:**
 
    - **`/chip`**: 
        > * `/chip_<i>.bam`: compressed binary version of a SAM file (the output of mapping).
        > * `/sample_chip_<i>_fastqc.html`: the quality metrics of the reads. Information about the quality score distribution across the reads and useful information such as adapter contamination or other overrepresented sequences, among others.
 
    - **`/input`:**
        > * `/input_<i>.bam`.
        > * `/sample_input_<i>_fastqc.html`.
 
     - **`/replica_results`:** 
        > * `/<i>_peaks.narrowPeak`: peak file generated with macs2 for each replica. The workflow can be modified to study histones PTMs. In that case, the output of macs should be `<i>_peaks.broadPeak`.
        > * `/<i>_summits.bed`: peak summits locations for every peak. File needed to find the motifs at the binding sites. It will only be generated if running in .narrowPeak mode.
        > * `/<i>_model.r`: R script to produce a PDF image of the model based on the data. Once the script is run, the PDF will automatically appear in  /results directory.
        > * `/<i>_peaks.xls`: a tabular file with information about called peaks.


 * **`experiment_name/results`:**

   - Motifs finding
    Motif analysis with [`HOMER`](http://homer.ucsd.edu/homer/index.html). By default, it performs de novo motif discovery as well as checking the enrichment of known motifs.
    
       > * `/knownResults.html`: formatted output of known motif finding.
       > * `/knownResults/ `: directory with `known<#>.motif` and `known<#>.logo.svg` for each motif.
       > * `/knownResults.txt`: statistics about known motif enrichment. By default, it is opened in Text.editor. For optimized view, open manually in Excel.

       > * `/homerResults.html`: formatted output of de novo motif finding.
       > * `/homerResults/ `: directory with `homer<#>.motif` and `homer<#>.logo.svg` for each motif.
       > * `/homerMotifs.motifs<motif_length>`: de novo motif finding, separated by motif length, and representing separate runs of the algorithm.
       > * `/homerMotifs.all.motifs`: all the `homerMotifs.motifs<motif_length>` files.

       > * `/seq.autonorm.tsv`: autonormalization statistics for lower-order oligo normalization.
       > * `/motifFindingParameters.txt`: command used to execute `findMotifsGenome.pl`.


    - ChIP-seq Analysis information:
       > * `/Rplots.pdf`: every graphical representation. This is, histogram of ChIP peaks over Chromosome(s), pieplot of the distribution of TF binding regions (or epigenetic mark) and a plotDistToTSS, which shows the distribution of genomic loci relative to TSS.
 

   - Regulome:
The regulome is defined as the global set of genes regulated by a transcription factor under certain circumstances, and it can be inferred from the cistrome, the set of genetic positions where the transcription factor binds under those certain conditions, by the association of ChipSeq peaks to target genes via the NDG criteria
      > * `/regulome.txt`: list of the genes predicted to be regulated by the TF (or histone modification).


   - GO terms and Kyoto Encycopedia of Genes and Genomes (KEGG) Analysis:
GO terms enrichment are calculated for biological process (BP). For molecular function or cellular components analysis, customize the R script (*chip.R*) by changing the enrichGO function argument to the desired parameter. 
KEGG analysis is at the level of complete metabolic or regulatory pathways, rather than at the level of specific biological functions or processes.
     > * `/kegg_terms.tsv`: Table separated by tab with the results of the KEGG analysis  
     > * `/go_terms.tsv`:Table separated by tab with the results of the GO terms 
     > * `/plots_go_bp.pdf`: Plots repressenting GO terms.
     > * `/kegg_images/`: directory with the biochemical pathways in .png and .xml format generated with [`pathview`](https://bioconductor.org/packages/release/bioc/html/pathview.html). Toolset for pathway-based data integration and visualization. It maps and renders user data on relevant pathway graphs. 


## 6. Case study: ChIPing put into practice!

We studied a transcription factor known as PRR5 with this repository. PRR5 acts as a transcripcional repressor in *Arabidopsis thaliana*, more specifically as a regulating factor of the circadian clock expression. The circadian clock of plants regulates a wide range of processes, such as hypocotyle elongation before sunrise or cold-stress response. Therefore, PRR5 allows the expression of different genes in throughout the day, by regulating plants' circadian clock.

In our previous study, we performed ChIP-seq analysis to get further insight of how PRR5 acts.

The experimental design consists of two biological replicates (samples), each one with a ChIP and a control (input). Reads belong to the chromosome 1 of *Arabidopsis thaliana*. The chromosome 1 FASTA file (genome.fa) along its annotation file (annotation.gtf) as well as the samples are found within `my_test` folder. 

Before sample processing, it is necessary to generate a workspace and build an index (starting with **chipipe.sh**). After that, the firt step of the sample processing is the processing control. The result of this analysis can be found in the .html of the quality control of the different samples (`fastqc.html`). This output show basic statistics, per base sequence quality, per base quality scores, per sequence GC content, and so on. For both samples, a green check mark appears in all sections, so the quality analysis was successful. 

Once verified the quality of the samples, reads are mapped to the reference genome is done (**sample_proc.sh**). Mapping is successful.

Peak calling, the next step in our workflow, is a computational method used to identify chromatin regions that have been enriched with aligned reads as a consequence of PRR5 binding ([`Learn about ChIP-seq!`](https://www.illumina.com/techniques/sequencing/dna-sequencing/chip-seq.html#:~:text=ChIP%2DSeq%20identifies%20the%20binding,exonuclease%20to%20trim%20unbound%20oligonucleotides)).

The output of the peak calling is `merged_2.narrowPeak`. After merging (**peak_call.sh**), the total number of peaks is 902. The number of peaks of each replic is the same as the total, in this example, replicates 1 and 2 are identical. This means that PRR5 binds to chromatin in 902 loci in chromosome 1.

TFs bind to specific DNA sequences (motifs). In few words, there are TF-associated sequences in the genome. Their analysis is useful to predict the binding of TFs and, consequently, regulatory effects. Aiming at this, we identified 223 known PRR5 binding motifs (`knownResults.html`). Additionally, 7 PRR5 binding motifs were found de novo (`homerResults.html`).

The final step was submitting the peaks' positions to an R script (**chip_bash**) that gave valuable info about where peaks are located along the genomic regions (promoters, intergenic, etc.) and a list of genes which might be putatively regulated by PRR5 (`regulome.txt`). However, we will not get the full potential of our ChIP-seq data unless we merge our Regulome and Differentially Expressed Genes data (RNA-seq), what would accurately tell which genes are truly regulated by PRR5 and how (up/down- regulation).

The GO terms file (`go_trms.tsv`) confirmed that PRR5 is involved in different processes related to plants circadian's clock. Among them, the response to light stimulus and cold aclimatation are highlighted. The KEGG terms file (`kegg_terms.tsv`) links PRR5 and plant hormone signal transduction, MAPK signaling pathway and plant photosyntesis.

Taken all this together, it falls into place that PRR5 is a main regulator of the circadian clock.

![](https://www.activemotif.com/uploads/images/web_site/services-hp-dec2015/chip-seq-services/chip-service-main-banner-web_2.png)
