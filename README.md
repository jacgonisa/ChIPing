# ChIPing: A pipeline for ChIP-seq data analysis 

Authors: 
* González Isa, Jacob (jacobgisa17@gmail.com) [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/jacobgisa.svg?style=social&label=Follow%20%40JacobGIsa)](https://twitter.com/jacobgisa)

* López Asensio, Pilar
* Pérez López, José Ignacio

### Course: Bioinformatics and Genomic Analysis
### BSc in Biochemistry
### University of Seville
 

## 1. Introduction 

This repository consists of three bash scripts (**chipipe.sh**, **sample_proc.sh** and **peak_call.sh**) and one R script (**chip_bash.R**). It runs in Unix environment.

ChIPing aims at analysing any given ChIP-seq samples of Transcription Factors (TFs) in *Arabidopsis thaliana* model organism. Each **sample** has a ChIP and control/input sequencing data. The number of samples is users' choice! 


## 2. Dependencies

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

> - **installation_directory.** Where you install the package.
> - **working_directory.** Where your analysis are saved.
> - **experiment_name.** The name the folder that contains everything.
> - **number_replicas.** The number of replicates.
> - **path_genome.** The path to access the genome file of the organism; e.g. /home/user/user_genomes/organism_genome.fa
> - **path_annotation.** The path to access the genome's annotation file of the organism; e.g. /home/user/user_annotations/organism_anno.gtf
> - **path_sample_chip_i.** (i= 1,2,3...) The path to access the ChIP-seq data of each sample; e.g. /home/user/user_experiment/sample_chip_i.fq.gz. If you have paired end files, you must write both paths in the same row, separated by space.
> - **path_sample_input_i** (i= 1,2,3...) The path to access the input data of each sample; e.g. /home/user/user_experiment/sample_input_i.fq.gz. If you have paired end files, you must write both paths in the same row, separated by space.
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
 * **`experiment_name/samples/replica_i`:**
 
    - **`/chip`**: 
        > * `/chip_i.bam`: compressed binary version of a SAM file.
        > * `/sample_chip_i_fastqc.html`: the quality metrics of the reads. Information about the quality score distribution across the reads and useful information such as adapter contamination or other overrepresented sequences, among others.
 
    - **`/input`:**
        > * `/input_i.bam`.
        > * `/sample_input_i_fastqc.html`.
 
     - **`/replica_results`:** 
        > * `/i_peaks.narrowPeak`: peak file generated with macs2 for each replica. The workflow can be modified to study histones PTMs. In that case, the output of macs should be a i_peaks.broadPeak.
        > * `/i_summits.bed`: peak summits locations for every peak. File needed to find the motifs at the binding sites. It will only be generated if running in .narrowPeak mode.
        > * `/i_model.r`: R script to produce a PDF image of the model based on the data. Once the script is run, the PDF will automatically appear in  /results directory.
        > * `/i_peaks.xls`: a tabular file with information about called peaks.


 * **`experiment_name/results`:**

   - Motifs finding
    Motif analysis with [`HOMER`](http://homer.ucsd.edu/homer/index.html). By default, it performs de novo motif discovery as well as checking the enrichment of known motifs.
    
       > * `/knownResults.html`: formatted output of known motif finding.
       > * `/knownResults/ `: directory with known<#>.motif and known<#>.logo.svg for each motif.
       > * `/knownResults.txt`: statistics about known motif enrichment. By default, it is opened in Text.editor. For optimized view, open manually in Excel.

       > * `/homerResults.html`: formatted output of de novo motif finding.
       > * `/homerResults/ `: directory with homer<#>.motif and homer<#>.logo.svg for each motif.
       > * `/homerMotifs.motifs<motif_length>`: de novo motif finding, separated by motif length, and representing separate runs of the algorithm.
       > * `/homerMotifs.all.motifs`: all the homerMotifs.motifs<motif_length> files.

       > * `/seq.autonorm.tsv`: autonormalization statistics for lower-order oligo normalization.
       > * `/motifFindingParameters.txt`: command used to execute findMotifsGenome.pl.


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
      > * `/kegg_images/`: directory with the pathways without marked enzymes in png and xml format generated with [`pathview`](https://bioconductor.org/packages/release/bioc/html/pathview.html)). Toolset for pathway-based data integration and visualization. It maps and renders user data on relevant pathway graphs. 


## 7. Case study

As an example, a case of study is PRR5, which acts as a transcripcional repressor in Arabidopsis thaliana, more specifically as a regulating factor of the circadian clock expression. The circadian clock of plants regulates a wide range of processes, including different processes regulated by cell cycle; sucj as hypocotyle elongation before sunrise, or cold-stress response proteins. In this way, PRR5 allows the expression of different genes in different day moments, by the regulation of plants circadian clock.

In the present study, a ChIP-seq analysis to study th main proteins of the PRR family target genes is carried out.

The samples to process are a chip and a control sample (input), with two replics, with reads for the chromosome 1 of Arabidopsis thaliana. The chromosome 1 FASTA file (genome.fa), and the chromosome 1 GTF file (annotation.gtf), reference genome and annotation, respectively, are also given. 

Before sample processing (sample_proc.sh), it is necessary to generate a workspace and to build an index. After that, the firt step of the sampke processing is the sample processing control. The result of this analysis can be found in the html of the quality control of the different samples (fastqc.html). This outputs show: basic statistics, per base sequence quality, per base quality scores, per sequence GC content, etc. For both samples, a green check mark appears in all sections, so the quality analysis results are good. 

*USCAR EL ARCHIVO DE LOS DATOS DEL MAPEO* Once verified the quality control is good, the reads mapping is done. A summary of the results of the mapping can be found on the mapping_value.txt, found on the results samples. In this case, the number of reads and overall alignment of the chip sample was 1561512 and 99.96%, respectively. For the control sample, those values were 235150 and 99.50%. As mentioned before, the alignment can be viewed on IGV. As PRR5 is related to the circadian clock, it would not be unusual to find that one of its targets is LHY -locus: AT1G01060-. When bam, and bam.bai files for chip and control samples are loaded in IGV, on the region of the mentioned locus, unspecific binding is found for control file, while for the chip file, can be found real binding and background noise.

Peak calling (peak_call.sh), the next step in our workflow, is a computational method used to identify areas in the genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing experiment. MACS2 peak calling algorithm was used.

Finally, using the toolbox HOMER, 242 known binding motifs were identified as PRR5 binding motifs, in addition to 4 PRR5 binding motifs found de novo. One of these binding motifs is a CCT motif.

With all this information, it becomes clear that PRR5 regulates different genes that participate in processes related to the circadian clock, like flowering, hypocotyl elongation, and even mediating cold-stress response.
