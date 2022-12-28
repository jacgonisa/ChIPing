# ChIPing: A pipeline for ChIP-seq data analysis 

Authors: 
* González Isa, Jacob
* López Asensio, Pilar
* Pérez López, José Ignacio

### Course: Bioinformatics and Genomic Analysis
### BSc in Biochemistry
### University of Seville
 

## 1. Introduction 

This repository consists of three bash scripts (**chipipe.sh**, **sample_proc.sh** and ) and one R script (**chip_bash.R**). It runs in Unix environment.

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
   * Submit via SGE chip_sample_processing.sh and control_sample_processing for each sample.
 2. **chip_sample_processing.sh** and **control_sample_processing.sh** 
    * Quality control ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).
    * Map to reference genome ([`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)).
    * Generate sorted bam file ([`samtools`](http://www.htslib.org)).
    * Submit via SGE peak_determination.sh for each sample.
 3. **peak_determination.sh**
    * Peak determination ([`masc2 callpeak`](https://github.com/macs3-project/MACS)).
    * Peak annotation by submitting *target_genes.R* for each sample.
    * Homer-motifs-finding ([`findMotifsGenome.pl`](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)).
    * Overlapping target genes for replicates, GO and KEGG enrichment by submitting *exp_analysis.R* for each experiment.
 4. **target_genes.R**
    * Install packages if neccesary ([`BiocManager`](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)).
    * Read arguments.
    * Read peak file ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)).
    * Definition of promoter region ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)).
    * Peak annotation ([`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html),[`DO.db`](http://bioconductor.org/packages/release/data/annotation/html/DO.db.html)).
    * Extract target genes from annotation.
 5. **exp_analysis.R**
    * Install packages if neccesary ([`BiocManager`](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)).
    * Read arguments.
    * Read target gene files for each replicate.
    * Extraction of overlapping genes in all replicates ([`VennDiagram`](https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf)).
    * Gene ontology enrichment ([`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html),[`TxDb.Athaliana.BioMart.plantsmart28`](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html),[`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)).
    * KEGG pathway enrichment ([`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html),[`TxDb.Athaliana.BioMart.plantsmart28`](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html),[`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html),[`pathview`](https://bioconductor.org/packages/release/bioc/html/pathview.html)).
    
    
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
 * **genome:**
   - Reference Genome 
   - Index of Reference Genome
 * **annotation:**
   - Reference Annotation
 * **samples:**
   - **Sample_i:**
  
      - **chip**: 
        > * **BAM file:** compressed binary version of a SAM file.
        > * **FastQC:** This html file reports the quality metrics of the reads. It provides information about the quality score distribution across the reads and useful information such as adapter contamination or other overrepresented sequences, among other issues.
        > * **stats.txt:** This file contains **bowtie2** alignment stats.
 
      - **input:**
        > * **BAM file:** compressed binary version of a SAM file.
        > * **FastQC:** This html file reports the quality metrics of the reads. It provides information about the quality score distribution across the reads and useful information such as adapter contamination or other overrepresented sequences, among other issues.
        > * **stats.txt:** This file contains **bowtie2** alignment stats.
 
      - **replica_results:** 
        > * **i_peaks.narrowPeak:** peak file generated with macs2 for each replica. The workflow can be modified to study histones PTMs. In that case, the output of macs should be a i_peaks.broadPeak.
        > * **i_summits.bed:** file in BED format containing the peak summits locations for every peak. T This file is needed to find the motifs at the binding sites. This file will only be generated if running in narrow peak mode.
        > * **i_model.r**: R script to produce a PDF image of the model based on the data. Once the script is run, the PDF will automatically appear in the results directory.
        > * **i_peaks.xls:** a tabular file which contains information about called peaks with useful information such as chromosome name, start position of peakend among other relevant details.
 * **results**
    ##### Motifs finding
     Motif analysis carried out with [`HOMER`](http://homer.ucsd.edu/homer/index.html). In this pipeline is used [`findMotifsGenome.pl`](http://homer.ucsd.edu/homer/ngs/peakMotifs.html), which manage all the steps for discovering motifs in genomic regions. By default, this will perform *de novo* motif discovery as well as check the enrichment of known motifs.


<details markdown="1">
    <summary>Output files</summary>
    
`knownResults/ directory`: contains files for the knownResults.html webpage, including known<#>.motif files for use in finding specific instance of each motif as well as an image.svg for each motif.   

`homerResults/ directory`: contains files for the homerResults.html webpage, including motif<#>.motif files for use in finding specific instance of each motif as well as an image.svg for each motif.
    
`homerMotifs.motifs<motif_length>`: output files from the de novo motif finding, separated by motif length, and represent separate runs of the algorithm.

`homerMotifs.all.motifs`: file composed of all the homerMotifs.motifs<motif_length> files.

`knownResults.html`: formatted output of known motif finding.

`knownResults.txt`: text file containing statistics about known motif enrichment. By default is opened in Text.editor. For optimized view, open manually in Excel.

`homerResults.html`: formatted output of de novo motif finding.

`seq.autonorm.tsv`: autonormalization statistics for lower-order oligo normalization.

`motifFindingParameters.txt`: command used to execute findMotifsGenome.pl.

</details>

##### ChIP-seq Analysis information:
    
<details markdown="1">
    <summary>R Script Output files</summary>
    
    `Rplots.pdf`: This file contains graphical representations: Histogram of ChIP peaks over Chromosome(s), pieplot of the distribution of TF binding regions (or epigenetic mark) as well as a plotDistToTSS which shows the distribution of genomic loci relative to TSS.
 
</details>

##### Regulome:

The regulome is defined as the global set of genes regulated by a transcription factor under certain circumstances, and it can be inferred from the cistrome, the set of genetic positions where the transcription factor binds under those certain conditions, by the association of ChipSeq peaks to target genes via the NDG criteria

<details markdown="1">
    <summary>Output files</summary>
    
    `Regulome.txt`: File with a list of the genes predicted to be affected by the TF binding or histone modification in study.
 
</details>

##### GO terms Analysis:


 As for the GO terms analysis, chiptube calculates the GO terms enrichment for biological processes (bp), molecular functions (mf) and cellular components (cc), and all the information is saved as tables in tsv format, as well as in plots that are represented in a pdf file for each one of the three categories (goplots, barplots, dotplots and cnetplots). Finally, as for GO terms enrichment, KEEG pathways enrichment information is saved as a table in a tsv file, and the proper pathways are shown as png files in this directory, while the xml and png (without marked enzymes) files are collected in kegg_images directory.
 
## 7. Case study

As an example, the repressor effect of Arabidopsis thaliana PRR5 protein, as a regulating factor of circadian clock expression, is studied. Circadian clock of plants regulates a wide range of processes, such as hypocotyl elongation before sunrise, or cold-stress response proteins. This way, the circadian clock allows the expression of different genes in different day moments, regulating their expression temporarily. 

Recently, some transcription factors that directly regulate clock-output genes have been discovered. One of these transcription factors is the family of genes called pseudo-response regulator (PRR). However, it remains still unknown which genes are regulated by this family of genes, and how are they regulated. Hence, in the present study, a ChIP-seq analysis to determine one of the main proteins of de PRR family target genes is carried out.

The samples to process are a chip and a control sample, with reads for the chromosome 1 of Arabidopsis thaliana. The chromosome 1 FASTA file, and the chromosome 1 GTF file, reference genome and annotation, respectively, are also given. 

As previously explained, the first step on the sample processing, after generating workspace and building index, is the sample quality control. On the HTML of the quality control, different sections can be found: basic statistics, per base sequence quality, per base quality scores, per sequence GC content, etc. For both samples, the quality control analysis results are good. This is really easy to know: generally, in the HTML, if a green check mark appears in the different sections, that means the quality analysis went okay.

After the quality analysis, the reads mapping is done. A summary of the results of the mapping can be found on the mapping_value.txt, found on the results samples. In this case, the number of reads and overall alignment of the chip sample was 1561512 and 99.96%, respectively. For the control sample, those values were 235150 and 99.50%. As mentioned before, the alignment can be viewed on IGV. As PRR5 is related to the circadian clock, it would not be unusual to find that one of its targets is LHY -locus: AT1G01060-. When bam, and bam.bai files for chip and control samples are loaded in IGV, on the region of the mentioned locus, unspecific binding is found for control file, while for the chip file, can be found real binding and background noise.

Even though the mapping can be viewed, with this visualisation it cannot be distinguished which reads belong to real PRR5 binding from those that are background noise. To remove the background, and retain significant reads, the peaks calling is carried out. The results of the peaks calling can also be viewed on IGV. Again, a peak can be found on the LHY locus. This way, it can be concluded that PRR5 probably regulates LHY gene.

These positions where the transcription factor binds significantly, on a condition given, are known as the cistrome. However, to determine which genes are regulated by PRR5, its regulome, the files from the alignment are run on the R script mentioned before. Once the regulome is determined, this can be used to make a GSEA, and identify relevant GO terms or KEGG Pathways. Some of the enrichments are related to the karrikin metabolism, or the photosynthesis.

Finally, using the toolbox HOMER, 242 known binding motifs were identified as PRR5 binding motifs, in addition to 4 PRR5 binding motifs found de novo. One of these binding motifs is a CCT motif.

With all this information, it becomes clear that PRR5 regulates different genes that participate in processes related to the circadian clock, like flowering, hypocotyl elongation, and even mediating cold-stress response.
