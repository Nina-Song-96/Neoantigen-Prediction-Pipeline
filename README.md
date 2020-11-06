# Neoantigen Prediction Pipeline

  * Author: Nina Song

## Overview of the Project

  * Acccumulating abberrations in tumor cells lead to present tumor-specific mutabt peptides (neoantigens) on MHC molecules on their surface, which can distinguish tumor cells from normal cells to the immune system and are thus targets for cancer immunotherapy. Technology Platfroms and algorithms enable indentification of tumor neoantigens, such as the use of next-generation Sequencing to screen the exome and transcriptome of a tumor. The process of neoantigen identification entails indentification of tumor-specific protein sequnces and prioritization of candidate neoantigens according to their abundance, presence on patient HLA molecules and ability to elicit a T-cell response. [<sup>1</sup>](#refer-anchor-1)
  * In this project, a computational pipeline, pVACseq, is used to predict the neoantigens present in melanoma with a focus on high-throughput sequencing data. Personalized Variant Antigens by Cancer Sequencing (pVAC-Seq) is a flexible, streamlined computational workflow for identification of neoantigens that integrates tumor mutation and expression data (DNA- and RNA-seq).[<sup>2</sup>](#refer-anchor-2) pVAC-seq is available at https://github.com/griffithlab/pVAC-Seq.

## Data Demonstration

 * Metastatic melanoma cell line (COLO829) is required in this project.
 * FASTQs for tumor and germline exomes of COLO829 is available in trgn.usc.edu.

## Proposed Analysis and Tools &Software

### Prepared Input Data

#### Alignment and variant calling

| Tool                                 | Function                                      |
| ------------------------------------ | --------------------------------------------- |
| BWA (version 0.5.9)                  | Used for alignment with default parameters    |
| Picard MarkDuplicates (version 1.46) | Deduplicates the resulting alignments         |
| HLAminer (version 1)                 | Generate HLA class I haplotype of the patient |

#### Somatic Mutation Detection and List Variants

| Tool                           | Function                                     |
| ------------------------------ | -------------------------------------------- |
| Strelka (version 1.0.10)       | Mutation Detection                           |
| Variant Effect Predictor (VEP) | Properly format a list of annotated variants |

### Perform Epitope Prediction

#### Generate FASTA file of peptide sequences

 * Generate a FASTA file of 17-21 mers

#### Run epitope prediction software

| Tool        | Function                                                     |
| ----------- | ------------------------------------------------------------ |
| NetMHC v3.4 | Predict high affinity peptides that bind to the HLA class I molecule |

#### Parse and filter the ouput

 * Filter: Binding-based

#### Filter neoepitope candidates

* Sequencing-based filter: Depth-based Filters and Gene Expression

### Results: Prioritized list of neoantigen vaccine candidates

# Milestone Update

## Input Data Preparation 

### Check the data and unzip the fastq files

* FASTQs for tumor and germline exomes of COLO829 on [TRGN.usc.edu](http://TRGN.usc.edu)

* _C1_ stands for germline samples (normal), and _T1_ stands for tumor samples

```{bash}
gunzip FDAV_COLO829m0x100_3_CL_Whole_C1_KHS5X_K12375_C59M4ACXX_GATGAATC_L001_R1_001.fastq.gz
gunzip FDAV_COLO829m100x0_3_CL_Whole_T1_KHS5X_K12373_C59M4ACXX_GACAGTGC_L001_R1_001.fastq.gz
```



### Alignment using BWA

* BWA is a software package for mapping DNA sequences against a large reference genome, such as the human genome.
* Download reference genome database from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz. In this project, hg38 is used as the reference. BWA first needs to construct the FM-index for the reference genome (the **index** command). **BWA-MEM** supports long reads and chimeric alignment, and is generally recommended as it is faster and more accurate
* Sample usage: `./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz`

```{bash}
git clone https://github.com/lh3/bwa.git
cd bwa; make
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
./bwa index hg38.fa
./bwa mem -t 4 hg38.fa ~/510Final/data/colo829/FDAV_COLO829m0x100_3_CL_Whole_C1_KHS5X_K12375_C59M4ACXX_GATGAATC_L001_R1_001.fastq ~/510Final/data/colo829/FDAV_COLO829m0x100_3_CL_Whole_C1_KHS5X_K12375_C59M4ACXX_GATGAATC_L001_R2_001.fastq | gzip -3 > normal-aln-pe-L001.sam.gz
./bwa mem -t 4 hg38.fa ~/510Final/data/colo829/FDAV_COLO829m100x0_3_CL_Whole_T1_KHS5X_K12373_C59M4ACXX_GACAGTGC_L001_R1_001.fastq ~/510Final/data/colo829/FDAV_COLO829m100x0_3_CL_Whole_T1_KHS5X_K12373_C59M4ACXX_GACAGTGC_L001_R2_001.fastq | gzip -3 > tumor-aln-pe-L001.sam.gz
```

* Here, two `normal-aln-pe-L00x.sam.gz` and two `tumor-aln-pe-L00x.sam.gz` are the result of alignment to the reference genome. _L00X_ stands for different lanes from sequencing machine. They will be merged later on.
* All of SAM files should be converted into sorted BAM files for further analysis. These steps are completed with `samtools`

```{bash}
gunzip normal-aln-pe-L001.sam.gz
samtools view -b -S normal-aln-pe-L001.sam > normal-aln-pe-L001.bam
samtools sort normal-aln-pe-L001.bam > normal-aln-pe-L001.sorted.bam
samtools index normal-aln-pe-L001.sorted.bam
```



### Merge read groups and mark duplicates per sample (aggregation + dedup)

* Picard MarkDuplicates is a useful tool to deduplicates the resulting alignments. After clone the repo and run Picard jar with all dependencies included, those resulting jar will be in `bulid/libs`.

* Sample usage of picard MarkDuplicates: 

  ```
  java -jar picard.jar MarkDuplicates \
        I=input.bam \
        O=marked_duplicates.bam \
        M=marked_dup_metrics.txt
  ```

  Set `REMOVE_DUPLICATES` to `true`, in order to remove the duplicates instead of just mark them.

* Here is a little discussion on when to merge reads from different lanes to be more accuate. In this project, they will be merged during this step. https://gatk.broadinstitute.org/hc/en-us/articles/360035889471-How-should-I-pre-process-data-from-multiplexed-sequencing-and-multi-library-designs-

```{bash}
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar
java -jar build/libs/picard.jar
java -jar build/libs/picard.jar MarkDuplicates -INPUT normal-aln-pe-L001.sorted.bam -INPUT normal-aln-pe-L002.sorted.bam -OUTPUT normal.merged.rmdup.bam -METRICS_FILE normal.merged.rmdup.txt -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT
java -jar build/libs/picard.jar MarkDuplicates -INPUT tumor-aln-pe-L001.sorted.bam -INPUT tumor-aln-pe-L002.sorted.bam -OUTPUT tumor.merged.rmdup.bam -METRICS_FILE tumor.merged.rmdup.txt -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT
```

* Till here, `normal.merged.rmdup.bam` and `tumor.merged.rmdup.bam` will be generated. Index those file for further analysis (Somatic mutation detection)

```
samtools faidx hg38.fa
samtools index tumor.merged.rmdup.bam
samtools index normal.merged.rmdup.bam
```



### Set up HLAminer for HLA-typing

* HLAminer is a pipeline for predicting HLA from shotgun sequence data. In this project, `HPTASRwgs_classI.sh` is used for HLA Predictions by Targeted Assembly of Shotgun Reads (longer but more accurate. 
* perl module `Bio::SearchIO` must be installed to use HPTASR
* Edit the fullpath location of bwa and other software dependencies in the shell scripts in the ./bin/ folder, as needed. eg, `vi parseXMLblast.pl` to change the first line to `/usr/bin/perl`
* Note: For HLAminer-1.4, user needs to update `HLA_ABC_GEN.fasta` by running the `updateHLAgenomic.sh` in ./database, and change the bwa directory in that script. Otherwise, it will result `[blastall] FATAL ERROR: contig1|size437|read60|cov13.73: Database ../database/HLA_ABC_GEN.fasta was not found or does not exist`
* Important: Place the fullpath location of normal and tumor fasta files in the "patient.fof" file, and make sure that `patient.fof` and `ncbiBlastConfig.txt` are in the working directory.

```#Instal HLAminer
wget https://www.bcgsc.ca/platform/bioinfo/software/hlaminer/releases/1.4/HLAminer_1-4.tar.gz
tar xzvf HLAminer-1.4.tar.gz
cd HLAminer-1.4/HLAminer_v1.4/

./HPTASRwgs_classI.sh
```

* Better: run the test demo to make sure everything is installed correctly. (To save your life since actual running really takes a long time:/)
  1. Copy ./test-demo/ eg. cp -rf test-demo foo
  2. In folder "foo", edit the patient.fof file to point to your NGS RNAseq data. Ensure all paths are ok.
  3. For HLA Predictions by Targeted Assembly of Shotgun Reads: execute ./HLAminer/foo/HPTASRrnaseq.sh 



### Somatic mutation detection by Strelka

* Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs

* Strelka is run in two steps: (1) configuration (specifying input data and options) and (2) workflow execution (specifying parameters on how strelka is executed). 

* Sample usage:

  ```bash
  # configuration
  ${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
      --normalBam normal.bam \
      --tumorBam tumor.bam \
      --referenceFasta hg38.fa \
      --runDir demo_somatic
  # execution on a single local machine
  demo_somatic/runWorkflow.py -m local -j 8
  ```

  * `--exome` is added for exome and amplicon inputs 


```{bash}
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2
# run demo to check successful installation
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash

~/510Final/strelka/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam normal.merged.rmdup.bam --tumorBam tumor.merged.rmdup.bam --referenceFasta hg38.fa --runDir ~/510Final/strelka/ --exome
python /home/ninasong/510Final/strelka/runWorkflow.py -m local -j 8
```

* After using Skrelka, `somatic.snvs.vcf` and `somatic.indels.vcf` are generated. "PASS" filter is applied to simply pass out the low quality reads by `awk`

* Both PASS.cvf are loaded to IGV to visualize the results of mutation detection. In this project, COLO829 represent hallmark mutation types including amplified *BRAF* V600E and a *CDK2NA* small deletion. After checking on IGV,

  BRAF V600E: Both BAM and VCF show the SNV on ch7:140753336 which will change amino acid 600 from Valine to Glutamic acid.

  CDK2NA del: Both BAM and VCF show the 2-bp deletion at ch9:21970953. 

```{bash}
awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' somatic.snvs.vcf > somatic.snvs.PASS.vcf
#sort the file
grep '^#' somatic.snvs.PASS.vcf > somatic.snvs.PASS.sorted.vcf && grep -v '^#' somatic.snvs.PASS.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> somatic.snvs.PASS.sorted.vcf
#merge both snvs and indel
java -jar $GATK MergeVcfs -I somatic.snvs.PASS.vcf - I somatic.indels.PASS.vcf - O somatic.merged.PASS.vcf.gz
```

* Use GATK MergeVcfs to merge both PASS files. (edit the ‘.bashrc’ file to include the path where GATK is located, thus, java -jar $GATK can be simply used in the further)



### Anotation of Variants Using VEP

* Conda is used here to emulate a Python environment. Our server has Python 3.8.3 as well as Python 2.7.5 installed, those two version will conflict with some of the features of VEP and pavctools in the following steps (after failed million times I got this conclusion;/). After testing, Python version 3.5-3.7 show no conflict with those feature.  The following code shows how to install Conda and use it to set up environment.

```bash
wget https://repo.continuum.io/archive/Anaconda3-2020.07-Linux-x86_64.sh
bash Anaconda3-2020.07-Linux-x86_64.sh
source ~/.bashrc
conda create -n pVACseq python=3.5
conda deactivate #to exit this version as need
```

*  **VEP** (Variant Effect Predictor) predicts the functional effects of genomic variants. In addition to the standard VEP annotations, pVAC-Seq also requires the annotations provided by the Downstream and Wildtype VEP plugins. In order to use the `plugin` feature in VEP, pvacseq needs to be installed first. 
*  During the step of `perl INSTALL.pl`, it may ask user to customize `cache` files for further usage. In this project, "homo sapiens" should be chose with the same version of VEP (eg, homo_sapiens_vep_98_GRCh38 ) 
*  MHCflurry and IEDB binding prediction tools are also installed in this step.


```{bash}
#install pvactools and pvacseq first
pip3 install pvactools
pip3 install pvacseq
pip list
#to check if installed both pvacseq and pvactools successfully
pvacseq/pvactools run -h

#download and install the VEP command line tool
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl

#In this project I chose VEP version 98 instead of 101 (the newest version) since this is the only one could be installed on our server with successfully running cache for anonation
#During perl INSTALL.pl, there will be some options for user to choose. In this example, I used the following chose since I was working on human samples:
y  #Do you want to install any cache files (y/n)?   
y  #Cache directory /home/xiyuliu/.vep does not exists - do you want to create it (y/n)?   
320  #The following species/files are available; which do you want (can specify multiple separated by spaces or 0 for all)?   
 - downloading ftp://ftp.ensembl.org/pub/release-98/variation/indexed_vep_cache/homo_sapiens_vep_98_GRCh38.tar.gz
 - unpacking homo_sapiens_vep_98_GRCh38.tar.gz
y  #Do you want to install any FASTA files (y/n)?   
85  #FASTA files for the following species are available; which do you want (can specify multiple separated by spaces, "0" to install for species specified for cache download)?  
 - downloading Homo_sapiens.GRCh38.dna.toplevel.fa.gz
 - downloading Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai
 - downloading Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi
y  #Do you want to install any plugins (y/n)?   
y  #Plugins directory ../.vep/Plugins does not exists - do you want to create it (y/n)?   
0  #The following plugins are available; which do you want (can specify multiple separated by spaces or 0 for all)?   


#check if cache is now working for annotation
./vep -i examples/homo_sapiens_GRCh38.vcf --cache
git clone https://github.com/Ensembl/VEP_plugins.git
pvacseq install_vep_plugin ~/510Final/VEP/VEP_plugins

#Install MHCflurry which require the same environment as pvactools
pip install mhcflurry
mhcflurry-downloads fetch
#check if mhc installed successfully
mhcflurry-predict -h

#Install IEDB
#Use the command 'python src/predict_binding.py' to get started
wget https://downloads.iedb.org/tools/mhci/3.0/IEDB_MHC_I-3.0.tar.gz
tar -zxvf IEDB_MHC_I-3.0.tar.gz
cd mhc_i
./configure

wget https://downloads.iedb.org/tools/mhcii/3.1/IEDB_MHC_II-3.1.tar.gz
tar -zxvf IEDB_MHC_II-3.1.tar.gz
cd mhc_ii
./configure.py
```

* Sample Usage: 

  ```bash
  perl variant_effect_predictor.pl \
  --input_file <input VCF> --format vcf --output_file <output VCF> \
  --vcf --symbol --terms SO --plugin Downstream --plugin Wildtype \
  [--dir_plugins <VEP_plugins directory>]
  ```


```{bash}
#copy the require file into a certain direction before operation
find /home/ninasong/ -name Wildtype.pm
cp /home/ninasong/510Final/VEP/VEP_plugins/Wildtype.pm /home/ninasong/.vep/Plugins/
#run VEP
vep \
--input_file ~/510Final/strelka/results/variants/somatic.merged.PASS.vcf \
--output_file ~/510Final/VEP/ensembl-vep/somatic.annotated.vcf \
--format vcf --vcf --symbol --terms SO --tsl \
--hgvs --fasta ~/510Final/strelka/hg38.fa \
--offline --cache \
--plugin Downstream --plugin Wildtype
```

* After this step, `somatic.annotated.vcf` will be a result of the annotations provided by the Downstream and Wildtype VEP plugins.



### Adding coverage data to your VCF

* pVACseq is able to parse coverage information directly from the VCF. The expected annotation format is outlined below.

```{bash}
#install vt
#vt is a variant tool set that discovers short variants from Next Generation Sequencing data
git clone https://github.com/atks/vt.git  
cd vt 
make 
make test

#Building and installing CMake without root
Since we will be installing things to your home directory, you should add $HOME/bin to your PATH if you haven’t already. This ensures that your shell knows where to look for binaries (cmake, kallisto, etc.). For this current session, run the following from your terminal:

export PATH=$HOME/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lib/:$LD_LIBRARY_PATH
Afterwards, place the same code into your shell startup file (e.g. one of ~/.bashrc, ~/.zshrc, etc.).

We’ll assume all downloads are in ~/Downloads. Make sure to change to the downloads directory before downloading each one of the source archives (*.tar.gz files)

#Download the latest “Unix/Linux Source” *.tar.gz file from CMake downloads page 
rsync -avz /Users/ninasong/Downloads/cmake-3.19.0-rc2.tar.gz ninasong@trgn.usc.edu:~/510Final/VEP/cmake 
tar -xf cmake*.tar.gz
cd cmake*
./configure --prefix=$HOME
make
make install

#download bam-readcount
git clone --recursive git://github.com/genome/bam-readcount.git
cd bam-readcount/
cmake ~/510Final/VEP/bam-readcount/
make
export PATH=~/510Final/VEP/bam-readcount/bin:$PATH
bam-readcount

```

## Perform epitope prediction
* At first , the results showed "The TSV file is empty. Please check that the input VCF contains missense, inframe indel, or frameshift mutations."
* For VCFs created by Strelka which does not output a GT field for its calls. Try to modify with vcf-genotype-annotator to vcf-genotype-annotator.When running this tool specify which GT value to use (0/1) and it will append another sample to the VCF with that genotype. Then we can then use the new dummy sample name when running pVACseq and it will process all the variants in the VCF.
```
vcf-genotype-annotator ~/510Final/pvacseq/test/somatic.annotated.genotype.vcf TUMOR {0/1}
```

* Run epitope prediction software
```
pvacseq run \
~/510Final/pvacseq/test/somatic.annotated.genotype.vcf \
TUMOR \
HLA-A*01:01,HLA-B*40:02,HLA-C*03:03 \
MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
~/510Final/pvacseq/test/prediction_1 \
-e 8,9,10 \
--normal-sample-name NORMAL \
--top-score-metric=lowest -d full 
```
