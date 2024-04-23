### FDA / PDA project
### Nanopore data analysis

### Last updated: MLR 11/18/2023

##Step 1: Download fastq files from NCBI
# Make file with all accession numbers and save as SraAccList.txt

#Install SRA Toolkit

PATH=$PATH:/storage/home/mlr355/work/Programs/sratoolkit.3.0.7-ubuntu64/bin

#Move to scratch - move space, but only available for a month
cd /storage/home/mlr355/scratch

#Download data from NCBI - downloads and fastq made by NCBI from fast5 submitted data

prefetch --option-file SraAccList.txt
fasterq-dump SRR*/SRR*.sra

#Move fastq files to folder
mkdir NanoFastq
cp *.fastq NanoFastq

#Remove extra data to save space
rm -r SRR* 

#Rename fastq files from SRR to SampleId for easier tracking
#Make a file named rename.txt with SRR and SampleID names in the first and second column, respectively, separated by a comma

awk -F',' 'system("mv " $1 " " $2)' ../rename.txt

##Step 2: Check quality of Nanopore sequences
#Load conda to install software
module load anaconda3
conda activate bioinf #Made environment to deposit all bioinformatic software

#Quality check of fastq files
pip install nanoplot

NanoPlot -t 4 --huge -o NanoPlotQC --fastq {fastq file} --verbose

#Move NanoPlot results to /work folder
cp -r NanoPlot* ../../work/FDA/Nanopore


##Step 3: Assemble reads into contigs using MetaFlye
#Install Flye assembler
conda install -c bioconda flye

#Use Flye assembler to assemble into contig - increase threads to 8 for faster assembly
cd NanoFastq

for f in *.fastq; do
	name="${f%.*}"
	file_output="${PATH:$PATH/}_MetaFlye"
flye --nano-raw ${f} --meta --out-dir ${file_output} -t 8
done

#Install AGB to visualize assembly graphs
conda create -c almiheenko -c bioconda -n AGB agb

#Run AGB - cannot make it work
conda activate AGB
agb.py --graph <graphviz file> -a Flye

##Step 4: Check quality of the assembled contigs
#Check quality of assembly
conda install -c bioconda quast

python /storage/home/mlr355/.conda/envs/bioinf/bin/metaquast.py {contigs} --max-ref-number 0

##Step 5: Polish the contigs with Medaka - one round of polishing is recommended by Nanopore

#Install medaka polisher
pip install medaka

#Run one polishing round with medaka
medaka_consensus -i *.fastq {basecalls} -d assembly.fa -o medaka_out -t 4 -m r941_min_fast_g303

#Step 6: Check quality of polished contigs
#Check quality of polished assembly
metaquast.py {contigs} --max-ref-number 0

#Rename all assembly and consensus files to include the sample name and move to work directory.
#Save the full MetaFlye and Medaka output in /storage/group/jzk303/default/data/mlr355/FDA



####ADDING STEPS FROM METAGENOMIC DATA ANALYSIS BOOK CHAPTER 12 - PP 235-259 (basing off on the long read analysis 
#### Last updated MLR 12/11/2023
#Step 7: Frame shifting correction using DIAMOND +MEGAN-LR
#DIAMOND is used to perform a frame-shifting aware DNA-to-protein alingment of long read data against NCBI-NR database

#Download reference database from https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nb.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

#Index reference database
diamond makedb --in nr.gz -d nr -p 40

#Align and index long read assemblies to reference database
diamond blastx -d nr.dmnd -q {assembly.fasta} -o {alignment.daa} -f 100 -F 15 --range-culling --top 10 -p 40

#Use meganizer to convert daa file into Megan format
daa-meganizer -i {.daa file} -mdb megan-map-Feb2022.db --longReads

#Open files in Megan GUI

#Compare files: File> Compare
#Select meganized files to compare

#To make taxonomic relative abundance, select Rank (Genus or Species)
#Click Show > Select Stacked barplot, then click on %

#For functional comparison
#Click COG - this will open functional tree of COGs
#To change rank level to subfunction, click Tree>Collapse at Level > "2"







Binning of contigs into MAGs - grouping contigs into draft metagenomes
#Align reads to polished contigs with minimap2
minimap2 -ax map-ont {polished} {fastq} > aligned.sam

#Convert SAM to BAM file
samtools view -bS aligned.sam > aligned.bam

#Sort BAM file 
samtools sort aligned.bam -o sorted_aligned.bam

#Install MetaBAT2
conda install -c bioconda metabat2

#Run MetaBAT2 - check docs
runMetaBat.sh -d -v -m 2000 consensus.fasta {sorted_aligned.bam}

#Step 8: Check quality of MAGs - according to standards established by Genomics Standards Consortium, high-quality MAGs is minimally consistent with holding completeness >90% and <5% contamination. 

#Install CheckM

#Run CheckM on all samples - assuming all samples MAGs are in the same directory
checkm lineage_wf -t 8 -x {bins_dir} {checkm_output}

#Step 9: Taxonomic annotation of MAGs with GTDB-Tk
#Install GTDB-Tk and download the database

#Run Taxonomic classification
gtdbtk classify_wf --genome-dir {bin_dir} --out_dir {gtdbtk_dir} --cpus 8 -x fasta

#Step 10: Functional annotation with Prokka


#Step 11: Quantifying genes of interest present in contigs


#Step 12: Summarizing functional annotation



