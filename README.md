# Hormone-network
The workflow and all the codes used for 'Transcription factor activity in cross-regulation between seedling hormone responses'

## Summary of data
- Plant materials: Three day old etiolated Arabidopsis whole seedlings
1. ChIP-seq data was generated from wild type (Col-0) and transgenetic lines (Col-0 ANAC055::ANAC055-YPet, Col-0 BES1::BES1-YPet, Col-0 EDF1::EDF1-YPet, Col-0 EDF2::EDF2-YPet, Col-0 EDF3::EDF3-YPet, Col-0 EIN3::EIN3-YPet, Col-0 ERF1::ERF1-YPet, Col-0 MYC2::MYC2-YPet, Col-0 MYC3::MYC3-YPet, Col-0 OBP2::OBP2-YPet, Col-0 RAP2.6L::RAP2.6L-YPet, Col-0 STZ::STZ-YPet, Col-0 TCP3::TCP3-YPet and Col-0 TGA5::TGA5-YPet) seedlings treated with air / hormone (BR, ET, JA, SA) for 2 h in this study.
 ```
The raw data for ChIP-seq can be found here:
/data/group/lewseylab/project/lynn/12_own_ChIPdata
```
2. Time series RNA-seq data for BR, SA and SL/KAR was generated from wild type (Col-0) seedlings treated with corresponding hormone for 0, 15 min, 30 min, 1 h, 2 h, 4 h, 8 h, 12 h, and 24 h in this study.
```
The raw data for RNA-seq can be found here:
SA: /data/group/lewseylab/project/lynn/04sa_february1
BR: /data/group/lewseylab/project/lynn/05bl_february1
SL/KAR: /data/group/lewseylab/project/lynn/06col-0_sl_february2
```
3. Time series ET RNA-seq data raw reads were downloaded from Sequence Read Archive (SRA, https://www.ncbi.nlm.nih.gov/sra) with accession number SRA063695.
4. Time series ABA RNA-seq raw reads and ChIP-seq raw reads for ABF1, 3, 4 were downloaded from GEO with accession number GSE80568.
5. Time series JA RNA-seq raw reads were downloaded from GEO with accession number GSE133408.
6. Validation RNA-seq data was generated from wild type (Col-0) and mpk6 mutant seedlings treated with mock (water) / hormone (ABA, ET, JA, SA) for 1 h.
```
The raw data for wild type and mpk6 mutant seedlings can be found here:
/data/group/lewseylab/project/lynn/21_mpk6_RNAseq_data
```
7. Validation proteomics and Phosphoproteomics data was generated from the same samples harvested for generating validation RNA-seq data.
```
The raw data of proteome can be found here:
/data/group/lewseylab/project/lynn/20_Proteomics_Data
```

## Aims
- Evaluation of transcriptome changes and transcriptional regulation during responses to multiple hormones
- Identification of convergence points in hormone cross-regulation network

## Tools used
List of tools used: FastQC V0.11.5, Trimglore V0.4.4, Bowtie2 V2.2.9, samtools V1.3.1, PhantomPeakQualTools v.2.0, MACS V2.1.0, JBrowse, bedtools V2.26.0, ChIPpeakAnno, HiSat2 V2.0.5, Htseq V0.8.0, Salmon v0.8.1, edgeR 3.28.1, Signalling and dynamic regulatory events modelling (SDREM), ClueGO v2.5.7 in Cytoscape v3.8.0, 3D RNA-seq tools, TSIS R package

## Data identifier, methods and original code

### The identifier of newly generated data in this study
1. Time series RNA-seq data: All the raw data (clean data) and processed data (count/TPM matrix.csv for RNA-seq data) have been uploaded to GEO with accession number: GSE182617 
2. Hormone-responsive ChIP-seq data and validation RNA-seq data: All the raw data (clean data) and processed data (count/TPM matrix.csv for RNA-seq data; .narrowPeak for ChIP-seq data) have been uploaded to GEO with accession number: GSE220957
3. Proteomics and Phosphoproteomics data: The mass spectrometry proteomics data (raw and search files) have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier: PXD039958

### Methods - take SA data for example
The summary of experimental design and analysis pipeline can be found in workflow_Hormone_network.pdf

#### RNA-seq data analysis
1.check quality and perform trimming of raw data
```
#produce quality control report and trim the raw reads via trimgalore!
module load fastqc/0.11.5
module load trimgalore/0.4.4
echo "Starting at: $(date)"
trim_galore --phred33 --gzip --fastqc *.fastq.gz
echo "Finished at: $(date)"
```
2.mapping (generate sorted.bam)/quasi-mapping
```
module load hisat-gcc/2.0.5
module load samtools-gcc/1.3.1
echo "Starting at: $(date)"
for i in `cat /data/group/lewseylab/project/lynn/04sa_february1/02cleandata/sample_list.txt`
do
hisat2 -x /data/group/lewseylab/project/lynn/reference_index/Ara_index -U /data/group/lewseylab/project/lynn/04sa_february1/02cleandata/${i}_trimmed.fq.gz --rna-strandness R -S /data/group/lewseylab/project/lynn/04sa_february1/03hisat/${i}.sam
samtools view -S /data/group/lewseylab/project/lynn/04sa_february1/03hisat/${i}.sam -b > /data/group/lewseylab/project/lynn/04sa_february1/03hisat/${i}.bam
samtools sort -n /data/group/lewseylab/project/lynn/04sa_february1/03hisat/${i}.bam -o /data/group/lewseylab/project/lynn/04sa_february1/03hisat/${i}_sorted.bam
done
echo "Finished at: $(date)"

module load salmon-gcc/0.8.1
echo "Starting at: $(date)"
  for i in `cat /data/group/lewseylab/project/lynn/14_SA_AS_data/00_TPM_Quantification/sample_list.txt`
do
echo "processing sample ${i}"
salmon quant -i /data/group/lewseylab/project/lynn/13_JA_AS_data/new_transcriptome_Arabidopsis/AtRTDv2_QUASI_index  -l A \
-r /data/group/lewseylab/project/lynn/04sa_february1/02cleandata/${i}_trimmed.fq.gz \
-p 8 \
--gcBias \
-o /data/group/lewseylab/project/lynn/14_SA_AS_data/00_TPM_Quantification/${i}_quant
done
echo "Finished at: $(date)"
```

3.Count (this step has been done for all hormone datasets but only kept for ET dataset at the end, see workflow_hormone_network.pdf for details)
```
module load htseq/0.8.0
echo "Starting at: $(date)"
for i in `cat /data/group/lewseylab/project/lynn/04sa_february1/02cleandata/sample_list.txt`
do
htseq-count -r name -f bam /data/group/lewseylab/project/lynn/04sa_february1/03hisat/${i}_sorted.bam /data/group/lewseylab/project/lynn/00germ_may17/Annotations/Araport11_GFF3_genes_transposons.201606.gtf > /data/group/lewseylab/project/lynn/04sa_february1/04htseq/${i}.count
done
echo "Finished at: $(date)"
```

#### ChIP-seq data analysis
1.check quality and perform trimming of raw data
```
#fastqc same as RNA-seq analysis

#quality measures for ChIP-seq via phantompeakqualtools
module load samtools-gcc/1.3.1
module load R-gcc7/3.6.0
module load phantompeakqualtools/1.1
echo "Starting at: $(date)"
for i in `cat ./sample_list1.txt`
do
Rscript data/group/lewseylab/project/lynn/13_JA_AS_data/00JA_paper_review/01_histone_frip/phantompeakqualtools/run_spp.R -rf -c=/data/group/lewseylab/project/lynn/12_own_ChIPdata/01_ETH_treatment/01_rap2.6l/${i}.confident.bam  -savp=/data/group/lewseylab/project/lynn/12_own_ChIPdata/01_ETH_treatment/01_rap2.6l/02phantom/${i}.plot.pdf  -out=/data/group/lewseylab/project/lynn/12_own_ChIPdata/01_ETH_treatment/01_rap2.6l/02phantom/${i}.score.txt
done
echo "Finished at: $(date)"
```

2.mapping via bowtie2
```
module load bowtie2-intel/2.2.9
module load samtools-gcc/1.3.1
echo "Starting at: $(date)"
for i in `cat ./sample_list.txt`
do
   bowtie2 -p 8 -x /data/group/lewseylab/project/lynn/reference_index/Bowtie2_ChIp_Ara_index -U ./${i}_trimmed.fq.gz -S ${i}.sam
   samtools view -S ${i}.sam -b > ${i}.bam
   samtools sort ${i}.bam -o ${i}.sorted.bam
   samtools view -b -q 10 ./${i}.sorted.bam > ./${i}.confident.bam
   samtools index ./${i}.confident.bam
   rm ./${i}.sam
done
echo "Finished at: $(date)"
```
3.peak calling
```
module load  macs/2.1.0.20150420
macs2 callpeak -t ./rap2.6l_air_rep1.confident.bam -c ../00_mock/final_mock_air.confident.bam -f BAM -g 1.19e8 -n rap2.6l_air_rep1 -B --outdir ./01final_peaks -q 0.05
```
4.handle biological replicates
```
#Keep only confident peaks--log 10 (q value) >= 15 or q < 10-15
for i in `cat ./sample_list.txt`
do
   awk -F "\t" '{if( $9 >= 15) print $0}' ../${i}_peaks.narrowPeak > ./cut_${i}_peaks.narrowPeak
done

#merge biological replicates
module load  bedtools-gcc/2.26.0
echo "Starting at: $(date)"
 bedtools intersect -a ./cut_rap2.6l_air_rep1_peaks.narrowPeak -b ./cut_rap2.6l_air_rep3_peaks.narrowPeak -f 0.5 -r  > ../01overlap_peaks/50cut_13_airnr.bed
 bedtools intersect -a ./cut_rap2.6l_air_rep1_peaks.narrowPeak -b ./cut_rap2.6l_air_rep4_peaks.narrowPeak -f 0.5 -r  > ../01overlap_peaks/50cut_14_airnr.bed
 bedtools intersect -a ./cut_rap2.6l_air_rep3_peaks.narrowPeak -b ./cut_rap2.6l_air_rep4_peaks.narrowPeak -f 0.5 -r > ../01overlap_peaks/50cut_34_airnr.bed
 cat ../01overlap_peaks/50cut_13_airnr.bed ../01overlap_peaks/50cut_14_airnr.bed ../01overlap_peaks/50cut_34_airnr.bed > ../01overlap_peaks/50cut_rap2.6l_airnr_fv.bed
 bedtools intersect -a ./cut_rap2.6l_eth_rep1_peaks.narrowPeak -b ./cut_rap2.6l_eth_rep4_peaks.narrowPeak -f 0.5 -r  > ../01overlap_peaks/50cut_rap2.6l_ethnr_fv.bed
echo "Finished at: $(date)"
```

### Original code
1. Differentially expressed genes (DEGs) analysis via edgeR (Scripts: edgeR)
2. Peak annotations via ChIPpeakAnno (Scripts: Peak_annotation)
3. Input files for running SDREM (gene expression matrix; TF-target interactions; protein-protein interactions) preparation (Scripts: edgeR (in the last section prepare gene expression inputs for SDREM) & Prepare TF-genes inputs for SDREM & Prepare PPI inputs for SDREM)
4. TF family distribution and enrichment analysis (Scripts: TF_family_enrichment)
5. Hub target genes identification and expression distribution plot generation (Scripts: Hub_target_genes_identification)
6. Gene Ontology enrichment analysis (Scripts: Gene_Ontology_analysis)
7. Differentially alternative spliced (DAS) genes and DEGs start time clarification (Scripts: Calculate_DAS_first_appear)
8. Differentially abundant proteins and phosphopeptides analysis (Scripts: modified_Pipeline_Script & modified_ProteomeDiscoverer_PoissonSeq_Pipeline)

