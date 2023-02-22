# Hormone-network
The workflow and all the codes used for 'Transcription factor activity in cross-regulation between seedling hormone responses'
## Summary of data
- Plant materials: Three day old etiolated Arabidopsis whole seedlings
1. ChIP-seq data was generated from wild type (Col-0) and transgenetic lines (Col-0 ANAC055::ANAC055-YPet, Col-0 BES1::BES1-YPet, Col-0 EDF1::EDF1-YPet, Col-0 EDF2::EDF2-YPet, Col-0 EDF3::EDF3-YPet, Col-0 EIN3::EIN3-YPet, Col-0 ERF1::ERF1-YPet, Col-0 MYC2::MYC2-YPet, Col-0 MYC3::MYC3-YPet, Col-0 OBP2::OBP2-YPet, Col-0 RAP2.6L::RAP2.6L-YPet, Col-0 STZ::STZ-YPet, Col-0 TCP3::TCP3-YPet and Col-0 TGA5::TGA5-YPet) seedlings treated with air / hormone (BR, ET, JA, SA) for 2 h in this study.
 ```
The raw data for ChIP-seq can be found here:
/data/group/lewseylab/project/lynn/12_own_ChIPdata
2. Time series RNA-seq data for BR, SA and SL/KAR was generated from wild type (Col-0) seedlings treated with corresponding hormone for 0, 15 min, 30 min, 1 h, 2 h, 4 h, 8 h, 12 h, and 24 h in this study.
```
The raw data for RNA-seq can be found here:
SA: /data/group/lewseylab/project/lynn/04sa_february1
BR: /data/group/lewseylab/project/lynn/05bl_february1
SL/KAR: /data/group/lewseylab/project/lynn/06col-0_sl_february2
3. Time series ET RNA-seq data raw reads were downloaded from Sequence Read Archive (SRA, https://www.ncbi.nlm.nih.gov/sra) with accession number SRA063695.
4. Time series ABA RNA-seq raw reads and ChIP-seq raw reads for ABF1, 3, 4 were downloaded from GEO with accession number GSE80568.
5. Time series JA RNA-seq raw reads were downloaded from GEO with accession number GSE133408.
6. Validation RNA-seq data was generated from wild type (Col-0) and mpk6 mutant seedlings treated with mock (water) / hormone (ABA, ET, JA, SA) for 1 h.
The raw data for wild type and mpk6 mutant seedlings can be found here:
/data/group/lewseylab/project/lynn/21_mpk6_RNAseq_data
7. Validation proteomics and Phosphoproteomics data was generated from the same samples harvested for generating validation RNA-seq data.
The raw data of proteome can be found here:
/data/group/lewseylab/project/lynn/20_Proteomics_Data
## Aims
- Evaluation of transcriptome changes and transcriptional regulation during responses to multiple hormones
- Identification of convergence points in hormone cross-regulation network
## Tools used
List of tools used: FastQC V0.11.5, Trimglore V0.4.4, Bowtie2 V2.2.9, samtools V1.3.1, PhantomPeakQualTools v.2.0, MACS V2.1.0, JBrowse, bedtools V2.26.0, ChIPpeakAnno, HiSat2 V2.0.5, Htseq V0.8.0, Salmon v0.8.1, edgeR 3.28.1, Signalling and dynamic regulatory events modelling (SDREM), ClueGO v2.5.7 in Cytoscape v3.8.0, 3D RNA-seq tools, TSIS R package
## Data identifier, analysis pipeline and original code
- The identifier of newly generated data in this study
1. Time series RNA-seq data: All the raw data (clean data) and processed data (count/TPM matrix.csv for RNA-seq data) have been uploaded to GEO with accession number: GSE182617 
2. Hormone-responsive ChIP-seq data and validation RNA-seq data: All the raw data (clean data) and processed data (count/TPM matrix.csv for RNA-seq data; .narrowPeak for ChIP-seq data) have been uploaded to GEO with accession number: GSE220957
3. Proteomics and Phosphoproteomics data: The mass spectrometry proteomics data (raw and search files) have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier: PXD039958
- Main analysis pipeline
1. The summary of all analysis processes can be found in workflow_Hormone_network.pdf
- Original code
1. Differentially expressed genes (DEGs) analysis via edgeR (Scripts: edgeR)
2. Peak annotations via ChIPpeakAnno
3. Input files for running SDREM (gene expression matrix; TF-target interactions; protein-protein interactions) preparation (Scripts: edgeR (in the last section prepare gene expression inputs for SDREM) & Prepare TF-genes inputs for SDREM & Prepare PPI inputs for SDREM)
4. TF family distribution and enrichment analysis (Scripts: TF_family_enrichment)
5. Hub target genes identification and expression distribution plot generation (Scripts: Hub_target_genes_identification)
6. Gene Ontology enrichment analysis
7. Differentially alternative spliced (DAS) genes and DEGs start time clarification (Scripts: Calculate_DAS_first_appear)
8. Differentially abundant proteins and phosphopeptides analysis (Scripts: modified_Pipeline_Script & modified_ProteomeDiscoverer_PoissonSeq_Pipeline)

