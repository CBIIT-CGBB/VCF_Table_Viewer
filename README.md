# Runx1 WES Variant Calling

## Pipeline Overview
### Germline calling
* GATK haplotypecaller
* GATK joint calling
* variant recalibration (VQSR)
* freebayes, LoFreq

### Germline filtering
* PASS + QUAL > 100
* freebayes: (QUAL>20, DP>8, QUAL/AO>10, SAF>0, SAR>0, RPR>1, RPL>1)
* LoFreq: kept if not in haplotypecaller or freebayes (?)   
* After haplotypecaller/freebayes/LoFreq merging:
    - region filtering
    - vcf normalization
    - Population filtering:
        - 1000G VAF < 1%
        - ExAC VAF < 1%
        - gnomAD3 popmax < 1%
* Manual review in IGVtools
* Re-evaluated depth, VAF, DP4 counts with reads seq length >= 80, no secondary alignment, mapping quality > 50

### Germline annotation
* SnpEff
* SnpSift

### Somatic calling
* Mutect2, Manta + Strelka, LoFreq, MuSE
* Compare with germline, produce combined list of variants
* For single sample: Use Mutect2 with PoN
    - if any variant found in unaffected family members -> germline

### Somatic filtering
* region filtering
    - "if the surrounding sequence met the criteria like: BBBBB, ABABABAB,
       ABCABCABC, ABCDABCD, mutations will be removed"  (Not sure what this means...)
* "depth >=10, VAF>=0.03, and alt-allele-reads>=3" in somatic sample, 
   and "depth >=8, VAF<0.01" in fibroblast control."
* at least 1 supporting read from each strand
* collect all variants in CHIP and leukemia driver genes regardless of filtering
    - if VAF>0.01 and alt reads >= 3 and good alignment -> keep

### Somatic annotation

### Prediction of deleterious variants
* Germline (Likely)Pathogenic in CLINVAR
    - include Fanconi anemia variants of uncertain or conflicting pathogenicity
* Nonsense and frameshift variants
* Splice donor/acceptor/region variants
* Splice AI DS_AG/DS_AL/DS_DG/DS_DL > 0.5
* 4 of 6 predictors predict deleterious outcomes
    - CADD_phred > 20
    - SIFT_pred='D'
    - FATHMM_pred='D'
    - Polyphen2_HDIV_pred='D' or 'P'
    - MutationTaster_pred='D' or 'A'
    - MetaSVM_pred='D'

    Can use REVEL, an ensemble method that uses scores from MutPred, FATHMM v2.3, VEST 3.0, PolyPhen-2, SIFT, PROVEAN, MutationAssessor, MutationTaster, LRT, GERP++, SiPhy, phyloP, and phastCons. REVEL score above 0.5, corresponding to a sensitivity of 0.754 and specificity of 0.891. 

    Can use DANN in place of CADD, but not clear what the score cutoff should be (0-1.0).  Higher is more likely to be deleterious.

## Gene Lists
<u>Myeloid Driver Gene List</u>

CSF3R	MPL	    NRAS	PTEN	SMC3	HRAS	WT1	KMT2A	CBL	ETV6	KRAS	PTPN11	FLT3	IDH2	TP53	SRSF2	SETBP1	CALR	JAK3	CEBPA	CBLC	DNMT3A	SF3B1	IDH1	ASXL1	GNAS	RUNX1	U2AF1	MYD88	CBLB	GATA2	PDGFRA	KIT	TET2	FBXW7	NPM1	IKZF1	CUX1	BRAF	EZH2	RAD21	JAK2	CDKN2A	ABL1	NOTCH1	ZRSR2	BCOR	KDM6A	GATA1	SMC1A	ATRX	STAG2	BCORL1	PHF6	ADCY5	ANK2	APOB	SHROOM2	ATP2B3	DST	C5	CACNA1B	CACNA1E	CALR	RUNX1	RUNX1T1	CBFB	CBL	    SCARB1	CD74	LRBA	CEBPA	CHD4	COL12A1	VCAN	DDX11	DNAH5	DNAH9	DNMT3A	DNMT3B	DOCK2	DRD2	DSCAM	GPR183	CELSR3	MEGF8	EGFR	ETV6	EZH2	FLG	FLT1	FLT3	GAS6	GATA2	GJB3	GRID1	GRIK2	GRIK4	GRM3	GRM8	HIVEP1	HNRNPK	TNC	IDH1	IDH2	ITPR3	JAK1	JAK2	JAK3	KCNA4	KCNH2	KCNQ2	KDR	KIT	KRAS	KRT19	MAP1B	MAP2	ADAM11	MEFV	MYC	MYH4	MYO5B	MYOC	NF1	    NPM1	NRAS	NTRK3	DDR2	P2RY2	PKHD1	PPP1R3A	PTCH1	PTPN11	PTPRN	RAD21	RBBP4	RFC3	RYR1	RYR3	SCN1A	SRSF2	TRA2B	SHC1	SI	SLC12A3	STRN	NR2E1	TP53	TUBA3C	TYK2	U2AF1	KDM6A	WT1	USP9X	SMC1A	DYSF	CUL1	STC2	CADPS	EED	SYNGAP1	FCGBP	PRPF4B	CACNA1G	BSN	    TOP3B	PKD2L1	SMC3	DCLK1	MTA2	NRXN3	ABCG2	CROCC	MAGI2	THRAP3	MED12	ZBTB33	EDIL3	SPEG	SEMA3A	SCML2	DLC1	PRPF8	STAG2	PTPRT	AKAP13	SCAF8	DIS3	RIMS1	SPEN	SMG1	HECW1	ATP10B	PSME4	MTUS2	SF3B1	KIAA0240	SUZ12	FLRT2	FKBP8	SETBP1	GIGYF2	LRIT1	DNAI1	FOXP1	C10orf28	CECR2	PCDHB1	COL5A3	PLCE1	WAC	    DDX41	KDM3B	LRP1B	CNTN5	PCDHB18	HYDIN	TET2	DCHS2	TMEM104	BCOR	PHIP	ATG16L1	C10orf118	MTMR8	CACNA2D3	PCDHA13	PCDHA6	KCNK13	NMUR2	VARS2	FAM40B	PLEKHH1	ALPK3	KCNT1	ZNF687	KIAA1529	RNF213	NLRC4	MLL3	LRRC4	SEMA4A	IKZF4	CSMD1	PRAMEF2	FAM65A	DYNC2H1	E2F8	TRPM3	TET1	ABTB1	KIAA1683	EPPK1	CRISPLD1	FAM57B	SYT15	HMCN1	PHF6	PDCD2L	TTBK1	KIF2B	LNX1	TCEAL3	CNTNAP4	NAV1	PKHD1L1	MUC16	PKD1L2	CSMD3	GBP4	ARAP2	ZC3H18	MYOM3	GPR112	KCNU1	FAM47A	TCEAL6	CCDC67	XIRP1	BMPER	ASXL1	CMYA5	UNC5B	OR8B12	PHACTR1	CADM2	NALCN	BOD1L	SLC39A5	KSR2	FAM154B	KIAA1267	EPHA10	FRYL	ILDR1	KRT79	FAM5C	FREM2	OR13H1	FAM70B	GSTK1	GALNTL4	C17orf97	CUEDC1	OR11H12	CLEC18B	MUC5B	MROH5	SBF1P1	MIR142	EEF1A1P29	HSP90B3P	IGHG3	ST13P13

<u>ACMG Gene List</u>

ACTA2	ACTC1	ACVRL1	APC	    APOB	ATP7B	BMPR1A	BRCA1	BRCA2	BTD	    CACNA1S	CASQ2	COL3A1	DSC2	DSG2	DSP	    ENG	    FBN1	FLNC	GAA	    GLA	    HFE	    HNF1A	KCNH2	KCNQ1	LDLR	LMNA	MAX	    MEN1	MLH1	MSH2	MSH6	MUTYH	MYBPC3	MYH11	MYH7	MYL2	MYL3	NF2	    OTC	    PALB2	PCSK9	PKP2	PMS2	PRKAG2	PTEN	RB1	    RET	    RPE65	RYR1	RYR2	SCN5A	SDHAF2	SDHB	SDHC	SDHD	SMAD3	SMAD4	STK11	TGFBR1	TGFBR2	TMEM127	TMEM43	TNNI3	TNNT2	TP53	TPM1	TRDN	TSC1	TSC2	TTN	    VHL	    WT1

<u>Integrative Oncogenomics Cancer Driver Genes</u>

https://www.nature.com/articles/s41568-020-0290-x

568 genes: IntOGen-DriverGenes.tsv
Leukemia subset (61 genes): IntOGen-DriverGenes_AML.tsv

<u>Clinical consequences of clonal hematopoiesis of indeterminate potential (CHIP) gene list</u>

'DNMT3A','TET2','ASXL1',
'TP53','JAK2','SF3B1','GNB1','CBL','SRSF2','PPM1D',
'GNAS','BRCC3','CREBBP','NRAS','RAD21','SETDB1','U2AF1','SETD2',
'BCL11B', 'BCOR', 'BCORL1', 'BIRC3', 'BRAF',
'CARD11', 'CD58', 'CD79B', 'CNOT3', 'CUX1', 'DDX3X', 'EP300', 'ETV6', 'EZH2', 'FAM46C',
'FBXW7', 'FLT3', 'FOXP1', 'HIST1H1C', 'IDH2', 'IZKF1', 'JAK2', 'JARID2', 'KMT2D',
'KDM6A', 'KIT', 'KLHL6', 'KRAS', 'LUC7L2', 'MPL', 'MYD88', 'NOTCH1', 'NOTCH2', 'PDSS2',
'PHF6', 'PIK3CA', 'PRDM1', 'PRPB40B', 'PTPN11', 'RIT1', 'RPS15', 'SETDB1', 'SF3A1',
'SMC1A', 'SMC3', 'STAG1', 'STAG2', 'STAT3', 'SUZ12', 'TBL1XR1', 'TET1', 'TNFAIP3', 'TNFRSF14', 'ZRSR2'

<u>CHIP74 List</u>

(Note: CHIP 74 are genes from TOPMed study.)

ASXL1	ASXL2	BCOR	BCORL1	BRAF	BRCC3	CBL	CBLB	CEBPA	CREBBP	CSF1R	CSF3R	CTCF	CUX1	DNMT3A	EED	EP300	ETNK1	ETV6	EZH2	FLT3	GATA1	GATA2	GATA3	GNA13	GNAS	GNB1	IDH1	IDH2	IKZF1	IKZF2	IKZF3	JAK1	JAK2	JAK3	KDM6A	KIT	KRAS	LUC7L2	KMT2A	KMT2D	MPL	NF1	NPM1	NRAS	PDS5B	PDSS2	PHF6	PHIP	PPM1D	PRPF40B	PRPF8	PTEN	PTPN11	RAD21	RUNX1	SETBP1	SETD2	SETDB1	SF1	SF3A1	SF3B1	SRSF2	SMC1A	SMC3	STAG1	STAG2	SUZ12	TET2	TP53	U2AF1	U2AF2	WT1	ZRSR2

<u>Combine lists in Final.ipynb:</u>

sensilist = set(chip74+leudrive)  
len(sensilist)  
mllist = set(list(sensilist)+myeloidpanel)  

## Sample Info

ss = somsamplesheet('../snakemake_wes_somatic/samplesheet.all.tsv')

seqdb = pd.read_csv('../../../00.data/Seq.db.tsv',sep='\t')

```python
notenrolled = ['RISC_PL_1','RISC_PL_2','EMKFPD01','EMKFPD02','EMKFPD06']+\
            ["FPD_2203","FPD_2404","FPD_8703","FPD_8957"] + \
            ["FPD_7194","FPD_5109"] +["FPD_2588"]# Not enrolled | S424= | L56S | p.Gln335His

transplants= {'FPD_0876', 'FPD_5002', 'FPD_0151', 'FPD_5718', 'FPD_7759','FPD_6070','FPD_9068','FPD_9168','FPD_2657', 'FPD_9068'}


malignantlist = ['FPD_0850','FPD_6070','FPD_5718','FPD_7399','FPD_7194','FPD_7759','FPD_9169','FPD_0028', 'FPD_0151','FPD_7658',
                'FPD_5002', 'FPD_0876', 'FPD_9068', 'FPD_3971','FPD_6256','FPD_8898','FPD_6070',
                'FPD_5534', 'FPD_1430', 'FPD_9168','FPD_2657']
# FPD_9236 MGUS ?????? # conclusion: not malignancy

transplants += ['9168_BM221','FPD2657_PB211','FPD9068_BM221']
transplants=set(transplants)
print('Transplant:',transplants)
```
Transplant: {'FPD_7759_BM3', 'FPD_0151_PB1', 'FPD9068_BM221', 'EMKFPD12', 'FPD_5002_BM1A', '9168_BM221', 'FPD2657_PB211', 'FPD_7759A_BM', 'FPD_7759_BMB', 'FPD_9068A_BM'}


## Link Individuals
```python
id2fid = pd.read_csv('link.individual.tsv',sep='\t',index_col=0)
id2fid['fam'] = [i.split('.')[0] for i in id2fid['fid']]
id2fid
```

## Load annotations

```
def loadanno(hdf):
    clndf = pd.read_csv('clinic.rev.txt', sep='\t',index_col=0).fillna('.')
    clndf = clndf.loc[hdf.columns]
```

## Load mutation tables
Somatic:
    - skip not enrolled
    - skip transplants
    - read likely for those with no Skin fibroblast sample
    - read true for those with Skin fibroblast sample
Germline:
    - skip not enrolled

## Maftools

Manually create .vcf files, then convert to .maf using vcf2maf.pl, then manually merge the maf files in cohorts:

fpd 129
fpdmalig 16
fpdnonmalig 113
ctrl 0
truenmsomfpd 62
truemsomfpd 11
true 73

/data/RUNX1/project/1.FPDMM/01.WES/current/Genotype/maf/Maftools:
cli.ctrl.txt      cli.fpdnonmalig.txt  cli.truenmsomfpd.txt  merge.ctrl.maf  merge.fpdmalig.maf     merge.true.maf
cli.fpdmalig.txt  cli.fpd.txt          cli.true.txt          merge.fpd.maf   merge.fpdnonmalig.maf  merge.truenmsomfpd.maf

It says "Go to mafV2/Maftools.ipynb", but can't find this file.

## SpliceAI
Add SpliceAI annotation

!module load snpEff &&
nohup  snpSift annotate /fdb/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz input.sorted.vcf > input.sorted.snv.vcf 2>/dev/null &
nohup  snpSift annotate /fdb/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz input.sorted.vcf > input.sorted.indel.vcf 2>/dev/null & 

DP4:
forward reference count, reverse reference count, forward alternate count and reverse alternate count.

## Pick germline cancer driver gene
Following variant calling, rare variants were further enriched by the application of three-steps: 
    i) Variant with MAF < 1% in the gnomAD 
    ii) Variants class, including missense, protein-truncating, and regulatory;
    iii) Mutation effects, i.e., variant results in protein truncation 
         predicted to be deleterious from 4/6 prediction tools or Path/LP/VUS in clinvar
        a. CADD>20
        b. SIFT_pred: D
        c. FATHMM_pred: D
        d. Polyphen2_HDIV_pred: D/P
        e. MutationTaster_pred: A/D
        f. MetaSVM_pred: D#

Convoluted code for getting predictions, selecting variants with 4/6 deleterious

Find out recurrent germline variants happened in how many families' FPD patients

## Annotate genes
using Harmonizome (/data/RUNX1/data/Harmonizome/checkall/glist.dbmining.tsv)

## List of .tsv output files in 01.WES/current/Genotype
```shell
[sierkml@helix Genotype]$ wc -l *.tsv
        276 BinaryMatrix.tsv
         61 IntOGen-DriverGenes_AML.tsv
        568 IntOGen-DriverGenes.tsv
        293 link.individual.tsv
        242 link.samples.tsv
        171 Log.multifammut.tsv
     305277 Masterdf.addfreebayes.tsv
     367028 Masterdf.addlofreq.tsv
     546140 Masterdf.CallableListed.tsv
     603870 Masterdf.Full.tsv
      42902 Masterdf.GeneSigFinal.tsv
      59027 Masterdf.InGeneSig.tsv
     100351 Masterdf.InGene.tsv
          1 Masterdf.LowVafSigFinal.tsv
    1165313 Masterdf.Merge.tsv
     546140 Masterdf.Normalized.tsv
     239184 Masterdf.PopFilted.tsv
      46111 Masterdf.Prefinal.af003.tsv
         41 Masterdf.Prefinal.lowvaf.tsv
      50007 Masterdf.Prefinal.tsv
     493010 Masterdf.QC.tsv
         89 newsamples.tsv
       6677 PredD.all.tsv
        143 PredD.HemeCaFamRatio.tsv
        127 PredD.sensi.tsv
       1429 Somatic.All.tsv
         87 Somatic.chip74.tsv
         99 Somatic.chipleudrive.tsv
     153422 SpliceAI.output.tsv
```

## Blacklist and whitelist regions

Simplify the below to just use GIAB regions: https://github.com/genome-in-a-bottle/genome-stratifications?tab=readme-ov-file

Also SEQC2: https://sites.google.com/view/seqc2/home/identify-high-confidence-regions

<u>Whitelist</u>
1. IDT Exome targets + 150 bp flank
2. kgmappable (?)
3. genome-in-a-bottle (GIAB) mappable

<u>Blacklist</u>
1. Encode Data Analysis Center (DAC) blacklist
2. Duke blacklist
3. Low Complexity Region (LCR) blacklist
4. Homopolymer blacklist
5. GEM mappable 

Notes: 
    - need to check on overlap of these lists - Duke/DAC and homopolymer/LCR should be the same/similar?
    - GEM mappable comes from old code that is not updated any more, should probably avoid

```python
# White region
wdic = {'exomebed': '/data/RUNX1/pipeline/wxs_phase1/ref/IDT_xGen_Exome_Research_Panel/hg38/xGen_Exome_Research_Panel.targets.withflank150.hg38.bed',
        'kgmappable': '/data/RUNX1/data/MutationCallable/hg38/source/20120824_strict_mask.hg38.bed',
        'giabmappable': '/data/RUNX1/data/MutationCallable/hg38/source/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed'
       } # Black region
bdic = {'blackdac': '/data/RUNX1/data/GMultiMap/Blacklist.DAC/hg38.ENCFF419RSJ.bed',
        'blackduke': '/data/RUNX1/data/GMultiMap/Blacklist.Duke/lists/hg38-blacklist.v2.bed',
        'blacklcr': '/data/RUNX1/data/MutationCallable/hg38/source/LCR-hs38.bed',
        'homopolymer': '/data/RUNX1/data/Homopolymer/Homo_sapiens_assembly38.fasta/Merge.Homopolymer.min7.ext1.bed',
        'mappable': '/data/RUNX1/data/GMultiMap/1.GEM.mappable/hg38.150.mappability.mappability.gz.multihit.bed'
       } # White region
wbeddic = {i: pybedtools.example_bedtool(wdic[i]) for i in wdic.keys()}
bbeddic = {i: pybedtools.example_bedtool(bdic[i]) for i in bdic.keys()}

bedmerge = (wbeddic['exomebed']+wbeddic['kgmappable'].cat(wbeddic['giabmappable']))
print(bedmerge.to_dataframe().shape[0])
for i in bbeddic.keys():
    print(i)
    bedmerge = bedmerge.subtract(bbeddic[i])
    print(bedmerge.to_dataframe().shape[0])
bedmerge.to_dataframe().to_csv('bedmerge.bed',sep='\t',index=None,header=None)
```
Output:
```
163086
blackdac
163086
blackduke
162748
blacklcr
187004
homopolymer
233887
mappable
238232
```
Should be using BedTool.total_coverage(), since we care about the number of bases not number of regions

## Variant Filtering
<u>Tools available for filtering VCF files</u>

* vcftools (https://vcftools.github.io/man_latest.html#OUTPUT%20OPTIONS)
* vcf-prioritise (https://github.com/duncanMR/vcf-prioritise)
    - uses ANNOVAR and VPOT-nf (https://github.com/duncanMR/VPOT-nf)
* vembrane (https://github.com/vembrane/vembrane)
    - python tool for filtering VCFs with conditionals
* PyVCF (https://pyvcf.readthedocs.io/en/latest/)

## Cancer mutation signatures
https://intogen.readthedocs.io/en/v2023/methods_driver_discovery.html

dNdScv
CBaSE
OncodriveCLUSTL
HotMAPS
smRegions
OncodriveFML

Molecular Oncology Almanac
https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-inputs.md

From Snakefile in current/snakemake_wespipe:
vcffilter -f "QUAL > 20 & DP > 8 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" 

