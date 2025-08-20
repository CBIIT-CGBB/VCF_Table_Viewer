# VCF Table Viewer 2
# April 2025 Run of Sarek pipeline

library(shiny) 
library(shinydashboard)
library(shinyFiles)
library(bslib)
library(fs)

library(tidyr)
library(dplyr)
library(ggplot2)

library(DT)
library(data.table) # for melt
library(stringr)

library(vcfR)
library(igvShiny)
library(rtracklayer)
library(GenomicAlignments)

##############################################################
####-------------Customization Section--------------------####
##############################################################

# needed for reading in the legend HTML file vcf_field_descriptions.html
addResourcePath("tmpuser", getwd()) 

global <- reactiveValues(sarekDir = "./",  
                         vcfDir = "./example_data/vcfs/filtered/", 
                         bamDir = "/mnt/BW-Data/recalibrated/")

### Sarek directory
#sarekDir_default <- "./"

### Samplesheet from sarek run
#sampleSheet_default <- paste0(sarekDir_default, "/samplesheet.csv")

### VCF file location (all filtered VCFs should be in the same directory)
#vcfDir_default <- "./vcfs/filtered/"

### location of bam files for IGV ###
#bamDir_default <- paste0(global$sarekDir,"/preprocessing/recalibrated/")
# if (Sys.info()['nodename'] != "NCI-02295810-ML") {
#   # Posit::Connect server
#   # -> The Biowulf file system is mounted at /mnt/BW-Data/ on appshare-dev
#    <- "/mnt/BW-Data/recalibrated/"
# }

sampleSheet <- read.csv("./example_data/samplesheet.csv",)
sampleSheet <- sampleSheet |> arrange(patient)
subject_list <- unique(sampleSheet$patient)

# lists of important genes to highlight
gene_lists <- read.csv("./example_data/Gene_lists.txt", header = T, sep = "\t")

### list of callers for dropdown
# TODO: get from vcf directory
callers <- c("haplotypecaller", "deepvariant", "mutect2", "strelka", "manta")

### list of filtering levels
# Left out due to size constraints:
#  - Annotation: full annotated VCF produced by sarek
#  - Region: filter VCF by GIAB mappable region
# 
# 1. Population:        filter by population allele frequency < 0.01
# 2. Mutation:          filter by significant mutations
# 3. ML Driver Genes:   filter by myeloid cancer driver genes
# 4. Genes of Interest: filter all genes of interest out of region-filtered VCF
filters <- c("Population", "Mutation", "ML Driver Genes", "Genes of Interest")
filterNames <- c("ann.rtgfilt.popfilt", "ann.rtgfilt.popfilt.sigmut", 
                 "ann.rtgfilt.popfilt.sigmut.genesmut", "ann.rtgfilt.allgenes")
names(filterNames) <- filters

##############################################################

# TODO:
# 1. merge calls from multiple callers into one table
# 2. Add tab with summary/plotting tools:
#   General overview of the WES cohort:
#     - how many patients sequenced,
#     - how many patients with at least one Bone Marrows (BM)
#     - how many have paired skin biopsy/BM
#     - how many have multiple timepoints (e.g. 1 time point- 2-4 time points, > 4)
# 
#   What are the most recurrent somatic mutations in Myeloid Malignancy genes?
#   Are there recurrent somatic mutations in NON-Myeloid genes and what are they?
#   Are there recurrent germline mutations?

printf <- function(...) print(noquote(sprintf(...))) # used with igvShiny


# for highlighting mutation severity columns
color_gradient <- function(dt, column_name, gradient_colors = c("#FF6666", "#DDDDDD")) {
  col_func <- colorRampPalette(gradient_colors)
  dt %>% 
    formatStyle(column_name, 
                backgroundColor = styleEqual(
                  sort(unique(dt$x$data[[column_name]]), decreasing = TRUE),
                  col_func(length(unique(dt$x$data[[column_name]])))
                )
    ) 
}

# user interface 
ui <- dashboardPage(
      
      skin = "purple",
      dashboardHeader(title = "VCF Table Viewer", titleWidth = 300),

      dashboardSidebar(
            collapsed = TRUE,
            width = 300,
            tags$head(tags$style(HTML('
                   /* Sidebar overall padding */
                   .main-sidebar, .left-side {
                   padding-top: 20px;
                   padding: 20px; /* Adjust as needed for padding around content */
                   }'
           ))),
           p("text "), # just used for spacer, otherwise "Input Files" overlaps with header for some reason
           h3("Input files"),
           shinyDirButton("sarek_dir", "Sarek output directory", "Select a folder", style="width:200px"),
           verbatimTextOutput("sarekDir", placeholder = TRUE),
           hr(style = "border-top: 1px solid #ccc; margin: 10px 0;"), # Add a styled horizontal rule
           
           p("Load the samplesheet used to run the sarek pipeline to provide the patients and samples.
             (Default: samplesheet.csv in the current directory. Required columns: patient, sample, status)"),
           fileInput("samplesheet", label = "Select samplesheet (.csv format):", 
                     accept = ".csv"),
           #verbatimTextOutput("stylesheet_file", placeholder = TRUE),
           hr(style = "border-top: 1px solid #ccc; margin: 10px 0;"), # Add a styled horizontal rule
           
           p("Specify the directory with the post-pipeline filtered VCF files."),
           shinyDirButton("vcf_dir", "Select the VCF file directory", "Select a folder", style="width:200px"), # 
           verbatimTextOutput("vcfDir", placeholder = TRUE),
           hr(style = "border-top: 1px solid #ccc; margin: 10px 0;"), # Add a styled horizontal rule
           
           p("Specify the directory with the BAM files. Default: preprocessing/recalibrated in the sarek output directory."),
           shinyDirButton("bam_dir", "Select the BAM file directory", "Select a folder", style="width:200px"),
           verbatimTextOutput("bamDir", placeholder = TRUE),
           hr(style = "border-top: 1px solid #ccc; margin: 10px 0;"), # Add a styled horizontal rule
           
           p("Load the lists of genes to highlight. (Currently equires columns 'Leudrive' and	'ACMG73'"),
           fileInput("gene_list", label = "Select the gene lists (tab-delimited, with headers):")
           #         accept = c(".txt", ".tsv")#,
           

      ), # dashboardSidebar
                
      dashboardBody(

                fluidRow(column(3, selectInput("subjectID", "Subject", 
                                     choices = subject_list)),
                         column(3, selectInput("caller", "Caller",
                                               choices = callers)),
                         column(3, selectInput("filterLevel", "Filter",
                                               choices = filters,
                                               selected = "ML Driver Genes"))
                         ),
                
                fluidRow(column(12, align= "right", htmlOutput("numberOfVariants"))),

                fluidRow(column(12, 
                                
                        tabBox( title = "",
                                id = "main",
                                width = 12,
                                    tabPanel("Table",
                                      div(dataTableOutput("dataTable"))
                                    ),
                                    tabPanel("Legend",
                                             htmlOutput("legend")
                                    ),
                                    tabPanel("BAM Viewer",
                                             fluidRow(
                                                      column(3, selectInput("variantList", "Selected Variants",
                                                                            choices = c(""))),
                                                      column(2, actionButton("removeTracks", "Remove Tracks"),
                                                             align = "left",
                                                             style = "margin-top: 25px;"),
                                                      column(7, htmlOutput("bamfile"))
                                             ),
                                             igvShinyOutput('igvShiny_0')
                                    ),
                                    tabPanel("Plots",
                                             h5("Allele Frequency from Mutect2 VCFs"),
                                             fluidRow(
                                                      column(4, actionButton("createPlot", "Generate plot for selected rows")),
                                                      column(4, downloadButton("downloadPlot", "Download Plot"))),
                                             plotOutput("plot")
                                    )
                                             #htmlOutput("viewer"))
                        ) # tabBox
                )) # fluidRow(column(12,
      ) # dashboard body
) # ui

# server is where all calculations are done, tables are pre-rendered
server <- function(input, output, session) {
  
  # Issues with merging VCFs/variants:
  # 1. Need to keep GT sections from multiple samples, and they may differ between callers
  #.   e.g. haplotypecaller: FPD_0028_FPD_0028_SK211E 0/1:53,41:94:99:1246,0,1560
  #         mutect2:         FPD_0028_FPD_0028_BM211E 1|0:65,4:0.07:69:27,1:29,3:57,4:1|0:178436129_G_GA:178436129:54,11,2,2
  #                          FPD_0028_FPD_0028_SK211E 0|0:80,3:0.018:83:29,3:29,0:72,3:1|0:178436129_G_GA:178436129:51,29,2,1
  #
  # 2. Need to create a column that lists callers, e.g.:
  #   a. read in haplotypecaller vcf
  #   b. read in 2nd caller (e.g. DeepVariant)
  #   c. go through HT df, if index is in DV df, add DV to list of callers
  #   d. repeat for other callers

  # Get the VCF and BAM file directories
  roots=c(wd='.', vol='/Volumes', mnt='/mnt')
  
  get_sarek_dir <- reactive({
    shinyDirChoose(input, 'sarek_dir', roots=roots) 
    output$sarekDir <- renderText(as.character(parseDirPath(roots=roots, input$sarek_dir)))
    
    sarekpath <- as.character(parseDirPath(roots=roots, input$sarek_dir))
    print(sarekpath)
    
    return(sarekpath)
  })
  
  get_vcf_dir <- reactive({
    shinyDirChoose(input, 'vcf_dir', roots=roots)  
    output$vcfDir <- renderText(as.character(parseDirPath(roots=roots, input$vcf_dir)))
      
    vcfpath <- as.character(parseDirPath(roots=roots, input$vcf_dir))
    print(vcfpath)
    
    return(vcfpath)
  })
  
  get_bam_dir <- reactive({
    shinyDirChoose(input, 'bam_dir', roots=roots) 
    output$bamDir <- renderText(as.character(parseDirPath(roots=roots, input$bam_dir)))
      
    bampath <- as.character(parseDirPath(roots=roots, input$bam_dir))
    print(bampath)
      
    return(bampath)
  })
  
  
  # sampleSheet <- reactive({
  #     req(global$sampleSheet_path)
  #     read.csv(global$sampleSheet_path)
  #     #sampleSheet <-ss |> arrange(patient)
  #     #return(sampleSheet)
  # })
  
  observeEvent(input$samplesheet, {
    sampleSheet <- read.csv(input$samplesheet$datapath)
    sampleSheet <- sampleSheet |> arrange(patient)
    updateSelectInput(session, "subjectID", choices = unique(sampleSheet$patient))
  })
  
  observeEvent(input$sarek_dir, {
    global$sarekDir <- input$sarek_dir
  })
  
  observeEvent(input$vcf_dir, {
    global$vcfDir <- input$vcf_dir
  })
  
  #observeEvent(input$bam_dir, {
  #})
    # if (is.integer(input$sarek_dir)) {
    #   sarekpath <- sarekDir_default
    #   output$sarekDir <- renderText(sarekDir_default)
    # } else {
    #   sarekpath <- get_sarek_dir()
    # }
    # 
    # samplesheet_default <- paste0(sarekpath, "/samplesheet.csv")
    # #ss_file <- "./samplesheet_3-31-25.csv"
    # sampleSheet <- read.csv(samplesheet_default)
    # 
    # choices_for_subjectID <- unique(sampleSheet$patient)
    # 
    # updateSelectInput(session, "subjectID", choices = choices_for_subjectID)
    
    ### list of individuals for dropdown
    #individual_list <- sort(unique(sampleSheet$patient))
    
  
  

  
  #-----------------------------------------------------------------------------
  #  generate variant dataframe
  #-----------------------------------------------------------------------------
  inputTable <- reactive({
    
    print("Getting input data frame")
    sarekpath <- get_sarek_dir()
    vcfpath <- get_vcf_dir()
    bampath <- get_bam_dir()
    #req(input$samplesheet)
    
    # if (is.integer(input$sarek_dir)) { # no selection made
    #   sarekpath <- glo
    #   output$sarekDir <- renderText(global$sarekDir)
    # } else {
    #   sarekpath <- get_sarek_dir()
    # }
    #sarekpath <- get_sarek_dir()
    #sampleSheet <- get_samplesheet()
    #vcfpath <- get_vcf_dir()
    
    print(paste("sarekpath:", global$sarekDir))
    print(paste("samplesheet:", global$sampleSheet))
    print(paste("vcfpath:", global$vcfDir))
    
    

    # if (is.integer(input$vcf_dir)) { # no selection made
    #   vcfpath <- vcfDir_default
    #   output$vcfDir <- renderText(vcfDir_default)
    # } else {
    #   vcfpath <- get_vcf_dir()
    # }
    
    subjID <- input$subjectID
    caller <- input$caller
    filterLevel <- input$filterLevel # region, population, mutation, driver, genes of interest
    
    inputDir <- paste0(global$vcfDir, caller) # sep = "/"
    #print(inputDir)

    inFile <- "NULL"
    
    # if (filterLevel == "Region") {
    #   filterName <- "ann.rtgfilt"
    # } else if (filterLevel == "Population") {
    #   filterName <- "ann.rtgfilt.popfilt"
    # } else if (filterLevel == "Mutation") {
    #   filterName <- "ann.rtgfilt.popfilt.sigmut"
    # } else if (filterLevel == "ML Driver Genes") {
    #   filterName <- "ann.rtgfilt.popfilt.sigmut.genesmut"
    # } else if (filterLevel == "Genes of Interest") {
    #   filterName <- "ann.rtgfilt.allgenes"
    # }
    filterName <- filterNames[filterLevel]
    # germline    
    if (caller == "haplotypecaller" || caller == "deepvariant") {
      
      germline_samples <- sampleSheet |> filter(patient == subjID & status == 0)
      
      if (dim(germline_samples)[[1]] != 0) {
        if (caller == "haplotypecaller") {
          inFile <- paste0(inputDir, "/", germline_samples[1,"sample"], ".", caller, ".filtered_snpEff_VEP.", filterName, ".vcf.gz")
        } else {
          inFile <- paste0(inputDir, "/", germline_samples[1,"sample"], ".", caller, "_snpEff_VEP.", filterName, ".vcf.gz")
        }
      }
    # somatic
    } else if (caller == "mutect2") { 
        inFile <- paste0(inputDir, "/", subjID, ".", caller, ".filtered_snpEff_VEP.", filterName, ".vcf.gz")
    
    # need to generate somatic_vs_germline pairs for strelka and manta output
    } else if (caller == "strelka") {
      # FPD_9978_BM241E_vs_FPD_9978_SK241E.strelka.somatic_merged_snpEff_VEP.ann.rtgfilt.allgenes.vcf.gz
        inFile <- paste0(inputDir, "/", subjID, ".", caller, ".filtered_snpEff_VEP.", filterName, ".vcf.gz")
    } else if (caller == "manta") {
      # FPD_9886_PB231E_vs_FPD_9886_SK221E.manta.somatic_sv_snpEff_VEP.ann.rtgfilt.vcf.gz
        inFile <- paste0(inputDir, "/", subjID, ".", caller, ".filtered_snpEff_VEP.", filterName, ".vcf.gz")
    }  
    
    print("Reading vcf file:")
    print(inFile)
    req(inFile)
    vcf=read.vcfR(inFile, checkFile = TRUE)
    
    # check that we have variants
    #  - getFIX returns a character vector if there is 1 variant, dataframe otherwise
    numVariants <- 0
    if (!is.vector(getFIX(vcf))) {
      if (dim(getFIX(vcf))[1] == 0) {
        output$numberOfVariants <- renderText({ paste("Number of variants in VCF: ", as.character(numVariants)) })
        return(NULL) # return null table instead of error message
      } else { numVariants <- dim(getFIX(vcf))[1] }
    } else {
      numVariants <- 1
    }
    # prints number of variants for given selections above the table
    output$numberOfVariants <- renderText({ paste("Number of variants in VCF: ", as.character(numVariants)) })

    
    # mutect2 FORMAT:
    # GT:AD:AF:DP:F1R2:F2R1:FAD:SB    0/1:124,17:0.065:141:47,0:27,2:100,12:51,73,0,17
    # 
    # haplotypecaller FORMAT:
    # GT:AD:DP:GQ:PL  0/1:53,41:94:99:1246,0,1560
    #
    # extract.gt(vcf, "DP", as.numeric=T) gets just the DP field
    # -> extract = F collects the rest of the genotype info except for the named field (e.g. "DP")
    # vcf@gt gets the whole gt field 
    
    # create data frame from VCF (fixed + genotype + INFO)
    if (!is.vector(getFIX(vcf))) {
      fixed <- getFIX(vcf)
    } else {
      fixed <- t(data.frame(getFIX(vcf)))
    }
    info <- INFO2df(vcf)
    my.vcf.df <- cbind(as.data.frame(fixed), vcf@gt, info)
    
    # split the INFO annotations into columns
    #cols = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|
    #Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|
    #SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|
    #gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|
    #gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|FREQS|
    #CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|CADD_phred|DANN_score|ExAC|FATHMM_pred|
    #Interpro_domain|LRT_pred|MetaSVM_pred|MutationTaster_pred|PROVEAN_pred|Polyphen2_HDIV_pred|Polyphen2_HVAR_pred|PrimateAI_pred|REVEL_score|SIFT_pred|dbSNP"
    #cols <- gsub(pattern='\\n',replacement="",x=cols)
    #cols <- gsub(pattern='\\s',replacement="",x=cols)
    #newcols <- strsplit(cols, "\\|")
    #newcols <- unlist(newcols)
    
    # get columns from vcf file directly
    ann <- grepl("ID=ANN", vcf@meta)
    txt <- vcf@meta[ann]
    cols <- strsplit(vcf@meta[ann], "Format: ") # get rid of leading text
    cols2 <- gsub("dbSNP\\\">","dbSNP" , cols[[1]][2]) # get rid of trailing text
    newcols <- strsplit(cols2, "\\|")
    newcols <- unlist(newcols)
    #print(length(newcols))    
    #print(newcols)

    #print(colnames(my.vcf.df))
    #print(dim(my.vcf.df))
    #print(my.vcf.df$ANN)
    
    # deal with multiple Annotations
    my.vcf.ANN.df <- my.vcf.df |> separate_wider_delim(ANN, delim=",", names = c("ANN"), too_many="drop")
    my.vcf.ANN.df <- my.vcf.ANN.df |> separate_wider_delim(ANN, delim = "|", names = newcols,  
                                                           too_many = "error", too_few = "error", 
                                                           names_repair = "universal") # 
    # rename AF columns
    if (caller == "haplotypecaller") {
      my.vcf.ANN.df <- my.vcf.ANN.df |> rename("AF...11" = "AF", "AF...70" = "AF_TG") # requires AF to be in columns 11 and 69
      #write.table(my.vcf.ANN.df, file="my.vcf.ANN.df.txt", quote=F)
    }
    #print(dim(my.vcf.ANN.df))
    #Warning: Expected 93 pieces. Additional pieces discarded in 7 rows [1, 2, 3, 4, 5, 6, 7].

    # get driver genes and ACMG genes
    my.vcf.ANN.df$Leudrive <- ifelse(my.vcf.ANN.df$SYMBOL %in% gene_lists$Leudrive &
                                       my.vcf.ANN.df$SYMBOL != '', 'Y', 'N')
    my.vcf.ANN.df$ACMG <- ifelse(my.vcf.ANN.df$SYMBOL %in% gene_lists$ACMG &
                                   my.vcf.ANN.df$SYMBOL != '', 'Y', 'N')
    
    #print(colnames(my.vcf.ANN.df))
    # create an index from chr,pos,ref,alt
    # Error in $<-.data.frame: replacement has 0 rows, data has 7
    my.vcf.ANN.df$index <- paste(my.vcf.ANN.df$CHROM, my.vcf.ANN.df$POS, my.vcf.ANN.df$REF, my.vcf.ANN.df$ALT, sep='.')
    
    # create aa mutation e.g. I255T
    #  Protein_position 70/393
    #  Amino_acids I/T
    
    my.vcf.ANN.df <- my.vcf.ANN.df %>%
      mutate(prot_pos = str_extract(Protein_position, "(^\\d+)/", group=1)) %>%
      mutate(AA1 = str_extract(Amino_acids, "(^\\w+)/", group=1)) %>% 
      mutate(AA2 = str_extract(Amino_acids, "/(\\w+$)", group=1)) %>%
      mutate(AA_mut = case_when(Consequence == "missense_variant" | Consequence == "frameshift_variant" 
                                ~ paste0(AA1, prot_pos, AA2))) %>%
      select(!prot_pos, !AA1, !AA2)
      
    # Extract 1 value from predictors that give values for each transcript:
    # PROVEAN_pred, PolyPhen2_HDIV_pred, PolyPhen2_HVAR_pred, REVEL_score
    my.vcf.ANN.df <- my.vcf.ANN.df %>% 
                     mutate(MT_pred = str_extract(MutationTaster_pred, "\\w"), .keep="unused", .after="MetaSVM_pred") %>%
                     mutate(FATHMM_pred = str_extract(FATHMM_pred, "\\w"), .keep="unused", .after="DANN_score") %>%
                     mutate(PROVEAN_pred = str_extract(PROVEAN_pred, "\\w"), .keep="unused", .after="MT_pred") %>%
                     mutate(PP2_HDIV_pred = str_extract(Polyphen2_HDIV_pred, "\\w"), .keep="unused", .after="PROVEAN_pred") %>%
                     mutate(PP2_HVAR_pred = str_extract(Polyphen2_HVAR_pred, "\\w"), .keep="unused", .after="PP2_HDIV_pred") %>%
                     mutate(REVEL_score = str_extract(REVEL_score, "\\d*\\.?\\d+"), .keep="unused")
    
    # order the columns logically  
    my.vcf.ANN.df <- my.vcf.ANN.df %>% 
      relocate(c("Allele", "FILTER", "SYMBOL", "AA_mut", "Leudrive", "ACMG", "IMPACT", "Consequence",  "Existing_variation"), .after=QUAL) %>%
      relocate(c("CLIN_SIG", "SIFT", "PolyPhen"), .after=SOMATIC) %>% 
      relocate(c("REVEL_score"), .after=DANN_score)
    my.vcf.ANN.df <- my.vcf.ANN.df %>% relocate(index)
    
    my.vcf.ANN.df <- my.vcf.ANN.df %>% mutate(across(c('DANN_score', 'CADD_phred', 'REVEL_score'), \(x) as.numeric(x))) %>%  
                                       mutate(across(c('DANN_score', 'CADD_phred', 'REVEL_score'), \(x) round(x, 3)))   %>%
                                       mutate_if(is.numeric, ~replace(., is.na(.), 0))

    # cut down # of columns to under 36 to avoid ajax error on Posit::Connect server
    # my.vcf.ANN.df <- my.vcf.ANN.df %>% 
    #   select(index, CHROM, POS, REF, ALT, IMPACT, FILTER, DP, Allele, SYMBOL, AA_mut, Leudrive, ACMG,
    #          Consequence, Gene, Feature, starts_with("FPD"), SWISSPROT, 
    #          VARIANT_CLASS, SIFT, PolyPhen, AF, MAX_AF, CLIN_SIG, 
    #          MT_pred, PROVEAN_pred, PrimateAI_pred, 
    #          FATHMM_pred, DANN_score, REVEL_score, CADD_phred) 
    #LRT_pred, MetaSVM_pred, MT_pred, PROVEAN_pred, PrimateAI_pred, FATHMM_pred, PP2_HDIV_pred, PP2_HVAR_pred, PrimateAI_pred, 
    return(my.vcf.ANN.df)
  })
  
  hide_columns <- function(df) {
    # get column indices for columns to begin the display hidden
    
    # TODO: make this easier to modify for default
    # select columns to show:
    show_cols_text <- "index,CHROM,POS,REF,ALT,FILTER,AS_FilterStatus,AS_SB_TABLE,DP,GERMQ,
    POPAF,Allele,Consequence,IMPACT,SYMBOL,AA_mut,Leudrive,ACMG,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,
    cDNA_position,CDS_position,Protein_position,Amino_acids,Existing_variation,VARIANT_CLASS,
    APPRIS,CCDS,ENSP,SWISSPROT,SIFT,PolyPhen,miRNA,AF,gnomADe_AF,MAX_AF,FREQS,CLIN_SIG,SOMATIC,
    CADD_phred,DANN_score,FATHMM_pred,LRT_pred,MetaSVM_pred,MT_pred,PROVEAN_pred,
    PP2_HDIV_pred,PP2_HVAR_pred,PrimateAI_pred,REVEL_score"
    
    # replace line returns and spaces
    show_cols_text <- gsub(pattern='\\n',replacement="",x=show_cols_text)
    show_cols_text <- gsub(pattern='\\s',replacement="",x=show_cols_text)
    TEMP <- scan(text=show_cols_text, what="", sep=",", quiet=TRUE) # split into separate items
    show_cols = dput(TEMP, file="tmp.txt") # add quotes
    
    # get the sampleIDs 
    df2 <- df %>% select(starts_with("FPD_")) # TODO: make this work for any samples
    sampleIDs <- names(df2)
    show_cols <- append(show_cols, sampleIDs)
    
    #print(paste("show_cols:", show_cols))
    hide_columns <- which(!(names(df) %in% show_cols)) - 1
  }
  
  #-----------------------------------------------------------------------------
  #  render data table
  #-----------------------------------------------------------------------------
  # TODO:
  #   1. Get column headers from VCF header
  #   2. Make selection of samples more generic
  #   3. Make selection of columns to highlight more flexible
  
  output$dataTable <- renderDT({
    
    df <- inputTable()
    
    if (is.null(df)) {
        validate("No variants found, select another filter")
        #output$numberOfVariants <- renderText({ "Number of variants: 0" })
    } #else {
      # Prints the number of variants above the table
      #output$numberOfVariants <- renderText({ paste("Number of variants: ", dim(df)[1]) })
    #}
    
    hide_cols <- hide_columns(df) # list of indices to hide initially

    dt <- DT::datatable(
      df, 
      rownames=FALSE, 
      #escape=FALSE, # need if using HTML
      extensions = c('FixedColumns','FixedHeader','Buttons'), # 'Select' to use datatables select instead of DT's
      #selection = 'none', if using datatables select
      options = list(
        dom = 'Blfrtip',
        fixedHeader=TRUE,
        fixedColumns=TRUE,
        pageLength = 100,
        lengthMenu = list(c(50, 100, -1), c('50','100','All')),
        autoWidth = TRUE,
        scrollX = TRUE,
        columnDefs = list(list(visible = FALSE, targets = hide_cols)
                          #list(targets = 0, width = '50px') # c(0,11) doesn't work for some reason
        ),
        buttons = list('colvis',
                    list(
                      extend = "excel",
                      text = 'Excel',
                      exportOptions = list(rows = '.selected') # only export selected rows
                    )),
        searchCols = list(  # initialize filters on each column
          NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL, #NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,  
          list(search = "PASS") # FILTER
          #strrep(NULL, 40)  # 54 columns
        )
      ),
      class = "display nowrap compact", # style
      filter = "top" # location of column filters
    # format columns
    ) %>% formatStyle('IMPACT',
                      backgroundColor = styleEqual(c("HIGH", "MODERATE"), c('#FA5F55', 'yellow'))
    ) %>% formatStyle(c('MT_pred', 'PROVEAN_pred',
                        'PrimateAI_pred','FATHMM_pred', 'PrimateAI_pred'), # 'LRT_pred','MetaSVM_pred', 'PP2_HDIV_pred', 'PP2_HVAR_pred',
                      backgroundColor = styleEqual(c("D"), c('#FA5F55'))
    ) %>% formatStyle(c('Leudrive','ACMG'),
                      backgroundColor = styleEqual(c("Y"), c('lightgreen'))
    ) %>% formatStyle('SYMBOL','ACMG',
                      fontWeight = styleEqual("Y", "bold")
    ) %>% formatStyle('SYMBOL', 'Leudrive',
                      fontWeight = styleEqual("Y", "bold")
    ) %>% formatStyle('REVEL_score',
                      backgroundColor  = styleInterval(c(0.5), c('white','#FA5F55'))
    ) %>% color_gradient("CADD_phred") %>%
          color_gradient("DANN_score") 
          #color_gradient("REVEL_score", gradient_colors = c("#FF6666", "#DDDDDD", "#6666FF"))
    #} else {
    #}
  }) # renderDT
  
  # explanation of vcf labels in table
  output$legend <- renderUI({
    tags$iframe(
      seamless="seamless",
      src="tmpuser/vcf_field_descriptions.html",
      width=800, 
      height=800)
  })
  
  
  #-----------------------------------------------------------------------------#
  # genome viewer
  #-----------------------------------------------------------------------------#

  observeEvent(input$dataTable_rows_selected, {
    x <- inputTable()[input$dataTable_rows_selected, ]
    
    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)
    
    # Update the selectInput menu
    updateSelectInput(session, "variantList", "Selected Variants", choices = x["index"], selected = c(""))
    
  })
  
  observeEvent(input$variantList, {
    
    req(input$variantList)
    
    if (Sys.info()['nodename'] == "NCI-02295810-ML") {
      # Posit::Connect server
      # -> The Biowulf file system is mounted at /mnt/BW-Data/ on appshare-dev
      bampath <- "/mnt/BW-Data/recalibrated/"
    }
    
    if (Sys.info()['nodename'] == "NCI-02295810-ML") {
      global$bamDir <- "/Volumes/sierk/runx/nf/sarek/sarek_april2025/preprocessing/recalibrated/"
    } 
    output$bamDir <- renderText(global$bamDir)
    
    x <- inputTable()
    variant <- x[x$index == input$variantList, ] # input$variantList
    chrom_pos <- paste0(variant$CHROM, ":", variant$POS)
    showGenomicRegion(session, id="igvShiny_0", chrom_pos) # chr21:10,397,614-10,423,341
    
    samples <- variant %>% select(starts_with("FPD_")) # FPD_0028_FPD_0028_SK211E
    #print(samples)
    sample_names <- colnames(samples)
    
    for (s in sample_names) {
      print(s)
      sampleID <- substr(s,10,24)
      print(paste0("sampleID: ", sampleID))
      
      subjID <- input$subjectID
      bamFile <- paste0(global$bamDir, sampleID, "/", sampleID, ".recal.bam")
      #bamFile <- "/Volumes/sierk/runx/nf/sarek_august2024/preprocessing/recalibrated/FPD_0271_BM221E/FPD_0271_BM221E.recal.bam"
      if (file.exists(bamFile)) {
        output$bamfile <- renderText({paste("loading bam file:", bamFile)})
        x <- readGAlignments(bamFile, param = Rsamtools::ScanBamParam(what="seq", which=GRanges(chrom_pos)))
        loadBamTrackFromLocalData(session, id="igvShiny_0", trackName=sampleID, data=x)
      } else {
        output$bamfile <- renderText({ paste("<font color=\"#FF0000\"><b>", "bam file missing: ", "</b></font>", bamFile) })
      }
    }

    #bamFileDisplay <- paste0(bamDir, "\n", sampleID, "/", sampleID, ".recal.bam")
    #output$bamfile <- renderText(bamFileDisplay)

  }, ignoreInit = TRUE) # don't display until user clicks on dropdown menu
  
  # from igvShinyDemo, potentially useful:
  observeEvent(input$removeTracks, {
    printf("---- removeUserTracks")
    removeUserAddedTracks(session, id="igvShiny_0")
  })
  
  observeEvent(input$igvReady, {
    printf("--- igvReady")
    containerID <- input$igvReady
    printf("igv ready, %s", containerID)
    # loadBedTrack(session, id=containerID, trackName="bed5 loaded on ready", tbl=tbl.bed5, color="red");
  })
  
  observeEvent(input$trackClick, {
    printf("--- trackclick event")
    x <- input$trackClick
    print(x)
  })
  
  observeEvent(input[["igv-trackClick"]], {
    printf("--- igv-trackClick event")
    x <- input[["igv-trackClick"]]
    print(x)
    attribute.name.positions <- grep("name", names(x))
    attribute.value.positions <- grep("value", names(x))
    attribute.names <- as.character(x)[attribute.name.positions]
    attribute.values <- as.character(x)[attribute.value.positions]
    tbl <- data.frame(name=attribute.names,
                      value=attribute.values,
                      stringsAsFactors=FALSE)
    dialogContent <- renderTable(tbl)
    html <- HTML(dialogContent())
    showModal(modalDialog(html, easyClose=TRUE))
  })
  
  observeEvent(input$getChromLocButton, {
    # printf("--- getChromLoc event")
    # sends message to igv.js in browser; currentGenomicRegion.<id> event sent back
    # see below for how that can be captured and displayed
    getGenomicRegion(session, id="igvShiny_0")
  })
  
  observeEvent(input$clearChromLocButton, {
    output$chromLocDisplay <- renderText({" "})
  })
  
  observeEvent(input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]], {
    newLoc <- input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]]
    #printf("new chromLocString: %s", newLoc)
    output$chromLocDisplay <- renderText({newLoc})
  })
  
  output$igvShiny_0 <- renderIgvShiny({
    cat("--- starting renderIgvShiny\n");
    genomeOptions <- parseAndValidateGenomeSpec(genomeName="hg38")
    x <- igvShiny(genomeOptions,
                  displayMode="SQUISHED",
                  tracks=list()
    )
    cat("--- ending renderIgvShiny\n");
    return(x)
  })
  
  #-----------------------------------------------------------------------------#
  
  #-----------------------------------------------------------------------------#
  # plot of allele frequencies
  #-----------------------------------------------------------------------------#
  # TODO:
  #   1. Enable for manta & strelka calls (these are currently not in a merged VCF like mutect2)
  #   2. Enable plotting of other values besides AF
  #   3. Provide error message if attempting to plot nonsupported callers
  
  observeEvent(input$createPlot, {

    subjID <- input$subjectID
    print(paste("Generating plot for", subjID, "..."))

    x <- inputTable()[input$dataTable_rows_selected, ]
    
    # Genotype formatting
    # haplotypecaller: 0/1:53,41:94:99:1178,0,1559
    #                  GT:AD:DP:GQ:PL
    # mutect2: GT:AD:AF:DP:F1R2:F2R1:FAD:SB
    #          0/1:240,8:0.005226:248:58,0:100,0:208,4:135,105,8,0
    #          SB “Per-sample component statistics which comprise the Fisher’s Exact Test to detect strand bias.”
    
    # need gene symbol + aa mutation if available, select out GT fields
    samples <- x %>% mutate(index = case_when(!is.na(AA_mut) ~ paste0(SYMBOL, "(", AA_mut, ")"), 
                                              .default = index)) %>% 
                     select(index, starts_with("FPD_")) 
                     
    print(samples)

    getAF <- ~as.numeric(unlist(strsplit(.x, ":"))[3])
    samples <- samples %>% rowwise() %>% mutate(across(starts_with("FPD_"), getAF)) %>%
                       rename_with(~str_remove(., "^FPD_[\\d]{4}_")) %>%
                       select(contains("SK"), sort(colnames(.)))
    #print(samples)
      
    factor(substring(x, 1, 2)) # orders the sample names
    
    samples_melt <- melt(as.data.table(samples), id = "index")
    #print(samples_melt)
    
    generatePlot <- function() {
      ggplot(data = samples_melt, aes(x = variable, y = value, color = index, group = index)) +
        geom_point(size = 2) + 
        geom_line() + 
        labs(title = paste(subjID, "Mutect2"), x = "Samples", y = "Allele Frequency", color = "Variants") +
        theme(text = element_text(size = 16),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.text = element_text(size=10))
    }
      
    if (input$caller == "mutect2") {
        
        output$plot <- renderPlot({
          generatePlot()  
        }) # renderPlot
        
        #observeEvent(input$downloadPlot, {
        output$downloadPlot <- downloadHandler(
          filename = "allele_freqs.pdf",
          content = function(file) {
            #pdf(file)
            #generatePlot()
            #dev.off()
            ggsave(file, plot=generatePlot(), dpi = 300, 
                   width = 8.5, height = 5.5, units = "in", device="pdf")
            #filename = function(){paste("input$plot3",'.png',sep='')},
            #content = function(file){
            #  ggsave(file,plot=data$plot)
          }
        )
    } # if input$caller == mutect2
    
  }) # observeEvent
  

      
} # server


# run the app
shinyApp(ui, server) # launch.browser = TRUE, options = list(width = 1600)


# all_cols <- c(CHROM,POS,ID,REF,ALT,QUAL,FILTER,AS_FilterStatus,AS_SB_TABLE,AS_UNIQ_ALT_READ_COUNT,CONTQ,DP,ECNT,GERMQ,
#                MBQ,MFRL,MMQ,MPOS,NALOD,NCount,NLOD,OCM,PON,POPAF,ROQ,RPA,RU,SEQQ,STR,STRANDQ,STRQ,TLOD,LOF,NMD,
#                Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,
#                Codons,Existing_variation,DISTANCE,STRAND,FLAGS,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID,CANONICAL,MANE_SELECT,MANE_PLUS_CLINICAL,TSL,APPRIS,CCDS,ENSP,
#                SWISSPROT,TREMBL,UNIPARC,UNIPROT_ISOFORM,GENE_PHENO,SIFT,PolyPhen,DOMAINS,miRNA,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,gnomADe_AF,gnomADe_AFR_AF,
#                gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,gnomADg_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,
#                gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_OTH_AF,gnomADg_SAS_AF,MAX_AF,MAX_AF_POPS,FREQS,
#                CLIN_SIG,SOMATIC,PHENO,PUBMED,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS,CADD_phred,DANN_score,ExAC,FATHMM_pred,
#                Interpro_domain,LRT_pred,MetaSVM_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,PrimateAI_pred,REVEL_score,SIFT_pred,dbSNP)

# You need to add "l" (small letter "L") to dom, that makes Blfrtip:
# B - Buttons
# l - Length changing input control
# f - Filtering input
# r - pRocessing display element
# t - Table
# i - Table information summary
# p - Pagination control
#
