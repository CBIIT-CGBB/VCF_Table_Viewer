library(vcfR)
library(tidyr)
library(dplyr)
library(shiny) #  Shiny web app
library(DT)    #  for data tables
library(data.table) # for melt
library(janitor)
library(stringr)
library(igvShiny)
library(rtracklayer)
library(GenomicAlignments)
library(ggplot2)

##############################################################
####-------------Customization Section--------------------####
##############################################################

# VCF file location (all filtered VCFs are in the same directory)
vcfDir <- "vcfs/filtered/"

# list of individuals for dropdown
individual_list <- read.csv("individual_list_unique.txt", header = F)

# matching of individuals to germline (SK) samples
sample_list <- read.csv("sample_list.txt", header = T, sep="\t")
# Individual	Haplotypecaller
# FPD_0028	FPD_0028_SK211E
# FPD_0151	-
# FPD_0271	FPD_0271_SK221E

# list of important genes to highlight
gene_lists <- read.csv("Gene_lists.txt", header = T, sep = "\t")

# location of bam files for IGV

# On a Mac:
# 1. In Finder, click command-k
# 2. Connect to server 'smb://hpcdrive.nih.gov/nextgen2/sierk/runx/nf/sarek_august2024'
# 3. In Terminal, type: 'ls -l /Volumes'
# 4. Change below to match what is in /Volumes (e.g. "/Volumes/sarek_august2024/preprocessing/recalibrated)
if (Sys.info()['nodename'] == "NCI-02295810-ML") {
  bamDir <- "/Volumes/sierk/runx/nf/sarek/sarek_august2024/preprocessing/recalibrated/"
} else {
  #The Biowulf file system is mounted at /mnt/BW-Data/ on appshare-dev
  bamDir <- "/mnt/BW-Data/recalibrated/"
}

# list of callers for dropdown
callers <- c("haplotypecaller", "deepvariant", "mutect2")

# list of filtering levels
filters <- c("Region", "Population", "Mutation", "ML Driver Genes")
# Region: filter VCF by GIAB mappable region
# Population: filter by population allele frequency < 0.01
# Mutation: filter by significant mutations
# ML Driver Genes: filter by myeloid cancer driver genes

##############################################################

# TODO:
# 1. merge calls from multiple callers into one table
printf <- function(...) print(noquote(sprintf(...)))

# match up individual ID with haplotypecaller samples
hap_dict <- setNames(as.character(sample_list$Haplotypecaller), sample_list$Individual)

addResourcePath("tmpuser", getwd()) # needed for reading in the legend HTML file
# vcf_field_descriptions.html

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
ui <- fluidPage(
                tags$head(
                    tags$style(HTML("table {table-layout: fixed;"))
                ),
                
                # leaving out "Annotated" and "freebayes" for now, too large of files...
                fluidRow(column(4, selectInput("subjectID", "Subject", 
                                     choices = individual_list)),
                         column(4, selectInput("caller", "Caller",
                                               choices = callers)),
                         column(4, selectInput("filterLevel", "Filter",
                                               choices = filters,
                                               selected = "ML Driver Genes"))
                         ),

                fluidRow(column(12, 
                                
                        tabsetPanel(type = "tabs",
                                    tabPanel("Table",
                                      div(dataTableOutput("dataTable"))
                                    ),
                                    tabPanel("Legend",
                                             htmlOutput("legend")
                                    ),
                                    tabPanel("BAM Viewer",
                                             fluidRow(#column(4, actionButton("createIGV", 
                                                      #                       "Generate IGV for selected rows")),
                                                      #column(4, selectInput("variantList",
                                                      #                      "Selected Variants",
                                                      #                      choices = c("test"),
                                                      #                      selected = "test"))
                                                      column(3, selectInput("variantList", "Selected Variants",
                                                                            choices = c(""))),
                                                      column(2, actionButton("removeTracks", "Remove Tracks"),
                                                             align = "left",
                                                             style = "margin-top: 25px;"),
                                                      column(7, htmlOutput("bamfile"))
                                             ),
                                             #div(dataTableOutput("viewer"))
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
                        ) # tabset panel
                ))

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
  
  inputTable <- reactive({
    
    subjID <- input$subjectID
    caller <- input$caller
    filterLevel <- input$filterLevel # annotated, region, population, mutation, driver
    # Note (11/22/24): removing full annotated vcf because it is too big and takes forever to load in the app
  
    
    # if (filterLevel == "Annotated") {
    #   inputDir <- paste0("./vcfs/annotation/",caller, "/")
    # } else {
    inputDir <- paste0(vcfDir, caller, "/")
    #}
    inFile <- "NULL"
    
    # .filtered_snpEff_VEP.ann.rtgfilt.popfilt.sigmut.vcf.gz
    #if (filterLevel == "Annotated") {
    #  filterName <- "ann"
    #} else 
    
    if (filterLevel == "Region") {
      filterName <- "ann.rtgfilt"
    } else if (filterLevel == "Population") {
      filterName <- "ann.rtgfilt.popfilt"
    } else if (filterLevel == "Mutation") {
      filterName <- "ann.rtgfilt.popfilt.sigmut"
    } else if (filterLevel == "ML Driver Genes") {
      filterName <- "ann.rtgfilt.allgenes"
    }
        
    if (caller == "haplotypecaller") {
      if (hap_dict[subjID] != "-") {
        inFile <- paste0(inputDir, hap_dict[subjID], ".", caller, ".filtered_snpEff_VEP.", filterName, ".vcf.gz")
      }
    } else if (caller == "mutect2") { 
      inFile <- paste0(inputDir, subjID, ".", caller, ".filtered_snpEff_VEP.", filterName, ".vcf.gz")
    } else if (caller == "freebayes") {
      inFile <- paste0(inputDir, subjID, ".", caller, ".filtered_snpEff_VEP.", filterName, ".vcf.gz")
    } else if (caller == "deepvariant") {
      inFile <- paste0(inputDir, hap_dict[subjID], ".", caller, ".filtered_snpEff_VEP.", filterName, ".vcf.gz")
    }
    
    print(inFile)
    req(inFile)
    vcf=read.vcfR(inFile, checkFile = TRUE)
    
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
    my.vcf.df <- cbind(as.data.frame(getFIX(vcf)), vcf@gt, INFO2df(vcf))
    
    # split the INFO annotations into columns
    cols = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|
    Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|
    SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|
    gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|
    gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|FREQS|
    CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|CADD_phred|DANN_score|ExAC|FATHMM_pred|
    Interpro_domain|LRT_pred|MetaSVM_pred|MutationTaster_pred|PROVEAN_pred|Polyphen2_HDIV_pred|Polyphen2_HVAR_pred|PrimateAI_pred|REVEL_score|SIFT_pred|dbSNP"
    cols <- gsub(pattern='\\n',replacement="",x=cols)
    cols <- gsub(pattern='\\s',replacement="",x=cols)
    newcols <- strsplit(cols, "\\|")
    newcols <- unlist(newcols)

    my.vcf.ANN.df <- my.vcf.df %>% separate(ANN, sep="\\|", into = newcols)

    # get driver genes and ACMG genes
    my.vcf.ANN.df$Leudrive <- ifelse(my.vcf.ANN.df$SYMBOL %in% gene_lists$Leudrive, 'Y', 'N')
    my.vcf.ANN.df$ACMG <- ifelse(my.vcf.ANN.df$SYMBOL %in% gene_lists$ACMG, 'Y', 'N')
    
    # create an index from chr,pos,ref,alt
    my.vcf.ANN.df$index <- paste(my.vcf.ANN.df$CHROM, my.vcf.ANN.df$POS, my.vcf.ANN.df$REF, my.vcf.ANN.df$ALT,sep='.')
    
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
    
    # select columns to show:
    show_cols_text <- "index,CHROM,POS,REF,ALT,FILTER,AS_FilterStatus,AS_SB_TABLE,DP,GERMQ,
    POPAF,Allele,Consequence,IMPACT,SYMBOL,AA_mut,Leudrive,ACMG,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,
    cDNA_position,CDS_position,Protein_position,Amino_acids,Existing_variation,VARIANT_CLASS,
    APPRIS,CCDS,ENSP,SWISSPROT,SIFT,PolyPhen,miRNA,AF,gnomADe_AF,MAX_AF,FREQS,CLIN_SIG,SOMATIC,
    CADD_phred,DANN_score,FATHMM_pred,LRT_pred,MetaSVM_pred,MT_pred,PROVEAN_pred,
    PP2_HDIV_pred,PP2_HVAR_pred,PrimateAI_pred,REVEL_score"

    # cut down to 36 columns to avoid ajax error on Posit::Connect server
    # show_cols_text_36 <- "index,CHROM,POS,REF,ALT,FILTER,DP,Allele,
    # Consequence,IMPACT,SYMBOL,AA_mut,Leudrive,ACMG,Gene,Feature,
    # Existing_variation,VARIANT_CLASS,
    # SWISSPROT,SIFT,PolyPhen,AF,MAX_AF,CLIN_SIG,
    # CADD_phred,DANN_score,FATHMM_pred,MT_pred,PROVEAN_pred,
    # PrimateAI_pred,REVEL_score"
    
    # replace line returns and spaces
    show_cols_text <- gsub(pattern='\\n',replacement="",x=show_cols_text)
    show_cols_text <- gsub(pattern='\\s',replacement="",x=show_cols_text)
    TEMP <- scan(text=show_cols_text, what="",sep=",") # split into separate items
    show_cols = dput(TEMP) # add quotes
    
    # get the sampleIDs 
    df2 <- df %>% select(starts_with("FPD_"))
    sampleIDs <- names(df2)
    show_cols <- append(show_cols, sampleIDs)
    
    #print(paste("show_cols:", show_cols))
    hide_columns <- which(!(names(df) %in% show_cols)) - 1
  }
  
  #-----------------------------------------------------------------------------
  #  render data table
  #-----------------------------------------------------------------------------
  
  output$dataTable <- renderDT({
    
    df <- inputTable()
    hide_cols <- hide_columns(df)

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
      bamFile <- paste0(bamDir, sampleID, "/", sampleID, ".recal.bam")
      #bamFile <- "/Volumes/sierk/runx/nf/sarek_august2024/preprocessing/recalibrated/FPD_0271_BM221E/FPD_0271_BM221E.recal.bam"
      if (file.exists(bamFile)) {
        output$bamfile <- renderText({paste("loading bam file:", bamFile)})
        x <- readGAlignments(bamFile, param = Rsamtools::ScanBamParam(what="seq", which=GRanges(chrom_pos)))
        loadBamTrackFromLocalData(session, id="igvShiny_0", trackName=sampleID, data=x)
      } else {
        #print(paste0("bam file missing: ", bamFile))
        output$bamfile <- renderText({ paste("<font color=\"#FF0000\"><b>", "bam file missing: ", "</b></font>", bamFile) })
        #output$bamfile <- renderText({ p("some text", style = "color:red") })
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
