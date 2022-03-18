
#title: "HeartExplorer website"
#author: "Neda R. Mehdiabadi"
#date: "18.03.2022"


#loading libraries
library(shinydashboard)
library(Seurat)
library(shinycssloaders)
library(dplyr)
library(ggplot2)
library(fst)
library(edgeR)
library(NMF)
library(shinyjs)
library(tableHTML)
library(htmlwidgets)


#loading data
load("data/datasets.Rdata")

#add logo and name to the header 
title <- tags$img(src = "weblogo.png", height = "50", width = "70", "Heart Explorer")

ui <- dashboardPage(
  dashboardHeader(title = title ,titleWidth = "300px"),
  dashboardSidebar(
    # Remove the sidebar toggle element
    tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
    sidebarMenu(
      menuItem("Overview", tabName = "about"),
      menuItem("Visualization", tabName = "visualization"),
      menuItem("Pathway enrichment analysis", tabName = "enrichment"),
      menuItem("Methods", tabName = "method"),
      menuItem("Data and code availability", tabName = "code")
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML('.skin-blue .main-header .logo {
                              background-color: black;
    }


        /* logo when hovered */
        .skin-blue .main-header .logo:hover {
                              background-color: black;
        }

        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: black;
        }
        .skin-blue .main-sidebar {
                                background-color:white;
                                }
        .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #128FC8;
        }
        .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                background-color: white;
                                color: black;
        }
        .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                                background-color: #D3D3D3;
                                }
        .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color:black;
        }'))),

    tabItems(
      # first tab content
      tabItem(tabName = "about",
              h2("Defining the Fetal Gene Program at Single Cell Resolution in Pediatric Dilated Cardiomyopathy.", align = "center"),
              h3("Authors:"),
              div(p("Neda R. Mehdiabadi, Choon Boon Sim, Belinda Phipson, Ravi K. R. Kalathur, Yuliangzi Sun, Celine J. Vivien, Menno ter Huurne, Adam T. Piers, James E. Hudson, Alicia Oshlack, Robert G. Weintraub, Igor E. Konstantinov, Nathan J. Palpant, David A. Elliott, Enzo R. Porrello"),style="text-align: justify;"),
              h3("Abstract"),
              h3("Background:"),
              div(p("A central dogma in cardiac biology is that the gene expression pattern observed in postnatal heart resemblance to those observed during fetal cardiac development in response to stress. The phenomenon of fetal gene re-activation in heart failure has been traditionally studied in cardiomyocytes, however, the extent to which the fetal gene program is recapitulated in other cardiac cell types is unknown. We present single nuclei RNA sequencing of apical left ventricle tissue from fetal (19-20 weeks, n=3; 27,760 nuclei), non-diseased (ND; 4-14 years, n=3; 16,964 nuclei) and early-onset DCM samples (5-10 years, n=4; 32,712 nuclei) to define the human fetal gene program in dilated cardiomyopathy (DCM), a common cause of heart failure in children and adults."),style="text-align: justify;"),
              h3("Methods and Results:"),
              div(p("We performed single nuclei RNA sequencing with heart tissues from 3 fetal, 3 non-diseased and 4 DCM (for sample characteristics see the Table, and for method see Fig.1A). Single nuclei RNA seq analysis revealed 7 cell clusters across fetal, ND and DCM samples. Also, we investigated that the fetal gene program is broadly re-engaged in cardiomyocytes and cardiac fibroblasts and is not restricted to specific cell sub-populations in DCM. Also, this study provided insights into the central dogma in cardiac biology that the postnatal heart adopts a fetal-like transcriptional state in response to stress. Here, we noted a proportion of genes (<10%) adopting a fetal-like expression pattern in both cardiomyocytes and cardiac fibroblasts, suggesting that multiple cardiac cell types redeploy developmental transcriptional networks in DCM."),style="text-align: justify;"),
              br(),
              div(img(src = "Fig1A.png",height="40%", width="60%"), style="text-align: center;"),
              p("Fig.1A Summary of experimental design showing total number of biological samples and nuclei sequenced, as well as the median number of genes and unique molecular identifiers (UMI) detected per nucleus.",style = "font-size:12px;"),
              tableHTML(fetalinfo,widths = rep(400, 4),caption ="Table.Demographic and Clinical Data of Samples Used in This Study.") %>%
                add_css_row(css = list('background-color', '#f2f2f2'),
                            rows = even(1:4)) %>%
                add_css_row(css = list('background-color', '#e6f0ff'),
                            rows = odd(1:4)) %>%
                add_css_column(css = list('text-align', 'center'), 
                               columns = names(fetalinfo)) %>%
                add_css_header(css = list('text-align', 'center'), headers = 1:4)%>%
                add_css_caption(css = list(c('text-align', 'font-size', 'color'), c('center', '20px', 'darkblue'))),
              br(),
              tableHTML(NDinfo,widths = rep(400, 4)) %>%
                add_css_row(css = list('background-color', '#f2f2f2'),
                            rows = even(1:4)) %>%
                add_css_row(css = list('background-color', '#e6f0ff'),
                            rows = odd(1:4)) %>%
                add_css_column(css = list('text-align', 'center'), 
                               columns = names(NDinfo)) %>%
                add_css_header(css = list('text-align', 'center'), headers = 1:4),
              br(),
              tableHTML(DCMinfo,widths = rep(320, 5)) %>%
                add_css_row(css = list('background-color', '#f2f2f2'),
                            rows = even(1:22)) %>%
                add_css_row(css = list('background-color', '#e6f0ff'),
                            rows = odd(1:22)) %>%
                add_css_column(css = list('text-align', 'center'), 
                               columns = names(DCMinfo)) %>%
                add_css_header(css = list('text-align', 'center'), headers = 1:5),
              h3("Conclusions:"),
              div(p("This work provides insights into the critical gene expression networks that underpin DCM disease pathogenesis in children."),style="text-align: justify;"),
              hr(),
              h4("Please address correspondence to:"),
              HTML(c("Associate Professor Enzo R. Porrello","<br/>","Murdoch Children's Research Institute", "<br/>","Parkville, Melbourne, VIC, 3052, Australia","<br/>","Tel: +61 3 9936 6140","<br/>","Email:","<a>", "enzo.porrello@mcri.edu.au","</a>","<br/>","<br/>","Associate Professor David E. Elliott","<br/>","Murdoch Children's Research Institute","<br/>","Parkville, Melbourne, VIC, 3052, Australia","<br/>","Tel: +61 3 9936 6668","<br/>","Email:","<a>","david.elliott@mcri.edu.au","</a>")),
              br(),
              br(),
              HTML(c("<p>If you tweet about this website, please use the #heartExplorer hashtag. <a href='https://twitter.com/intent/tweet?button_hashtag=heartExplorer&ref_src=twsrc%5Etfw' target = '_blank'><b>Tweet #heartExplorer</b></a>")),
              HTML("<p>For more information about our research projects, please check our <a href = 'https://www.cardioregen.org.au/' target = '_blank'><b>CardioRegen Project Website.</b></a></p>"),
              hr(),
              h4("If you encounter any bugs, please report it to:"),
              HTML(c("Email:","<a>","neda.rahmanimehdiaba@mcri.edu.au","</a>")),
              br(),
              div(h5("with support from"),style="text-align: center;"),
              div(img(src = "logo.png", height="20%", width="40%"), style="text-align: center;")
      ),
      #second tab content
      tabItem(tabName = "visualization",
              tags$head(tags$style(HTML('.content-wrapper {background-color:white;}.selectize-input {border-radius: 90px 90px 90px 90px;border-color: black;text-align:center}.selectize-input.dropdown-active {border-radius:90px;border-color:black;}.selectize-control.single{text-align: center;}'))),
              column(6,align = "center",offset = 3, selectInput("selectgene", "Select a gene:",choices =NULL,  multiple=FALSE, selectize = TRUE)),
              plotOutput("umap") %>% withSpinner(color="black"),
              br(),
              br(),
              br(),
              br(),
              plotOutput("feature"),
              br(),
              column(6,align = "center", offset = 3, selectInput("selectcluster", label = "Select a cell type:", choices = c("Cardiomyocytes"="CM","Endothelial cells" = "Endo","Epicardial cells"="EpC","Fibroblast"="Fib","Immune cells"="Immune","Neurons"="Neuron","Smooth muscle cells"="Smc"), selected = "Cardiomyocytes", multiple=FALSE, selectize=TRUE)),
              fluidRow(
                column(6,plotOutput("vlnclust")),
                column(6,plotOutput("logfcclust")),
              ),
              column(8,offset = 2, plotOutput("vlnclustbybiorep")),
              column(6,align = "center", offset = 3, selectInput("selectgroup", label = "Select a group:", choices = c("Fetal","ND","DCM") , selected = "Fetal", multiple=FALSE, selectize=TRUE)),
              fluidRow(
                column(8, offset = 2, plotOutput("vlngroup"))
              ),
              htmlOutput("text")
             
      ),
      
      # third tab content
      tabItem(tabName = "enrichment", useShinyjs(),
              tags$head(tags$style(HTML('.content-wrapper {background-color:white;}.selectize-input {border-radius: 90px 90px 90px 90px;border-color: black ;text-align:center}.selectize-input.dropdown-active {border-radius:90px;border-color: black;}.selectize-control.single{text-align: center;}.btn{display:block;border-radius: 10%;border: 1px solid black;}'))),
              fluidRow(
                column(5,offset=1, align = "center",radioButtons("group1", "Select a group:", 
                                      choices = c("Fetal","ND","DCM"), inline = TRUE, selected = "Fetal")),
                column(6,align = "center",radioButtons("group2", "Select a group:", 
                                        choices = c("Fetal","ND","DCM"), inline = TRUE, selected = "ND")),
                      ),
              column(6,align = "center", offset = 3, selectInput("selectclusterforenrichment", label = "Select a cell type:", choices = c("Cardiomyocytes"="CM","Endothelial cells" = "Endo","Epicardial cells"="EpC","Fibroblast"="Fib","Immune cells"="Immune","Neurons"="Neuron","Smooth muscle cells"="Smc"), selected = "Cardiomyocytes", multiple=FALSE, selectize=TRUE)),
              fluidRow(
              column(6, offset = 3, align = "center", actionButton("go", "Plot",icon("chart-bar"), 
                                                                   style="color: white; background-color: black; border-color: black; font-size:150%;height: 50px;width: 100px;"))
              ),
              fluidRow(
              column(5, plotOutput("barplot")  %>% withSpinner(color="black")),
              column(7, plotOutput("heatmap")  %>% withSpinner(color="black"))
              ),
              plotOutput("pathway"),
              br(),
              fluidRow(column(6, offset = 3, align = "center", HTML("<b>The list of up-regulated and down-regulated genes (FDR<0.05)</b>"))),
              fluidRow(column(10, offset = 1, DT::dataTableOutput("degtable"))),
              br(),
              fluidRow(column(6, offset = 3, align = "center", HTML("<b>Download the list of up-regulated and down-regulated genes (FDR<0.05)</b>"))),
              fluidRow(column(6, offset =3, align = "center", downloadButton("dl", "Download CSV",icon("download"),style="color: black; background-color: white; border-color: black; font-size:100%; height = 50px;width: 150px;"))),
              
      ),
      #Fourth tab content
      tabItem(tabName = "method",
              h3("Cardiac Nuclei Isolation"),
              div(p("Cardiac nuclei were isolated as previously described[1] with minor modifications. In brief, heart tissue samples were dissected and minced into approximately 1 mm^3 cube pieces and electrically homogenized (IKA) in 15 mL 
                    lysis buffer (0.32 M sucrose, 10 mM Tris-HCl (pH = 8), 5 mM CaCl2, 5 mM magnesium acetate, 2 mM EDTA, 0.5 mM EGTA, 1 mM DTT and 1X Complete Protease Inhibitor). For single nucleus RNA sequencing (RNA-seq), ~20-100 mg of human 
                    left ventricular tissue was used per sample whereas ~0.5-1 g of tissue was required for bulk RNA-seq of purified cardiomyocyte nuclei (mouse and human). The lysate was combined with another 15 mL lysis buffer and subsequently 
                    homogenized with 15-20 strokes using a 40 mL dounce tissue grinder (Wheaton). The cell lysate was then filtered through a 100 uM cell strainer followed by 70 uM and 40 uM cell strainers (BD Falcon) and then centrifuged to pellet 
                    nuclei at 1000xg (Beckman Coulter Allegra X-15R) for 5 mins. Nuclei pellets were then resuspended in 1 M sucrose buffer (1 M sucrose, 10 mM Tris-HCl (pH = 8), 5 mM magnesium acetate, 1 mM DTT and 1X Complete Protease Inhibitor) 
                    and the suspension was cushioned on top of 2X volume of the sucrose buffer, followed by centrifugation to pellet nuclei at 1000xg for 5 mins. All steps were conducted on ice. Isolated cardiac nuclei pellets were washed once in PBS 
                    (Thermo Fisher Scientific) and centrifuged to pellet nuclei at 1000xg. Nuclei pellets were then resuspended in PBS prior to FACS sorting."),style="text-align: justify;"),
              h3("Single Nucleus RNA-seq Library Preparation and Sequencing"),
              div(p("Isolated nuclei were stained with Hoechst 33342 (Thermo Fisher Scientific) prior to sorting on an Influx cell sorter (BD) with 70 uM nozzle and 60 psi pressure setting. Sorted nuclei were counted on a haemocytometer to calculate 
                    nuclei density and then loaded onto the Chromium Controller 4 (10X Genomics) for gel bead emulsion (GEM) formation with ~10,000 nuclei loaded per sample for library preparation. Following GEM formation, library preparation was 
                    conducted according to the manufacturer's recommended protocol using the Chromium Next GEM Single Cell 3' GEM, Library & Gel Bead Kit v3.1. Libraries were sequenced on the NovaSeq 6000 (Illumina) at 50,000 reads per nuclei resolution."),style="text-align: justify;"),
              h3("Bioinformatics Analysis for Single Nuclei RNA Sequencing"),
              div(p("Raw fastq reads for each sample were mapped, processed, and counted using Cell Ranger (v3.0.2). Following this, the counts were then aggregated together to create a table of unique molecular identifier (UMI) counts for 33,939 genes 
                    for each of the samples. All pre-processing and filtering steps of the datasets were subsequently carried out using the R statistical programming language (v3.6.0). The quality of the cells was assessed for each sample independently 
                    by examining the total number of cells, the distributions of total UMI counts, the number of unique genes detected per sample and the proportions of ribosomal and mitochondrial content per cell. Briefly, following the removal of mitochondrial, 
                    ribosomal genes as well as those genes that were not annotated, genes that had at least one count in at least 20 cells were selected for downstream analysis, assuming a minimum cluster size of 20 cells. All genes on the X and Y chromosomes were 
                    removed prior to clustering. For each sample, we performed SCTransform normalization[2], data integration of the biological replicates[3-5], data scaling and graph-based clustering separately, using the R package Seurat (v3.0.2). Data integration 
                    of the biological replicates for each group was performed using CCA[4] from the Seurat package with 30 dimensions and 3000 integration anchors followed by data scaling. Clustering of the cells was performed with 20 principal components (PCs) and 
                    an initial resolution of 0.3. Marker genes to annotate clusters were identified as significantly up-regulated genes for each cluster using moderated t-tests, accounting for the mean variance trend and employing robust empirical Bayes shrinkage of 
                    the variances, followed by TREAT tests specifying a log-fold-change threshold of 0.5 and false discovery rate (FDR) cut off <0.05, using the limma R package (v3.40.2). Moreover, heatmaps showing expression of previously published marker genes were 
                    used to aid in interpretation of the clusters. Visualization of the datasets was primarily carried out using nonlinear dimensionality reduction UMAP[6] plots."),style="text-align: justify;"),
              h3("References"),
              div(p("[1]Bergmann O, Jovinge S. Isolation of cardiomyocyte nuclei from post-mortem tissue.J Vis Exp. 2012; 65:e4205. doi: 10.3791/4205"),style="text-align: justify;"),
              div(p("[2]Hafemeister, C., Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (2019). doi: https://doi.org/10.1186/s13059-019-1874-1"),style="text-align: justify;"),
              div(p("[3]Butler A, Hoffman P, Smibert P, Papalexi E, Satija R. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat Biotechnol. 2018 Jun;36(5):411-420. doi: 10.1038/nbt.4096"),style="text-align: justify;"),
              div(p("[4]Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM 3rd, Hao Y, Stoeckius M, Smibert P, Satija R. Comprehensive Integration of Single-Cell Data. Cell. 2019 Jun 13;177(7):1888-1902.e21. doi: 10.1016/j.cell.2019.05.031"),style="text-align: justify;"),
              div(p("[5]Stuart, T., Satija, R. Integrative single-cell analysis. Nat Rev Genet 20, 257-272 (2019). https://doi.org/10.1038/s41576-019-0093-7"),style="text-align: justify;"),
              div(p("[6]Becht E, McInnes L, Healy J, Dutertre CA, Kwok IWH, Ng LG, Ginhoux F, Newell EW. Dimensionality reduction for visualizing single-cell data using UMAP. Nat Biotechnol. 2018 Dec 3. doi: 10.1038/nbt.4314"),style="text-align: justify;"),
      ),
      #Fifth tab content
      tabItem(tabName = "code",
              h3("Data and code availability"),
              div(p("All snRNA-seq raw fastq.gz files including sample details have been deposited to Gene Expression Omnibus under accession No. GSE185100."),style="text-align: justify;"),
              HTML(c("<p>Comprehensive bioinformatics analyses with the source code can be retrieved from <a href='https://neda-mehdiabadi.github.io/Fetal-Gene-Program-snRNAseq/' target = '_blank'><b> analysis website</b></a> and <a href='https://github.com/neda-mehdiabadi/Fetal-Gene-Program-snRNAseq' target = '_blank'><b> source code</b></a>, respectively.")),
              HTML(c("<p>Also, the code used to build this website is available on <a href='https://neda-mehdiabadi.github.io/Fetal-Gene-Program-snRNAseq/' target = '_blank'><b> GitHub/HeartExplorer</b></a>.")),
      )
    )
  )
)




server <- function(input,output,session) {
  
  ds <- c()
  csvFile <- data.frame()

  updateSelectizeInput(session, "selectgene", choices = rownames(cds),selected=sample(c("ANLN","TNNT2","ACTN2","NRXN1","NRXN3","DCN","VWF","PECAM1","MYH11","ACTA2","PTPRC","PDGFRB"),1), server = T)
  #tab "visualization", umap plot
  output$umap <- renderPlot({
    DimPlot(cds, reduction = "umap",label=FALSE,label.size = 6, cols = c("Er" = "brown4", "CM(Prlf)"= "darkorange2","CM"="#7570B3","Endo"="#E7298A","EpC"="#66A61E","Fib"="#E6AB02","Immune"="#A6761D","Neuron"="#666666","Smc"="#1F78B4"), split.by = "orig.ident")
  })
  
  #tab "visualization", feature plot
  output$feature <- renderPlot({
    gene= input$selectgene
    FeaturePlot(cds, features = gene, cols = c("lightgrey","blue4"), pt.size = 0.5,ncol = 3, split.by = "orig.ident", reduction = "umap")
  })
  
  #tab "visualization", violin plot for a gene split.by=orig.ident
  output$vlnclust <- renderPlot({
    Idents(cds) <- cds$BCT_origIdent
    gene <- input$selectgene
    cluster <- input$selectcluster
    chr <- tolower(substr(cluster,1,2))
    col <- which(substr(colnames(comparison),1,8) %in% c(paste("fdr.DF",chr,sep=""),paste("fdr.DN",chr,sep=""),paste("fdr.FN",chr,sep="")))
    max = max(cds@assays$SCT@data[gene,])
    p <- VlnPlot(cds, features = gene,pt.size = 0.1,  split.by = "orig.ident",idents = paste(levels(cds$orig.ident),cluster,sep = "."),y.max = max+2, split.plot = F)+xlab("")+NoLegend()
    index <- which(comparison$SYMBOL==gene)
    if(comparison[index,col[1]] >= 0.05){
      DF.adj.p <- "ns"
    }else{
      DF.adj.p <- signif(comparison[index,col[1]],digits=3)
    }
    if(comparison[index,col[2]] >= 0.05){
      DY.adj.p <- "ns"
    }else{
      DY.adj.p <- signif(comparison[index,col[2]],digits=3)
    }
    if(comparison[index,col[3]] >= 0.05){
      FY.adj.p <- "ns"
    }else{
      FY.adj.p <- signif(comparison[index,col[3]],digits=3)
    }
    p + annotate(geom = "segment", x = c(1,2,1), y = c(max+1,max+0.5,max+1.5), xend=c(2,3,3) ,yend =c(max+1,max+0.5,max+1.5))+
      annotate(geom = "text", x=c(2.5,1.5,2), y=c(max+0.7,max+1.2,max+1.7), label=c(DY.adj.p,FY.adj.p,DF.adj.p))+
      ggtitle(paste("Expression Level of ",gene," in ",cluster,sep = ""))+theme(plot.title = element_text(hjust = 0.5,vjust=-2, size = 15))
  })
  
  
  #tab "visualization", logFC plot for a gene split.by = orig.ident
  output$logfcclust <- renderPlot({
    Idents(cds) <- cds$Broad_celltype
    gene <- input$selectgene
    cluster <- input$selectcluster
    chr <- tolower(substr(cluster,1,2))
    col <- which(substr(colnames(comparison),1,10) %in% c(paste("LogFC.DF",chr,sep=""),paste("LogFC.DN",chr,sep=""),paste("LogFC.FN",chr,sep="")))
    subdata <- data.frame(matrix(NA, nrow = 3,ncol = 3))
    colnames(subdata) <- c("log2FC","Comparison","adjPvalue") 
    subdata$Comparison <- c("DCM vs Fetal","DCM vs ND","Fetal vs ND")
    index <- which(comparison$SYMBOL==gene)
    subdata$log2FC[1] <- comparison[index,col[1]]
    subdata$log2FC[2] <- comparison[index,col[2]]
    subdata$log2FC[3] <- comparison[index,col[3]]
    col <- which(substr(colnames(comparison),1,8) %in% c(paste("fdr.DF",chr,sep=""),paste("fdr.DN",chr,sep=""),paste("fdr.FN",chr,sep="")))
    for (i in 1:3) {
      if (comparison[index,col[i]] >= 0.05) {
        subdata$adjPvalue[i] <- "ns"
      } else {
        subdata$adjPvalue[i] <- signif(comparison[index,col[i]],digits=3)
      }
    }
    subdata <- as_tibble(subdata)
    ##abs(log2FC) < 0.006 will be visualized as well 
    index <- which(abs(subdata[,1]) < 0.006)
    subdata[index,1] <- sign(subdata[index,1])*0.006
    subdata <- mutate(subdata, condition=if_else(log2FC>0, "pos", "neg"))
    ggplot(subdata, aes(x=Comparison, y=log2FC, fill=condition)) +
      geom_bar(stat = "identity", width = 0.2) +
      geom_text(aes(label=adjPvalue), vjust=-0.5, size = 4)+
      scale_fill_manual(values = c("pos"="blue","neg"="red"))+
      theme_light(base_size=15) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=10),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10)
      )+
      xlab("") +
      ylab("Log2 Fold Change")+
      ylim(min(0,min(subdata$log2FC))-1,max(0,max(subdata$log2FC)+1))+
      ggtitle(paste(gene," in ",cluster,sep = ""))+
      theme(plot.title = element_text(hjust = 0.5,vjust=-3, size = 15, face = "bold"))
  })
  
  #tab "visualization", violin plot for a gene split.by=biorep
  output$vlnclustbybiorep <- renderPlot({
    Idents(cds) <- cds$BCT_biorep
    gene <- input$selectgene
    cluster <- input$selectcluster
    chr <- tolower(substr(cluster,1,2))
    max = max(cds@assays$SCT@data[gene,])
    p <- VlnPlot(cds, features = gene,pt.size = 0.1, idents = paste(levels(cds$biorep),cluster,sep = "."), y.max = max+2, split.plot = F)+xlab("")+NoLegend()
    p + ggtitle(paste("Expression Level of ",gene," in ",cluster,sep = ""))+theme(plot.title = element_text(hjust = 0.50,vjust=-2, size = 15))
  })
  
  #tab "visualization", violin plot for a gene split.by=Broad_celltype
  output$vlngroup <- renderPlot({
    Idents(cds) <- cds$BCT_origIdent
    gene <- input$selectgene
    group <- input$selectgroup
    if (group == "Fetal") {
      Cols = c("brown4","darkorange2","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#1F78B4")
    } else {
      Cols = c("#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#1F78B4")
    }
    p <- VlnPlot(cds, features = gene,idents =paste(group,levels(cds$Broad_celltype),sep = "."), split.plot = F, cols = Cols)+xlab("")+NoLegend()
    p + ggtitle(paste("Expression Level of ",gene," in ",group,sep = ""))+theme(plot.title = element_text(hjust = 0.47,vjust=-2, size = 15))
    })
  
  #tab "visualization", text indicating adj. P value for a gene split.by=Broad_celltype
  output$text <- renderUI({
    group <- input$selectgroup
    if(group=="Fetal"){
      ds <- read.fst("data/Fetal.fst")
    }else if (group=="ND"){
      ds <- read.fst("data/ND.fst")
    }else{
      ds <- read.fst("data/DCM.fst")
    }
  index <- which(ds$gene==input$selectgene)
   if (length(index) != 0){
       str <- c("According to the violin plot above:",'<br/>')
       HTML(c(str,paste(ds$gene[index]," is differentially expressed in ",ds$cluster[index], " with adjusted p-value:", signif(ds$p_val_adj[index],digits = 3),'<br/>')))
    } else {
      HTML(c("According to the violin plot above:","<br/>","The gene is not differentially expressed in any cell-types!"))
    }
  })
  
  #tab "pathway enrichment analysis"
  toListen <- reactive({
    list(input$group1,input$group2)
  })
  observeEvent(toListen(), {
    if(identical(input$group1, input$group2)) {
      disable("go")
    } else  { 
      enable("go")
    }
  })
  #tab "pathway enrichment analysis", bar plot for No. DEG
  clicked <- eventReactive(input$go, {
    cluster <- input$selectclusterforenrichment
    chr <- tolower(substr(cluster,1,2))
    colfdr <- which(substr(colnames(comparison),1,8) == c(paste("fdr.",label$short[label$long==input$group1],label$short[label$long==input$group2],chr,sep="")))
    collog <- which(substr(colnames(comparison),1,10) == c(paste("LogFC.",label$short[label$long==input$group1],label$short[label$long==input$group2],chr,sep="")))
    colLogCPM <- which(substr(colnames(sumexpr),1,2)==chr & substr(sub(".*\\.", "", colnames(sumexpr)),1,nchar(sub(".*\\.", "", colnames(sumexpr)))-1) %in% c(input$group1,input$group2))
    dt <- data.frame(matrix(NA,nrow = 2, ncol = 2))
    colnames(dt) <- c("Var1","Freq")
    dt$Var1[1:2] <- c("Up","Down")
    dt$Freq[1] <- dim(comparison %>% filter(comparison[,colfdr] < 0.05 & comparison[,collog]>0))[1]
    dt$Freq[2] <- dim(comparison %>% filter(comparison[,colfdr] < 0.05 & comparison[,collog] < 0))[1]
    barplot(dt$Freq,beside=TRUE,legend=TRUE,col=c("red","darkblue"),ylab="Number of significant genes (FDR<0.05)", ylim = c(0,(max(dt$Freq)+500)),main=paste(dt$Freq[1],"genes up-regulated in",input$group1,"\nand",dt$Freq[2],"genes up-regulated in",input$group2), cex.main = 1,cex.lab = 1)
    legend("topright", 
           legend = c("Up       (log2FC>0)", "Down (log2FC<0)"), 
           fill = c("red", "darkblue"),y=(max(dt$Freq)+100))
    up <- comparison %>% filter(comparison[,colfdr] < 0.05 & comparison[,collog]>0)
    up <- up[order(up[,collog],decreasing = T),]
    dn <- comparison %>% filter(comparison[,colfdr] < 0.05 & comparison[,collog] < 0)
    dn <- dn[order(dn[,collog],decreasing = F),]
    index.up <- match(up$SYMBOL,sumexpr$symbol)
    index.dn <- match(dn$SYMBOL,sumexpr$symbol)
    csvFile <- cbind(up[,c(1:5,collog,colfdr)],sumexpr[index.up,colLogCPM])
    csvFile <<- rbind(csvFile,cbind(dn[,c(1:5,collog,colfdr)],sumexpr[index.dn,colLogCPM]))
    colnames(csvFile)[6:7] <<- c(paste("log2FC(",input$group1,"/",input$group2,") in ",cluster, sep = ""),"FDR")
    colnames(csvFile)[8:dim(csvFile)[2]] <<- paste(colnames(csvFile)[8:dim(csvFile)[2]],"(log2CPM)", sep = "\t")
  }, ignoreNULL = FALSE)
  
  #tab "pathway enrichment analysis", pathway analysis
  clickedPath <- eventReactive(input$go, {
    cluster <- input$selectclusterforenrichment
    chr <- tolower(substr(cluster,1,2))
    coltstat <- which(substr(colnames(comparison),1,10) == c(paste("tstat.",label$short[label$long==input$group1],label$short[label$long==input$group2],chr,sep="")))
    c2.id <- ids2indices(Hs.c2,comparison$ENTREZID)
    reactome.id <-c2.id[grep("REACTOME",names(c2.id))]
    camera <- cameraPR(comparison[,coltstat],reactome.id)
    camera.dn <- camera[camera[,2]=="Down",]
    camera.up <- camera[camera[,2]=="Up",]
    camera.up <- camera.up %>% filter(camera.up$FDR < 0.05)
    camera.up <- camera.up[order(camera.up$NGenes,decreasing = T),]
    camera.dn <- camera.dn %>% filter(camera.dn$FDR < 0.05)
    camera.dn <- camera.dn[order(camera.dn$NGenes,decreasing = T),]
    nsets.up <- min(10,dim(camera.up)[1])
    nsets.dn <- min(10,dim(camera.dn)[1])
    if(nsets.dn !=0 & nsets.up!=0){
      all.cam <- rbind(camera.up[1:nsets.up,], camera.dn[1:nsets.dn,])
      scores <- -log10(all.cam$FDR)
      names(scores) <- rownames(all.cam)
      names(scores) <- gsub("REACTOME_","",names(scores))
      names(scores) <- substr(names(scores),1,250)
      par(mfrow=c(1,1))
      par(mar=c(5,60,3,2))
      barplot(scores[length(scores):1],horiz = T,las=1,col=c(rep("darkblue",nsets.dn),rep("red",nsets.up)),cex.names=0.9,
              cex.axis = 0.8,xlab="-log10(FDR)",cex.lab=0.9, xlim = c(0,max(scores)+2))
      abline(v= -log10(0.05),lty=2, col= "black")
      title(main = paste("Top", nsets.up,"up-regulated pathways in",input$group1,"\nFollowed by top",nsets.dn,"up-regulated pathways in",input$group2,"\n -log10(0.05):black dashed line"),adj = 0.5, line = 0.3, cex.main = 1)
    }else if(nsets.dn ==0 & nsets.up !=0 ){
      all.cam <- camera.up[1:nsets.up,]
      scores <- -log10(all.cam$FDR)
      names(scores) <- rownames(all.cam)
      names(scores) <- gsub("REACTOME_","",names(scores))
      names(scores) <- substr(names(scores),1,250)
      par(mfrow=c(1,1))
      par(mar=c(5,60,3,2))
      barplot(scores[length(scores):1],horiz = T,las=1,col=rep("red",each=nsets.up),cex.names=0.9,
              cex.axis = 0.8,xlab="-log10(FDR)",cex.lab=0.9, xlim = c(0,max(scores)+2))
      abline(v= -log10(0.05),lty=2, col= "black")
      title(main = paste("Top", nsets.up,"up-regulated pathways in",input$group1,"\nFollowed by 0 up-regulated pathways in",input$group2,"\n -log10(0.05):black dashed line"),adj = 0.5, line = 0.3, cex.main = 1)
    } else if(nsets.dn !=0 & nsets.up ==0 ){
      all.cam <- camera.dn[1:nsets.dn,]
      scores <- -log10(all.cam$FDR)
      names(scores) <- rownames(all.cam)
      names(scores) <- gsub("REACTOME_","",names(scores))
      names(scores) <- substr(names(scores),1,250)
      par(mfrow=c(1,1))
      par(mar=c(5,60,3,2))
      barplot(scores[length(scores):1],horiz = T,las=1,col=rep("darkblue",each=nsets.dn),cex.names=0.9,
              cex.axis = 0.8,xlab="-log10(FDR)",cex.lab=0.9, xlim = c(0,max(scores)+2))
      abline(v= -log10(0.05),lty=2, col= "black")
      title(main = paste("0 up-regulated pathways in",input$group1,"\nFollowed by top",nsets.dn,"up-regulated pathways in",input$group2,"\n -log10(0.05):black dashed line"),adj = 0.5, line = 0.3, cex.main = 1)
    }
  }, ignoreNULL = FALSE)
  
  #tab "pathway enrichment analysis", heatmap for top genes
  clickedheatmap <- eventReactive(input$go, {
    cluster <- input$selectclusterforenrichment
    chr <- tolower(substr(cluster,1,2))
    colfdr <- which(substr(colnames(comparison),1,8) == c(paste("fdr.",label$short[label$long==input$group1],label$short[label$long==input$group2],chr,sep="")))
    collog <- which(substr(colnames(comparison),1,10) == c(paste("LogFC.",label$short[label$long==input$group1],label$short[label$long==input$group2],chr,sep="")))
    up <- comparison %>% filter(comparison[,colfdr] < 0.05 & comparison[,collog]>0)
    up <- up[order(up[,collog], decreasing = T),]
    dn <- comparison %>% filter(comparison[,colfdr] < 0.05 & comparison[,collog]<0)
    dn <- dn[order(dn[,collog], decreasing = F),]
    top10.up <- up$SYMBOL[1:10]
    top10.up <- na.omit(top10.up)
    top10.dn <- dn$SYMBOL[1:10]
    top10.dn <- na.omit(top10.dn)
    col <- which(substr(colnames(sumexpr),1,2)==chr & substr(sub(".*\\.", "", colnames(sumexpr)),1,nchar(sub(".*\\.", "", colnames(sumexpr)))-1) %in% c(input$group1,input$group2))
    row.up <- which(sumexpr$symbol %in% top10.up)
    row.dn <- which(sumexpr$symbol %in% top10.dn)
    row <- c(row.up,row.dn)
    aheatmap(sumexpr[row,col],Rowv = NA,Colv = NA, labRow = sumexpr$symbol[row],
             fontsize=10,color="-RdYlBu",cexRow =1.2, cexCol = 1.2,labCol = colnames(sumexpr)[col],
             scale="row", main =paste("Top",length(top10.up),"genes up-regulated (FDR<0.05) in",input$group1,"\nFollowed by top",length(top10.dn),"genes up-regulated (FDR<0.05) in",input$group2,"\n**Row-Z Score values**"))
  
  }, ignoreNULL = FALSE)
    
  clickedtable <- eventReactive(input$go, {
    csvFile[,6:dim(csvFile)[2]] <-  signif(csvFile[,6:dim(csvFile)[2]],digits=3)
    csvFile %>% DT::datatable(rownames = FALSE,options = list(
      lengthMenu = c(10, 25, 50, 100), pageLength = 10, scrollX = TRUE,
      columnDefs = list(list(className = "dt-center", targets = "_all")),
      paging = TRUE, selection='none',initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'font-size': '90%','background-color': '#000', 'color': '#fff'});",
        "}")
    ) 
    ) %>% DT::formatStyle(columns = c(0:dim(csvFile)[2]), fontSize = '90%')
  }, ignoreNULL = FALSE)
  
  #tab "pathway enrichment analysis", bar plot for No. DEG
  output$barplot <- renderPlot({
    clicked()
  })
  
  #tab "pathway enrichment analysis", pathway analysis
  output$pathway <- renderPlot({
    clickedPath()
  })
  
  #tab "pathway enrichment analysis", heatmap for top genes
  output$heatmap <- renderPlot({
    clickedheatmap()
    
  })
  output$dl <- downloadHandler(
    filename = "gene-list.csv",
    content = function(file) {
    write.csv(csvFile, file, row.names = FALSE)
    }
  )
  #tab "pathway enrichment analysis", data table
  output$degtable <- DT::renderDataTable({
   clickedtable()
  })
  
}

shinyApp(ui, server)