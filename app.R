library(shiny)
library(tidyverse)
library(protoclust)
library(plotly)
library(DT)
library(readxl)
source("scripts/helpers.R")

# It is important that this piece of code go first, so that R understands that
# this is a shiny app
ui <- fluidPage(title = "Choosing gene expression",
                tabsetPanel(
                  tabPanel(title = "CDS",
                           sidebarLayout(
                             sidebarPanel(
                               selectInput(inputId="cluster_this_cds",
                                           label="What should be clustered?",
                                           choices=list("Log Fold Changes (with respect to Control)" = "logfc",
                                                        "Positive normalized counts (VST)" = "poscounts")),
                               selectInput(inputId = "c_method_cds",
                                           label = "Clustering Method:",
                                           choices = list("ward.D2", "minmax", "complete", "median", "centroid")),
                               selectInput(inputId = "d_metric_cds",
                                           label = "Distance Metric:",
                                           choices = list("cosine", "spearman", "pearson")),
                               numericInput(inputId = "k_clusters_cds",
                                            label = "Number of clusters:",
                                            value = 20, min=1),
                               actionButton(inputId = "cluster_param_cds", 
                                            label = "Submit", icon = NULL, width = NULL),
                               numericInput(inputId = "high_cds",
                                            label="Highlight Cluster:",
                                            value=0, min = 0),
                               actionButton(inputId = "high_param_cds",
                                            label = "Highlight")),
                             mainPanel(plotOutput("dendo_cds"),
                                       plotOutput("wss_cds"))
                           ),
                           sidebarLayout(
                             sidebarPanel(selectInput(inputId = "color_by_cds",
                                           label = "Color genes according to:",
                                           choices = list("organism", 
                                                          "prototype (if minmax is used)",
                                                          "Salmonella pre-selection",
                                                          "Campylobacter pre-selection"))
                             ),
                             mainPanel(plotOutput("cluster_wrap_cds"),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br(),
                                       br()),
                           ),
                           sidebarLayout(
                             sidebarPanel(numericInput(inputId = "cluster_interest_cds",
                                                       label = "Focus on Cluster:",
                                                       value=1, min=1),
                                          selectInput(inputId = "organism_interest_cds",
                                                      label = "Show:",
                                                      choices = c("Both organisms", "Only Salmonella", "Only Campylobacter")),
                                          textInput(inputId = "gene_interest_cds",
                                                    label="Highlight Gene (if more than one, separate with comma)",
                                                    value = "None"),
                                          actionButton(inputId="focus_param_cds",
                                                       label = "Submit")),
                             
                             mainPanel(plotlyOutput("focus_profile_cds"),
                                       DTOutput("cluster_members_cds"))
                           ),
                           "All CDS",
                           DTOutput("all_cds_table")
                           
                      ),
                  
                  
                  tabPanel(title = "sRNA",
                           sidebarLayout(
                             sidebarPanel(
                               selectInput(inputId="cluster_this_srna",
                                           label="What should be clustered?",
                                           choices=list("Log Fold Changes (with respect to Control)"="logfc",
                                                        "Positive normalized counts (VST)" = "poscounts")),
                               selectInput(inputId = "c_method_srna",
                                           label = "Clustering Method:",
                                           choices = list("ward.D2", "minmax", "complete", "median", "centroid")),
                               selectInput(inputId = "d_metric_srna",
                                           label = "Distance Metric:",
                                           choices = list("cosine", "spearman","pearson")),
                               numericInput(inputId = "k_clusters_srna",
                                            label = "Number of clusters:",
                                            value = 20, min=1),
                               actionButton(inputId = "cluster_param_srna", 
                                            label = "Submit", icon = NULL, width = NULL),
                               numericInput(inputId = "high_srna",
                                            label="Highlight Cluster:",
                                            value=0, min = 0),
                               actionButton(inputId = "high_param_srna",
                                            label = "Highlight")),
                             mainPanel(plotOutput("dendo_srna"),
                                       plotOutput("wss_srna"))
                           ),
                           sidebarLayout(
                             sidebarPanel(selectInput(inputId = "color_by_srna",
                                                                  label = "Color genes according to:",
                                                                  choices = list("organism", 
                                                                                 "prototype (if minmax is used)",
                                                                                 "Salmonella pre-selection"))
                           ),
                           mainPanel(plotOutput("cluster_wrap_srna"))
                           ),
                           sidebarLayout(
                             sidebarPanel(numericInput(inputId = "cluster_interest_srna",
                                                       label = "Focus on Cluster:",
                                                       value=1, min=1),
                                          selectInput(inputId = "organism_interest_srna",
                                                      label = "Show:",
                                                      choices = c("Both organisms", "Only Salmonella", "Only Campylobacter")),
                                          textInput(inputId = "gene_interest_srna",
                                                    label="Highlight Gene (if more than one, separate with comma)",
                                                    value = "None"),
                                          actionButton(inputId="focus_param_srna",
                                                       label = "Submit")),
                             
                             mainPanel(plotlyOutput("focus_profile_srna"),
                                       DTOutput("cluster_members_srna"))
                           ),
                           "All sRNAs",
                           DTOutput("all_srnas_table")
                           ) # This tab
                ) # Tabset
              ) #fluidui

# Here I will import the data to be used for the rest of the analysis. 
# Since we only really need to import this once, we can leave it outside of
# the server function. The data I call here was performed with clustering_usw

cds_log.fc <- read_tsv("shiny_data/cds_logFC.tsv") %>% 
  column_to_rownames("gene_name")

srna_log.fc <- read_tsv("shiny_data/srna_logFC.tsv") %>% 
  column_to_rownames("gene_name")

# Some average expression info
joint_cds.avg <- read_tsv("shiny_data/joint_cds_avg.tsv") %>% 
  column_to_rownames("gene_name")

cds_rowmeans <- rowMeans(joint_cds.avg)
cds_avgExp <- data.frame("gene_name"=names(cds_rowmeans),
                         "Average_Expression"=cds_rowmeans)


joint_srna.avg <- read_tsv("shiny_data/joint_srna_avg.tsv") %>% 
  column_to_rownames("gene_name")

srna_rowmeans <- rowMeans(joint_srna.avg)
srna_avgExp <- data.frame("gene_name"=names(srna_rowmeans),
                         "Average_Expression"=srna_rowmeans)

# Genomic features
salm_features <- read_tsv("shiny_data/merged_features_salmonella.tsv") %>% 
  mutate(organism="Salmonella",
         gene_name=paste0(locus_tag, "-Salmonella")) %>% 
  select(gene_name, locus_tag, symbol, name)

camp_features <- read_tsv("shiny_data/merged_features_campylobacter.tsv") %>% 
  mutate(organism="Campylobacter",
         gene_name=paste0(locus_tag, "-Campylobacter")) %>% 
  select(gene_name, locus_tag, symbol, name)

genomic_features <- bind_rows(salm_features, camp_features)


# Some pre-selected features for genes
salm.preselection <- read_tsv("shiny_data/pre_selected_salmonella.txt") %>% 
  rename("symbol"=ID) %>% 
  left_join(salm_features, by="symbol") %>% 
  filter(!is.na(locus_tag)) %>% 
  mutate(gene_name = paste0(locus_tag, "-Salmonella"))%>% 
  select(gene_name, symbol)

camp.preselection <- read_excel("shiny_data/Campy regulators and stress resp for Roberto.xlsx") %>% 
  rename("locus_tag" = `Locus tag 81-176`) %>% 
  left_join(camp_features, by="locus_tag") %>% 
  mutate(gene_name = paste0(locus_tag, "-Campylobacter")) %>% 
  select(gene_name, gene, `clone?`)



# Here we start our app

server <- function(input, output){
  
  # Default clustering parameters
  ## CDS
  clustering_parameters_cds <- reactiveValues(k_clusters = 20,
                                          c_method = "ward.D2",
                                          d_metric = "cosine",
                                          cluster_this = "logfc")
  ## sRNA
  clustering_parameters_srna <- reactiveValues(k_clusters = 20,
                                              c_method = "ward.D2",
                                              d_metric = "cosine",
                                              cluster_this = "logfc")
  
  # Initialize parameters
  ## CDS
  observeEvent(input$cluster_param_cds, {
    clustering_parameters_cds$k_clusters <- input$k_clusters_cds
    clustering_parameters_cds$d_metric <- input$d_metric_cds
    clustering_parameters_cds$c_method <- input$c_method_cds
    clustering_parameters_cds$cluster_this <- input$cluster_this_cds
  })
  
  ## sRNA
  observeEvent(input$cluster_param_srna, {
    clustering_parameters_srna$k_clusters <- input$k_clusters_srna
    clustering_parameters_srna$d_metric <- input$d_metric_srna
    clustering_parameters_srna$c_method <- input$c_method_srna
    clustering_parameters_srna$cluster_this <- input$cluster_this_srna
  })
  
  renderPrint(clustering_parameters_cds)
  
  # Clustering  
  ## CDS
  clust_cds <- reactive({
    cluster_counts(what_to_cluster=clustering_parameters_cds$cluster_this, 
                   distance_metric=clustering_parameters_cds$d_metric,
                   how_to_cluster=clustering_parameters_cds$c_method,
                   logfc_df=cds_log.fc,
                   posc_df=joint_cds.avg)
  })

  ## sRNA
  clust_srna <- reactive({
    cluster_counts(what_to_cluster=clustering_parameters_srna$cluster_this, 
                                    distance_metric=clustering_parameters_srna$d_metric, 
                                    how_to_cluster=clustering_parameters_srna$c_method,
                                    logfc_df=srna_log.fc, 
                                    posc_df=joint_srna.avg)
  })
  
  # Dendrogram Parameters
  # Here would go the clusters that we want to hightlight
  ## CDS
  dendro_high_cds <- reactiveValues(
    c_high = 0
    )
  
  ##sRNA
  dendro_high_srna <- reactiveValues(
    c_high = 0
  )
  
  # Update highlight
  ## CDS
  observeEvent(input$high_param_cds, {
    dendro_high_cds$c_high <- input$high_cds
  })
  
  ## sRNA
  observeEvent(input$high_param_srna, {
    dendro_high_srna$c_high <- input$high_srna
  })
  
  # Dendogram
  ## CDS
  output$dendo_cds <- renderPlot({
    
    plot_dendrogram(how_to_cluster=clustering_parameters_cds$c_method, 
                    number_of_clusters=clustering_parameters_cds$k_clusters, 
                    highlight_this_cluster=dendro_high_cds$c_high, 
                    clustered_object=clust_cds()$cluster_object)
    
  })
  
  ## sRNA
  output$dendo_srna <- renderPlot({
    
    plot_dendrogram(how_to_cluster=clustering_parameters_srna$c_method, 
                    number_of_clusters=clustering_parameters_srna$k_clusters, 
                    highlight_this_cluster=dendro_high_srna$c_high, 
                    clustered_object=clust_srna()$cluster_object)
  })
  
  # Within Sum Squares plot
  output$wss_cds <- renderPlot({
    plot_wss(clustered_object=clust_cds()$cluster_object, 
             what_was_clustered=clustering_parameters_cds$cluster_this, 
             how_was_it_clustered=clustering_parameters_cds$c_method, 
             number_of_clusters=clustering_parameters_cds$k_clusters, 
             logfc_df = cds_log.fc, 
             posc_df = joint_cds.avg)
  })
  
  output$wss_srna <- renderPlot({
    plot_wss(clustered_object=clust_srna()$cluster_object, 
             what_was_clustered=clustering_parameters_srna$cluster_this, 
             how_was_it_clustered=clustering_parameters_srna$c_method, 
             number_of_clusters=clustering_parameters_srna$k_clusters, 
             logfc_df = srna_log.fc, 
             posc_df = joint_srna.avg)
  })
  
  # Cluster Data Frame
  ## CDS
  cluster_info_cds <- reactive({
    cluster_longform(how_to_cluster=clustering_parameters_cds$c_method, 
                     number_of_clusters=clustering_parameters_cds$k_clusters, 
                     clustered_object=clust_cds()$cluster_object,
                     logfc_dataframe=cds_log.fc)
    
  })
  
 
  ## sRNA
  cluster_info_srna <- reactive({
    cluster_longform(how_to_cluster=clustering_parameters_srna$c_method, 
                     number_of_clusters=clustering_parameters_srna$k_clusters, 
                     clustered_object=clust_srna()$cluster_object,
                     logfc_dataframe=srna_log.fc)
    
  })
  
  # Cluster Profiles
  ## CDS
  output$cluster_wrap_cds <- renderPlot({
    cluster_profile(dataframe_for_plot = cluster_info_cds(), 
                    plot_according_to = input$color_by_cds,
                    salmonella_preselection = salm.preselection, 
                    campylobacter_preselection = camp.preselection)
  }, height = 750)
  
  ## sRNA
  output$cluster_wrap_srna <- renderPlot({
    cluster_profile(dataframe_for_plot = cluster_info_srna(), 
                    plot_according_to = input$color_by_srna,
                    salmonella_preselection = salm.preselection, 
                    campylobacter_preselection = camp.preselection)
  })
  
  # Focused Data Parameters
  ## CDS
  focus_parameters_cds <- reactiveValues(c_focus = 1,
                                         g_focus = c("None"),
                                         o_focus = "Both organisms")
  
  ## sRNA
  focus_parameters_srna <- reactiveValues(c_focus = 1,
                                          g_focus = c("None"),
                                          o_focus = "Both organisms")
  
  # Update Focused Parameters
  ## CDS
  observeEvent(input$focus_param_cds, {
    focus_parameters_cds$c_focus <- input$cluster_interest_cds
    focus_parameters_cds$g_focus <- strsplit(input$gene_interest_cds, ",")[[1]]
    focus_parameters_cds$o_focus <- input$organism_interest_cds
  })
  
  ## sRNA
  observeEvent(input$focus_param_srna, {
    focus_parameters_srna$c_focus <- input$cluster_interest_srna
    focus_parameters_srna$g_focus <- strsplit(input$gene_interest_srna, ",")[[1]]
    focus_parameters_srna$o_focus <- input$organism_interest_srna
  })
  
  
  # Focused Plot
  ## CDS
  output$focus_profile_cds <- renderPlotly({
    
    focused_profile(which_organism = focus_parameters_cds$o_focus, 
                    which_cluster = focus_parameters_cds$c_focus, 
                    which_gene = focus_parameters_cds$g_focus, 
                    dataframe_for_plot = cluster_info_cds())
  })
  
  ## sRNA
  output$focus_profile_srna <- renderPlotly({
    
    focused_profile(which_organism = focus_parameters_srna$o_focus, 
                    which_cluster = focus_parameters_srna$c_focus, 
                    which_gene = focus_parameters_srna$g_focus, 
                    dataframe_for_plot = cluster_info_srna())
  })
  
  
  # Focused Data
  ## CDS
  output$cluster_members_cds <- renderDT({
    
    cluster_table(dataframe_for_info = cluster_info_cds(), 
                  which_cluster = focus_parameters_cds$c_focus, 
                  distance_matrix = clust_cds()$distance_matrix,
                  genomic_features_dataframe = genomic_features,
                  salmonella_preselection = salm.preselection,
                  campylobacter_preselection = camp.preselection,
                  average_expression_information = cds_avgExp,
                  how_to_cluster = clustering_parameters_cds$c_method)
    
  })
  
  ## sRNA
  output$cluster_members_srna <- renderDT({
    cluster_table(dataframe_for_info = cluster_info_srna(), 
                  which_cluster = focus_parameters_srna$c_focus, 
                  distance_matrix = clust_srna()$distance_matrix,
                  genomic_features_dataframe = genomic_features,
                  salmonella_preselection = salm.preselection,
                  campylobacter_preselection = camp.preselection,
                  average_expression_information = srna_avgExp,
                  how_to_cluster = clustering_parameters_srna$c_method)
  })
  
  output$all_cds_table <- renderDT({
    
    complete_table(dataframe_with_info = cluster_info_cds(), 
                   genomic_features_dataframe = genomic_features)
    
  })
  
  output$all_srnas_table <- renderDT({
    
    complete_table(dataframe_with_info = cluster_info_srna(), 
                   genomic_features_dataframe = genomic_features)
    
  })
  
}

shinyApp(ui = ui, server = server)
