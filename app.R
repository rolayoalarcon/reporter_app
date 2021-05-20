library(shiny)
library(tidyverse)
library(protoclust)
library(plotly)
library(DT)

# It is important that this piece of code go first, so that R understands that
# this is a shiny app
ui <- fluidPage(title = "Choosing gene expression",
                tabsetPanel(
                  tabPanel(title = "CDS",
                           sidebarLayout(
                             sidebarPanel(
                               selectInput(inputId = "c_method_cds",
                                           label = "Clustering Method:",
                                           choices = list("minmax", "complete")),
                               selectInput(inputId = "d_metric_cds",
                                           label = "Distance Metric:",
                                           choices = list("cosine", "pearson")),
                               numericInput(inputId = "k_clusters_cds",
                                            label = "Number of clusters:",
                                            value = 20, min=1),
                               actionButton(inputId = "cluster_param_cds", 
                                            label = "Submit", icon = NULL, width = NULL)),
                             mainPanel(plotOutput("dendo_cds"))
                           ),
                           sidebarLayout(
                             sidebarPanel(selectInput(inputId = "color_by_cds",
                                           label = "Color genes according to:",
                                           choices = list("organism", 
                                                          "prototype (if minmax is used)",
                                                          "Salmonella pre-selection"))
                             ),
                             mainPanel(plotOutput("cluster_wrap_cds"))
                             
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
                           )
                           
                      ),
                  
                  
                  tabPanel(title = "sRNA",
                           sidebarLayout(
                             sidebarPanel(
                               selectInput(inputId = "c_method_srna",
                                           label = "Clustering Method:",
                                           choices = list("minmax", "complete")),
                               selectInput(inputId = "d_metric_srna",
                                           label = "Distance Metric:",
                                           choices = list("cosine", "pearson")),
                               numericInput(inputId = "k_clusters_srna",
                                            label = "Number of clusters:",
                                            value = 20, min=1),
                               actionButton(inputId = "cluster_param_srna", 
                                            label = "Submit", icon = NULL, width = NULL)),
                             mainPanel(plotOutput("dendo_srna"))
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
                           )
                           ) # This tab
                ) # Tabset
              ) #fluidui

# Here I will import the data to be used for the rest of the analysis. 
# Since we only really need to import this once, we can leave it outside of
# the server function. The data I call here was performed with clustering_usw

cds_rlog.fc <- read_tsv("shiny_data/cds_logFC.tsv") %>% 
  column_to_rownames("gene_name")

srna_rlog.fc <- read_tsv("shiny_data/srna_logFC.tsv") %>% 
  column_to_rownames("gene_name")


salm_features <- read_tsv("shiny_data/merged_features_salmonella.tsv") %>% 
  mutate(organism="Salmonella",
         gene_name=paste0(locus_tag, "-Salmonella")) %>% 
  select(gene_name, locus_tag, symbol, name, feature_interval_length, product_length)

camp_features <- read_tsv("shiny_data/merged_features_campylobacter.tsv") %>% 
  mutate(organism="Campylobacter",
         gene_name=paste0(locus_tag, "-Campylobacter")) %>% 
  select(gene_name, locus_tag, symbol, name, feature_interval_length, product_length)

genomic_features <- bind_rows(salm_features, camp_features)

sussane.df <- read_tsv("shiny_data/pre_selected_salmonella.txt") %>% 
  rename("symbol"=ID) %>% 
  left_join(salm_features, by="symbol") %>% 
  filter(!is.na(locus_tag)) %>% 
  mutate(gene_name = paste0(locus_tag, "-Salmonella"))%>% 
  select(gene_name, symbol)

# Columns must be genes genes!!
correlated_distance <- function(X, method=c("cosine", "pearson")){
  
  if(method=="cosine"){
    X.t <- t(X) # Rows are genes
    sim <- X.t / sqrt(rowSums(X.t * X.t))
    sim <- sim %*% t(sim)
    
  }else if(method=="pearson"){
    sim <- cor(X, method = "pearson")
  }
  
  dist.mat <- 1 - sim
  dist.obj <- as.dist(dist.mat)
  
  return(dist.obj)
}


# Here we start our app

server <- function(input, output){
  
  # Default clustering parameters
  ## CDS
  clustering_parameters_cds <- reactiveValues(k_clusters = 20,
                                          c_method = "minmax",
                                          d_metric = "cosine",
                                          df = cds_rlog.fc,
                                          df_rn = cds_rlog.fc %>% rownames_to_column("gene_name"))
  ## sRNA
  clustering_parameters_srna <- reactiveValues(k_clusters = 20,
                                              c_method = "minmax",
                                              d_metric = "cosine",
                                              df = srna_rlog.fc,
                                              df_rn = srna_rlog.fc %>% rownames_to_column("gene_name"))
  
  # Initialize parameters
  ## CDS
  observeEvent(input$cluster_param_cds, {
    clustering_parameters_cds$k_clusters <- input$k_clusters_cds
    clustering_parameters_cds$d_metric <- input$d_metric_cds
    clustering_parameters_cds$c_method <- input$c_method_cds
    clustering_parameters_cds$df <- cds_rlog.fc
    clustering_parameters_cds$df_rn <- cds_rlog.fc %>% rownames_to_column("gene_name")
  })
  
  ## sRNA
  observeEvent(input$cluster_param_srna, {
    clustering_parameters_srna$k_clusters <- input$k_clusters_srna
    clustering_parameters_srna$d_metric <- input$d_metric_srna
    clustering_parameters_srna$c_method <- input$c_method_srna
    clustering_parameters_srna$df <- srna_rlog.fc
    clustering_parameters_srna$df_rn <- srna_rlog.fc %>% rownames_to_column("gene_name")
  })
  
  # Clustering  
  ## CDS
  clust_obj_cds <- reactive({
    dist_obj <- correlated_distance(t(clustering_parameters_cds$df), 
                                    method = clustering_parameters_cds$d_metric)
    
    if(clustering_parameters_cds$c_method == "minmax"){
      cobj <- protoclust(dist_obj)
    }else if(clustering_parameters_cds$c_method == "complete"){
      cobj <- hclust(dist_obj, method = "complete")
    }
    cobj
  })
  
  ## sRNA
  clust_obj_srna <- reactive({
    dist_obj <- correlated_distance(t(clustering_parameters_srna$df), 
                                    method = clustering_parameters_srna$d_metric)
    
    if(clustering_parameters_srna$c_method == "minmax"){
      cobj <- protoclust(dist_obj)
    }else if(clustering_parameters_srna$c_method == "complete"){
      cobj <- hclust(dist_obj, method = "complete")
    }
    cobj
  })
  
  # Dendogram
  ## CDS
  output$dendo_cds <- renderPlot({
    plot(clust_obj_cds(), labels=rep("", nrow(clustering_parameters_cds$df)))
    rect.hclust(clust_obj_cds(), k = clustering_parameters_cds$k_clusters)
    })
  
  ## sRNA
  output$dendo_srna <- renderPlot({
    plot(clust_obj_srna(), labels=rep("", nrow(clustering_parameters_srna$df)))
    rect.hclust(clust_obj_srna(), k = clustering_parameters_srna$k_clusters)
  })
  
  
  # Cluster Data Frame
  ## CDS
  cluster_info_cds <- reactive({
    if(clustering_parameters_cds$c_method == "minmax"){
      cut <- protocut(clust_obj_cds(), k=clustering_parameters_cds$k_clusters)
    }else if(clustering_parameters_cds$c_method == "complete"){
      cut <- list("cl"=cutree(clust_obj_cds(), k=clustering_parameters_cds$k_clusters))
    }
    
    c.df <- data.frame("cluster"=cut$cl,
                               "gene_name"=names(cut$cl),
                               "organism"=sapply(names(cut$cl), function(x){tail(strsplit(x, "-")[[1]],n=1)}),
                               "prototype"="Not Prototype")
    
    if(clustering_parameters_cds$c_method == "minmax"){
      c.df[cut$protos, "prototype"] <- "Prototype"
    }
    
    c.df.expression <- clustering_parameters_cds$df_rn %>% 
      left_join(c.df, by="gene_name") %>% 
      gather(Ctrl, As, Bs, Hyp, Li, Nd, Ns, Oss, Oxs, Sp, Tm, Vic,
             key = "Condition", value="logFC") %>% 
      mutate(Condition=factor(Condition, levels = c("Ctrl", "As", "Bs", "Hyp", "Li",
                                                    "Nd", "Ns", "Oss", "Oxs", "Sp",
                                                    "Tm", "Vic")))
    
    c.df.expression
    
  })
  
  ## sRNA
  cluster_info_srna <- reactive({
    if(clustering_parameters_srna$c_method == "minmax"){
      cut <- protocut(clust_obj_srna(), k=clustering_parameters_srna$k_clusters)
    }else if(clustering_parameters_srna$c_method == "complete"){
      cut <- list("cl"=cutree(clust_obj_srna(), k=clustering_parameters_srna$k_clusters))
    }
    
    c.df <- data.frame("cluster"=cut$cl,
                       "gene_name"=names(cut$cl),
                       "organism"=sapply(names(cut$cl), function(x){tail(strsplit(x, "-")[[1]],n=1)}),
                       "prototype"="Not Prototype")
    
    if(clustering_parameters_srna$c_method == "minmax"){
      c.df[cut$protos, "prototype"] <- "Prototype"
    }
    
    c.df.expression <- clustering_parameters_srna$df_rn %>% 
      left_join(c.df, by="gene_name") %>% 
      gather(Ctrl, As, Bs, Hyp, Li, Nd, Ns, Oss, Oxs, Sp, Tm, Vic,
             key = "Condition", value="logFC") %>% 
      mutate(Condition=factor(Condition, levels = c("Ctrl", "As", "Bs", "Hyp", "Li",
                                                    "Nd", "Ns", "Oss", "Oxs", "Sp",
                                                    "Tm", "Vic")))
    
    c.df.expression
    
  })
  
  # Cluster Profiles
  ## CDS
  output$cluster_wrap_cds <- renderPlot({
    if(input$color_by_cds == 'organism'){
      feat_campy <- cluster_info_cds() %>% filter(organism=="Campylobacter")
      feat_salm <- cluster_info_cds() %>% filter(organism=="Salmonella")
  
      
      ggplot() +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data = feat_salm,
                  colour=alpha("grey", 0.5)) +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data=feat_campy,
                  colour=alpha("steelblue", 0.5)) +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~cluster) +
        ggtitle("CDS logFC by Cluster", subtitle = "Campylobacter in blue")
      
    }else if(input$color_by_cds == "prototype (if minmax is used)"){
      gprototypes <- cluster_info_cds() %>% filter(prototype=="Prototype")
      gnormal <- cluster_info_cds() %>% filter(prototype=="Not Prototype")
      
      ggplot() +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data = gnormal,
                  colour=alpha("grey", 0.5)) +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data=gprototypes,
                  colour="red") +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~cluster) +
        ggtitle("CDS logFC by Cluster", subtitle = "Prototype in red")
      
    }else if(input$color_by_cds == "Salmonella pre-selection"){
      feat_salm <- cluster_info_cds() %>% filter(organism=="Salmonella")
      salm.sussane <- feat_salm %>% 
        filter(gene_name %in% sussane.df$gene_name) %>% 
        left_join(sussane.df, by="gene_name")
      
      ggplot() +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data = feat_salm,
                  colour=alpha("grey", 0.5)) +
        geom_line(aes(x=Condition, y=logFC, group=gene_name, color=symbol), 
                  data=salm.sussane) +
        geom_text(aes(x=Condition, y=logFC + 0.2, label=symbol),
                  data = salm.sussane %>% filter(Condition=="Vic"),
                  size=1.8) +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = "none") +
        facet_wrap(~cluster) +
        ggtitle("CDS logFC by cluster", subtitle = "Pre-selection of reporters")
      }
  })
  
  ## sRNA
  output$cluster_wrap_srna <- renderPlot({
    if(input$color_by_srna == 'organism'){
      feat_campy <- cluster_info_srna() %>% filter(organism=="Campylobacter")
      feat_salm <- cluster_info_srna() %>% filter(organism=="Salmonella")
      
      
      ggplot() +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data = feat_salm,
                  colour=alpha("grey", 0.5)) +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data=feat_campy,
                  colour=alpha("steelblue", 0.5)) +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~cluster) +
        ggtitle("sRNA logFC by Cluster", subtitle = "Campylobacter in blue")
      
    }else if(input$color_by_srna == "prototype (if minmax is used)"){
      gprototypes <- cluster_info_srna() %>% filter(prototype=="Prototype")
      gnormal <- cluster_info_srna() %>% filter(prototype=="Not Prototype")
      
      ggplot() +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data = gnormal,
                  colour=alpha("grey", 0.5)) +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data=gprototypes,
                  colour="red") +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_wrap(~cluster) +
        ggtitle("sRNA logFC by Cluster", subtitle = "Prototype in red")
      
    }else if(input$color_by_srna == "Salmonella pre-selection"){
      feat_salm <- cluster_info_srna() %>% filter(organism=="Salmonella")
      salm.sussane <- feat_salm %>% 
        filter(gene_name %in% sussane.df$gene_name) %>% 
        left_join(sussane.df, by="gene_name")
      
      ggplot() +
        geom_line(aes(x=Condition, y=logFC, group=gene_name), data = feat_salm,
                  colour=alpha("grey", 0.5)) +
        geom_line(aes(x=Condition, y=logFC, group=gene_name, color=symbol), 
                  data=salm.sussane) +
        geom_text(aes(x=Condition, y=logFC + 0.2, label=symbol),
                  data = salm.sussane %>% filter(Condition=="Vic"),
                  size=1.8) +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = "none") +
        facet_wrap(~cluster) +
        ggtitle("sRNA logFC by cluster", subtitle = "Pre-selection of reporters")
    }
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
    
    if(focus_parameters_cds$o_focus == "Both organisms"){
      f.df <- cluster_info_cds() %>% 
        filter(cluster == focus_parameters_cds$c_focus) 
      
    }else if(focus_parameters_cds$o_focus == "Only Salmonella"){
      f.df <- cluster_info_cds() %>% 
        filter(cluster == focus_parameters_cds$c_focus) %>% 
        filter(organism == "Salmonella")
    }else if(focus_parameters_cds$o_focus == "Only Campylobacter"){
      f.df <- cluster_info_cds() %>% 
        filter(cluster == focus_parameters_cds$c_focus) %>% 
        filter(organism == "Campylobacter")
    }
    
    f.plot <- ggplot(f.df, aes(x=Condition, y=logFC, group=gene_name)) +
                         geom_line(color=alpha("darkgrey", 0.5))
    
    if(focus_parameters_cds$g_focus != "None"){
      f.plot <- f.plot + 
        geom_line(data=f.df %>% filter(gene_name %in% focus_parameters_cds$g_focus),
                  aes(x=Condition, y=logFC, group=gene_name, color=gene_name))
    }
    ggplotly(f.plot)
  })
  
  ## sRNA
  output$focus_profile_srna <- renderPlotly({
    
    if(focus_parameters_srna$o_focus == "Both organisms"){
      f.df <- cluster_info_srna() %>% 
        filter(cluster == focus_parameters_srna$c_focus) 
      
    }else if(focus_parameters_srna$o_focus == "Only Salmonella"){
      f.df <- cluster_info_srna() %>% 
        filter(cluster == focus_parameters_srna$c_focus) %>% 
        filter(organism == "Salmonella")
    }else if(focus_parameters_srna$o_focus == "Only Campylobacter"){
      f.df <- cluster_info_srna() %>% 
        filter(cluster == focus_parameters_srna$c_focus) %>% 
        filter(organism == "Campylobacter")
    }
    
    f.plot <- ggplot(f.df, aes(x=Condition, y=logFC, group=gene_name)) +
      geom_line(color=alpha("darkgrey", 0.5))
    
    if(focus_parameters_srna$g_focus != "None"){
      f.plot <- f.plot + 
        geom_line(data=f.df %>% filter(gene_name %in% focus_parameters_srna$g_focus),
                  aes(x=Condition, y=logFC, group=gene_name, color=gene_name))
    }
    ggplotly(f.plot)
  })
  
  
  # Focused Data
  ## CDS
  output$cluster_members_cds <- renderDT({
    f.df <- cluster_info_cds() %>% 
      filter(cluster == focus_parameters_cds$c_focus)
    
    
    g_feat <- genomic_features %>% 
      filter(gene_name %in% f.df$gene_name) %>% 
      mutate("Pre-selected"=if_else(gene_name %in% sussane.df$gene_name, "Pre-selected", ""))
    
    if(clustering_parameters_srna$c_method == "minmax"){
      gprototypes <- f.df %>% filter(prototype=="Prototype")
      
      g_feat <- g_feat %>% 
        mutate("Prototype" = if_else(gene_name %in% gprototypes$gene_name, "Prototype", ""))
    }
    
    g_feat
    
  })
  
  ## sRNA
  output$cluster_members_srna <- renderDT({
    f.df <- cluster_info_srna() %>% 
      filter(cluster == focus_parameters_srna$c_focus)
    
    
    g_feat <- genomic_features %>% 
      filter(gene_name %in% f.df$gene_name) %>% 
      mutate("Pre-selected"=if_else(gene_name %in% sussane.df$gene_name, "Pre-selected", ""))
    
    if(clustering_parameters_srna$c_method == "minmax"){
      gprototypes <- f.df %>% filter(prototype=="Prototype")
      
      g_feat <- g_feat %>% 
        mutate("Prototype" = if_else(gene_name %in% gprototypes$gene_name, "Prototype", ""))
    }
  })
  
}

shinyApp(ui = ui, server = server)
