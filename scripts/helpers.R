# This is going to be the helper functions for the app

library(tidyverse)
library(protoclust)
library(plotly)


# A function that decides the matrix and then runs the correlated distance function
cluster_counts <- function(what_to_cluster, distance_metric, how_to_cluster, logfc_df, posc_df){
  
  # Decide which df should be clustered
  if(what_to_cluster == "logfc"){
    df_input <- logfc_df
    expect_negatives <- TRUE
  }else if(what_to_cluster == "poscounts"){
    df_input <- posc_df
    expect_negatives <- FALSE
  }
  
  # Get the correct distances
  d.output <- correlated_distance(X=t(df_input), method = distance_metric, expect_neg_distance = expect_negatives)
  dist_obj <- d.output$dist_obj
  
  # Now perform clustering
  if(how_to_cluster == "minmax"){
    cobj <- protoclust(dist_obj)
  }else{
    cobj <- hclust(dist_obj, method = how_to_cluster)
  }
  
  return(list("cluster_object"=cobj, "distance_matrix"=d.output$dist_mat))
}

# Columns must be genes genes!!
correlated_distance <- function(X, method=c("cosine", "pearson"), expect_neg_distance=FALSE){
  
  if(method=="cosine"){
    X.t <- t(X) # Rows are genes
    sim <- X.t / sqrt(rowSums(X.t * X.t))
    sim <- sim %*% t(sim)
    diag(sim) <- 1 # R is a little inefficient sometimes
    
  }else{
    sim <- cor(X, method = method)
  }
  
  dist.mat <- 1 - sim
  
  if(expect_neg_distance){
    dist.mat <- sqrt(dist.mat * 0.5) # comes from https://arxiv.org/pdf/1208.3145.pdf
  }
  
  
  dist.obj <- as.dist(dist.mat)
  
  return(list("dist_obj"=dist.obj, "dist_mat"=dist.mat))
}


# Helping out to plot the deprogram
plot_dendrogram <- function(how_to_cluster, number_of_clusters, highlight_this_cluster, clustered_object){
  
  # First we assign genes to clusters 
  cut <-  get_cut(how_to_cluster = how_to_cluster, 
                  number_of_clusters=number_of_clusters,
                  clustered_object = clustered_object)
  
  
  
  # Then we plot according to the cluster we want to highlight
  if(highlight_this_cluster == 0){
    
    plot(clustered_object, labels=rep("", length(cut$cl)))
    rect.hclust(clustered_object, k = number_of_clusters) 
    
  }else{
    
    plot(clustered_object, 
         labels=sapply(cut$cl, function(x){if(x==highlight_this_cluster){return("*")}else{return("")}}))
    rect.hclust(clustered_object, k = number_of_clusters)
    
  }
}

my_wss <- function(agg.by, what_to_aggregate){
  x.SS <- aggregate(what_to_aggregate, by=list(agg.by), function(x) sum(scale(x, scale=FALSE)^2))  
  sum(x.SS[,-1])
}

plot_wss <- function(clustered_object, what_was_clustered, how_was_it_clustered, 
                     number_of_clusters, logfc_df, posc_df){
  
  if(how_was_it_clustered == "minmax"){
    hcuts <- bind_cols(lapply(2:50, function(x){cas <- protocut(clustered_object, k=x)$cl; data.frame("cluster_assignment" = cas)}))
    colnames(hcuts) <- as.character(2:50)
  }else{
    hcuts <- cutree(clustered_object, k=2:50)
  }
  
  if(what_was_clustered=="logfc"){
    input.df <- logfc_df
  }else if(what_was_clustered == "poscounts"){
    input.df <- posc_df
  }
  
  tibble(wss = apply(hcuts, 2, my_wss, what_to_aggregate=input.df),
         nk = 2:50) %>% 
    ggplot(aes(x=nk, y=wss)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = number_of_clusters, color="red", linetype="dashed") +
    labs(x="Number of clusters",
         y="WSS",
         title="Total Within-group Sum of Squared distance")
  
  
}

# Mini help for cut object
get_cut <- function(how_to_cluster, number_of_clusters, clustered_object){
  
  if(how_to_cluster == "minmax"){
    cut <- protocut(clustered_object, k=number_of_clusters)
  }else{
    cut <- list("cl"=cutree(clustered_object, k=number_of_clusters))
  }
  
  return(cut)
}

# Dataframe of clustering

cluster_longform <- function(how_to_cluster, number_of_clusters, clustered_object, logfc_dataframe){
  # Get assignment of clusters
  cut <- get_cut(how_to_cluster = how_to_cluster,
                 number_of_clusters = number_of_clusters,
                 clustered_object = clustered_object)
  
  # Begin creating the DF
  c.df <- data.frame("cluster"=cut$cl,
                     "gene_name"=names(cut$cl),
                     "organism"=sapply(names(cut$cl), function(x){tail(strsplit(x, "-")[[1]],n=1)}),
                     "prototype"="Not Prototype")
  
  # Assign prototypes
  if(how_to_cluster == "minmax"){
    c.df[cut$protos, "prototype"] <- "Prototype"
  }
  
  c.df.expression <- logfc_dataframe %>% 
    rownames_to_column("gene_name") %>% 
    left_join(c.df, by="gene_name") %>% 
    gather(Ctrl, As, Bs, Hyp, Li, Nd, Ns, Oss, Oxs, Sp, Tm, Vic,
           key = "Condition", value="logFC") %>% 
    mutate(Condition=factor(Condition, levels = c("Ctrl", "As", "Bs", "Hyp", "Li",
                                                  "Nd", "Ns", "Oss", "Oxs", "Sp",
                                                  "Tm", "Vic")))
  
  return(c.df.expression)
  
}


cluster_profile <- function(dataframe_for_plot, plot_according_to,
                            salmonella_preselection, campylobacter_preselection){
  if(plot_according_to == 'organism'){
    feat_campy <- dataframe_for_plot %>% filter(organism=="Campylobacter")
    feat_salm <- dataframe_for_plot %>% filter(organism=="Salmonella")
    
    
    ggplot() +
      geom_line(aes(x=Condition, y=logFC, group=gene_name), data = feat_salm,
                colour=alpha("grey", 0.5)) +
      geom_line(aes(x=Condition, y=logFC, group=gene_name), data=feat_campy,
                colour=alpha("steelblue", 0.5)) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~cluster) +
      ggtitle("CDS logFC by Cluster", subtitle = "Campylobacter in blue")
    
  }else if(dataframe_for_plot == "prototype (if minmax is used)"){
    gprototypes <- dataframe_for_plot %>% filter(prototype=="Prototype")
    gnormal <- dataframe_for_plot %>% filter(prototype=="Not Prototype")
    
    ggplot() +
      geom_line(aes(x=Condition, y=logFC, group=gene_name), data = gnormal,
                colour=alpha("grey", 0.5)) +
      geom_line(aes(x=Condition, y=logFC, group=gene_name), data=gprototypes,
                colour="red") +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~cluster) +
      ggtitle("CDS logFC by Cluster", subtitle = "Prototype in red")
    
  }else if(plot_according_to == "Salmonella pre-selection"){
    feat_salm <- dataframe_for_plot %>% filter(organism=="Salmonella")
    salm.sussane <- feat_salm %>% 
      filter(gene_name %in% salmonella_preselection$gene_name) %>% 
      left_join(salmonella_preselection, by="gene_name")
    
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
    
  }else if(plot_according_to == "Campylobacter pre-selection"){
    feat_camp <- dataframe_for_plot %>% 
      filter(organism=="Campylobacter")
    
    camp.sarah  <- feat_camp %>% 
      filter(gene_name %in% campylobacter_preselection$gene_name) %>% 
      left_join(campylobacter_preselection, by="gene_name") %>% 
      mutate(clone_candidate = if_else(is.na(`clone?`), "Not", "yes"))
      
    
    ggplot() +
      geom_line(aes(x=Condition, y=logFC, group=gene_name), data = feat_camp,
                colour=alpha("grey", 0.5)) +
      geom_line(aes(x=Condition, y=logFC, group=gene_name, color=clone_candidate),
                data=camp.sarah) +
      geom_text(aes(x=Condition, y=logFC + 0.2, label=gene),
                data = camp.sarah %>% filter(Condition=="Vic"),
                size=1.8) +
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "none") +
      scale_color_brewer(palette = "Dark2") +
      facet_wrap(~cluster) +
      ggtitle("CDS logFC by cluster", subtitle = "Genes of interest are green. Genes for cloning are brown")
    
  }
}


focused_profile <- function(which_organism, which_cluster, which_gene, 
                            dataframe_for_plot){
  
  if(which_organism == "Both organisms"){
    f.df <- dataframe_for_plot %>% 
      filter(cluster == which_cluster) 
    
  }else if(which_organism == "Only Salmonella"){
    f.df <- dataframe_for_plot %>% 
      filter(cluster == which_cluster) %>% 
      filter(organism == "Salmonella")
    
  }else if(which_organism == "Only Campylobacter"){
    f.df <- dataframe_for_plot %>% 
      filter(cluster == which_cluster) %>% 
      filter(organism == "Campylobacter")
  }
  
  f.plot <- ggplot(f.df, aes(x=Condition, y=logFC, group=gene_name)) +
    geom_line(color=alpha("darkgrey", 0.5))
  
  if(which_gene != "None"){
    f.plot <- f.plot + 
      geom_line(data=f.df %>% filter(gene_name %in% which_gene),
                aes(x=Condition, y=logFC, group=gene_name, color=gene_name))
  }
  ggplotly(f.plot)
}

cluster_table <- function(dataframe_for_info, which_cluster, distance_matrix,
                          genomic_features_dataframe,
                          salmonella_preselection,
                          campylobacter_preselection,
                          average_expression_information,
                          how_to_cluster){
  
  f.df <- dataframe_for_info %>% 
    filter(cluster == which_cluster)
  
  cluster_members <- f.df %>% 
    select(gene_name) %>% 
    unique()
  
  representative_info <- gene_distances(cluster_members$gene_name,
                                        distance_matrix)
  
  g_feat <- genomic_features_dataframe %>% 
    filter(gene_name %in% cluster_members$gene_name) %>% 
    mutate("Pre-selected"=if_else((gene_name %in% salmonella_preselection$gene_name) | (gene_name %in% campylobacter_preselection$gene_name), "Pre-selected", "")) %>%
    left_join(representative_info, by="gene_name") %>% 
    left_join(average_expression_information, by="gene_name")
  
  if(how_to_cluster == "minmax"){
    gprototypes <- f.df %>% filter(prototype=="Prototype")
    
    g_feat <- g_feat %>% 
      mutate("Prototype" = if_else(gene_name %in% gprototypes$gene_name, "Prototype", ""))
  }
  
  g_feat %>% 
    select(-c(species_max_distance,species_avg_distance))
}


# Some info on distances
species_distances <- function(species, mem.df, d){
  spec_members <- mem.df %>% 
    filter(organism == species) %>% 
    select(gene_name)
  
  if(length(spec_members$gene_name) <= 1){
    return(data.frame("gene_name"=spec_members$gene_name,
                      "species_max_distance"=0,
                      "species_avg_distance"=0))
  }
  
  spec_dist <- d[spec_members$gene_name, spec_members$gene_name]
  spec_max_distances <- apply(spec_dist, 1, max)
  spec_avg_distances <- rowSums(spec_dist) / (ncol(spec_dist) - 1)
  
  
  s.df <- data.frame("gene_name"=names(spec_max_distances),
                     "species_max_distance"=spec_max_distances,
                     "species_avg_distance"=spec_avg_distances)
  
  return(s.df)
  
}

# Not too sure about this
gene_distances <- function(members, dmat){
  
  if(length(members) <= 1){
    return(data.frame("gene_name"=members,
                      "cluster_max_distance"=0,
                      "cluster_avg_distance"=0,
                      "species_max_distance"=0,
                      "species_avg_distance"=0))
  }
  
  clus_dist <- dmat[members, members]
  clus_max_distances <- apply(clus_dist, 1, max)
  clus_avg_distances <- rowSums(clus_dist) / (ncol(clus_dist) - 1) # Se we can ignore the 0 of the diagonal
  
  o.df <- data.frame("gene_name"=names(clus_max_distances),
                     "cluster_max_distance"=clus_max_distances,
                     "cluster_avg_distances"=clus_avg_distances,
                     "organism"=sapply(names(clus_max_distances), function(x){tail(strsplit(x,"-")[[1]], n=1)}))
  
  sp.df <- bind_rows(lapply(unique(o.df$organism), species_distances, o.df, clus_dist))
  
  distance_info <- o.df %>% 
    left_join(sp.df, by="gene_name") %>% 
    select(-organism)
  
  return(distance_info)
}

complete_table <- function(dataframe_with_info, genomic_features_dataframe){
  dataframe_with_info %>% 
    select(gene_name, cluster, organism) %>% 
    distinct() %>% 
    left_join(genomic_features_dataframe, by="gene_name")
  
}