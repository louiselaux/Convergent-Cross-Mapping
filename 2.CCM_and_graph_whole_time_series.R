
# Load the library
library(dplyr)
library(igraph)
library(ggplot2)

# Data
env_window<- residuals_table_ccm

# Dataframe to store the result
final_dataframe <- data.frame(variable = NA, target = NA,
                              Sign_OK = NA, Convergence_OK=NA, Rho_Real = NA, Mean_Rho_Surr = NA,
                              Rho_Diff = NA, P_Value_Surr = NA,
                              Rho_20 = NA, Rho_50 = NA, Rho_60 = NA, Rho_80 = NA,
                              
                              mean_rho_20=NA, mean_rho_80=NA, p_value_conv=NA)



#add an empty list
rho_pairs <- list()
# Loop on each variable
for (columns in colnames(env_window)) {
  E <- embed_stats_tot$E_best[embed_stats_tot$variables == columns]  # Extract the best E for each causee variable
  
  # All possible causal effects
  targetlist <- setdiff(colnames(env_window), columns)
  
  # For each causal effect that might exist
  for (target in targetlist) {
    # Correlation between the variable and the target
    corr_value <- cor(env_window[[columns]], env_window[[target]], use = "complete.obs")
    
    # CCM
    datatable <- CCM(
      dataFrame = env_window,
      E = E,
      Tp = -1,
      columns = columns,
      target = target,
      libSizes = "20 50 60 80",
      sample = 100,
      noTime = TRUE,
      showPlot = FALSE,
      includeData = TRUE
    )
    
    # Extract rho value for each lib size
    dt <- as.data.frame(datatable$CCM1_PredictStat)
    
    rho_20 <- dt %>% filter(LibSize == "20") %>% pull(rho)
    rho_50 <- dt %>% filter(LibSize == "50") %>% pull(rho)
    rho_60 <- dt %>% filter(LibSize == "60") %>% pull(rho)
    rho_80 <- dt %>% filter(LibSize == "80") %>% pull(rho)
    
    # Mean of rho
    mean_rho_20 <- mean(rho_20, na.rm = TRUE)
    mean_rho_50 <- mean(rho_50, na.rm = TRUE)
    mean_rho_60 <- mean(rho_60, na.rm = TRUE)
    mean_rho_80 <- mean(rho_80, na.rm = TRUE)
    
    # Test t to compare rho 80 and rho 20
    if (length(rho_80) > 1 && length(rho_20) > 1) {
      t_test_result <- t.test(rho_50, rho_20, paired = TRUE)
      convergence_ok <- t_test_result$p.value < 0.05  # Convergence is significant if the p-value is lower than 0.05
      p_value_conv<- t_test_result$p.value
    } else {
      convergence_ok <- FALSE
    }
    
    # Store the rhos
    rho_pairs[[paste(columns, target, sep = "_")]] <- rho_80
    
    # Initialization of surrogate values
    rho_diff <- NA
    p_value_surr <- NA
    mean_rho_surr <- NA
    
    # If convergence ok, generate surrogates
    if (convergence_ok) {
      T_period <- 26
      ts_data <- ts(env_window[[columns]], frequency = T_period)
      surrogates <- SurrogateData(ts = ts_data, method = "seasonal", T_period = T_period, num_surr = 100)
      
      rho_surrogates <- numeric(100)
      
      for (i in 1:100) {
        surrogate_data <- surrogates[, i]
        env_window_sur <- env_window
        env_window_sur[[paste0(columns, "_surr")]] <- surrogate_data
        
        surrogate_ccm <- CCM(
          dataFrame = env_window_sur,
          E = E,
          Tp = -1,
          columns = paste0(columns, "_surr"),
          target = target,
          libSizes = "80",
          sample = 1,
          noTime = TRUE,
          showPlot = FALSE,
          includeData = TRUE
        )
        
        dt_surr <- as.data.frame(surrogate_ccm$CCM1_PredictStat)
        rho_surrogates[i] <- mean(dt_surr$rho, na.rm = TRUE)
      }
      
      mean_rho_surr <- mean(rho_surrogates, na.rm = TRUE)
      num_surr_higher <- sum(rho_80 > rho_surrogates, na.rm = TRUE)
      p_value_surr <- 1 - (num_surr_higher / 100)
      rho_diff <- rho_80 - mean_rho_surr
    }
    
    # Store the results
    final_dataframe <- rbind(final_dataframe, data.frame(
      variable = columns, target = target,
      Sign_OK = convergence_ok, Convergence_OK = convergence_ok,
      Rho_Real = rho_80, Mean_Rho_Surr = mean_rho_surr,
      Rho_Diff = rho_diff, P_Value_Surr = p_value_surr,
      Rho_20 = rho_20, Rho_50 = rho_50, Rho_60 = rho_60, Rho_80 = rho_80,mean_rho_20=mean_rho_20, mean_rho_80=mean_rho_80, p_value_conv= p_value_conv,
      #Increasing_Rho = (rho_80 > rho_60) & (rho_60 > rho_50) & (rho_50 > rho_20),
      stringsAsFactors = FALSE
    ))
  }
}

# Comparison of directions
major_effects <- data.frame(variable = NA, target = NA, p_value = NA, direction = NA, stringsAsFactors = FALSE)



for (pair in names(rho_pairs)) {
  vars <- strsplit(pair, "_")[[1]]
  
  if (length(vars) == 2) {
    X <- vars[1] #column
    Y <- vars[2] #target
    
    if (paste(Y, X, sep = "_") %in% names(rho_pairs)) {
      rho_xy <- unlist(rho_pairs[[pair]])  #
      rho_yx <- unlist(rho_pairs[[paste(Y, X, sep = "_")]])  #
      
      #
      if (length(rho_xy) > 1 && length(rho_yx) > 1) {
        t_test_direction <- t.test(rho_xy, rho_yx, paired = TRUE)
        
        if (t_test_direction$p.value < 0.05) {  #
          if (mean(rho_xy, na.rm = TRUE) > mean(rho_yx, na.rm = TRUE)) {
            major_effects <- rbind(major_effects, data.frame(variable = X, target = Y, p_value = t_test_direction$p.value, direction = paste(Y, "-->", X), stringsAsFactors = FALSE))
          } else {
            major_effects <- rbind(major_effects, data.frame(variable = Y, target = X, p_value = t_test_direction$p.value, direction = paste(X, "-->", Y), stringsAsFactors = FALSE))
          }
        }
      }
    }
  }
}


# Final filtering : keep only relation with convergence significant and  P-value_surr < 0.05
filtered_effects <- major_effects %>%
  inner_join(final_dataframe, by = c("variable", "target")) %>%
  filter(p_value_conv<0.00001, P_Value_Surr <= 0.011, p_value <0.05) #does not work

#
edges <- filtered_effects %>%
  select(target,variable)

graph <- graph_from_data_frame(edges, directed = TRUE)


#E(graph)$arrow.size <- edges$edge_size

# Plot of the graph
plot(graph, vertex.size = 18, vertex.label.color = "black",
     vertex.color = "lightblue",
     #edge.arrow.size = E(graph)$arrow.size,
     # edge.arrow.width = E(graph)$arrow.width,
     edge.arrow.size = 1,  #
     edge.arrow.width = 1,  #
     edge.color = "darkgrey",  #
     edge.width = 1,  #
     edge.lty = 1,  #
     edge.curved = 0.2,  #
     main = "Main causal relationship")


#Other try
#
filtered_effects <- major_effects %>%
  inner_join(final_dataframe, by = c("variable", "target")) %>%
  mutate(edge_style = ifelse(P_Value_Surr < 0.05, 1, 2)) %>%  #
  filter(Convergence_OK == TRUE)%>%mutate(edge_size=)  #

#
edges <- filtered_effects %>%
  select(target, variable, edge_style)

#
graph <- graph_from_data_frame(edges, directed = TRUE)

#
E(graph)$lty <- edges$edge_style  # Style des arêtes (1 = solide, 2 = pointillés)

#
plot(graph, vertex.size = 18, vertex.label.color = "black",
     vertex.color = "lightblue",
     edge.arrow.size = 1.5,  #
     edge.arrow.width = 2,  #
     edge.color = "darkgrey",  #
     edge.width = 2,
     edge.lty = E(graph)$lty,
     edge.curved = 0.2,
     main = "Main causal relationship (solid = significant, dashed = non-significant)")

#pred en haut
# pas convaincue
layout_tree <- layout_as_tree(graph, root = which(V(graph)$name %in% c("logpred","logsmallcarni")))

# Plot du graphe avec le layout
plot(graph, layout = layout_tree, vertex.size = 18, vertex.label.color = "black",
     vertex.color = "lightblue",
     edge.arrow.size = 1.5, edge.arrow.width = 2,
     edge.color = "darkgrey", edge.width = 2, edge.lty = E(graph)$lty,
     edge.curved = 0.2, main = "Main causal relationship")
