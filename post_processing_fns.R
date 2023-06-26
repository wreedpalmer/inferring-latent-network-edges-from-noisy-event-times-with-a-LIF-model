# function for converting matrix to data frame
# diag.rm=T indicates an input matrix with diagonal removed
matrix2df <- function(mat_X, val_name_to="value", row_name_to="row", col_name_to="col", diag.rm=F){
  lapply(
    1:ncol(mat_X),
    function(col_j){
      lapply((1:nrow(mat_X)),
             function(row_i){
               tibble(
                 !!sym(row_name_to) := ifelse(diag.rm & row_i >= col_j, row_i + 1, row_i),
                 !!sym(col_name_to) := col_j,
                 !!sym(val_name_to) := mat_X[row_i, col_j]
               )
             }) %>%
        bind_rows
    }
  ) %>%
    bind_rows
}

# function for list of same-sized matrices to data frame, with option of ignoring diagonal entries
# "step" variable indicates place in list
# diag.rm=T indicates an input matrix with diagonal removed
#### NOT EFFICIENT !!!!
matrixList2df <- function(mat_list, val_name_to="value", row_name_to="row", col_name_to="col", diag.rm=F){
  lapply(
    1:ncol(mat_list[[1]]),
    function(col_j){
      lapply((1:nrow(mat_list[[1]])),
             function(row_i){
               tibble(
                 !!sym(val_name_to) := vapply(mat_list, function(mat_X) mat_X[row_i, col_j], 0),
                 !!sym(row_name_to) := ifelse(diag.rm & row_i >= col_j, row_i + 1, row_i),
                 !!sym(col_name_to) := col_j
               ) %>%
                 mutate(step = row_number())
             }) %>%
        bind_rows
    }
  ) %>%
    bind_rows
}

### extracting / organizing estimates for evaluation / visualization
get_indices <- function(target_indices, source_index){
  ifelse(target_indices < source_index, target_indices, target_indices - 1)
}

get_W_est <- function(runList){
  n <- runList$model_params$n
  W_est <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    W_est[setdiff(1:n, i), i] <-
      runList$adam_update$values_out$W[get_indices(setdiff(1:n, i), i), i]
  }
  return(W_est)
}

get_W_paths_df <- function(runList, run_name = NULL, sim = FALSE, mat_list = NULL){
  if(is.null(mat_list)){
    mat_list = lapply(runList$adam_update$values_epoch, "[[",  "W")
  }
  
  out_df <- matrixList2df(
    mat_list = mat_list,
    val_name_to = "Wij_est",
    row_name_to = "source_node",
    col_name_to = "target_node",
    diag.rm = TRUE
  )
  
  if(!is.null(run_name)){
    out_df$run_name <- run_name
  }
  
  # add the known edge weights if dealing with simulation fit
  if (sim) {
    out_df <- out_df %>%
      left_join(
        matrix2df(
          runList$model_params$W_true, val_name_to="Wij_true", row_name_to="source_node", col_name_to="target_node", diag.rm=FALSE
        )
      ) %>%
      mutate(
        sqrd_error = (Wij_true - Wij_est) ^ 2
      )
  }
  
  return(out_df)
}

get_non_W_param_updates_df <- function(runList, run_name = NULL, sim = FALSE){
  
  out_df <- lapply(
    seq_along(runList$adam_update$values_epoch),
    function(update_step){
      runList$adam_update$values_epoch[[update_step]]$phi %>% t() %>%
        cbind("eta" = runList$adam_update$values_epoch[[update_step]]$eta) %>%
        as_tibble() %>% mutate(node_i = row_number(), step = update_step - 1)
    }
  ) %>%
    bind_rows() %>%
    mutate(sigma_tau_ratio = sigma / tau) %>%
    pivot_longer(cols = c(1:8, 11), names_to = "variable") %>%
    mutate(
      variable = factor(variable,
                        levels = c("eta", "sigma",
                                   setdiff(rownames(runList$adam_update$values_out$phi), "sigma"),
                                   "sigma_tau_ratio")
                        )
      )
  
  if(!is.null(run_name)){
    out_df$run_name <- run_name
  }
    
  # add known true values for eta and sigma if dealing with simulation fit
  if (sim) {
    out_df <- out_df %>%
      left_join(
        tibble(
          node_i = 1:runList$model_params$n,
          sigma = runList$model_params$sigma_vec,
          eta = runList$model_params$eta_vec
        ) %>%
          pivot_longer(cols = 2:3, names_to = "variable", values_to = "true_value")
      ) %>%
      mutate(
        sqrd_error = (true_value - value) ^ 2
      )
  }
  
  return(out_df)
}



### function for plotting adjacency matrix
getHeatMap <- function(A, title="", annotate=FALSE, black_white = T){
  ids = colnames(A)
  dimnames(A) <- NULL
  indices = 1:nrow(A)
  
  A.melt <- reshape2::melt(A[indices,indices]) %>%
    mutate(from = ids[indices][Var1],
           to = ids[indices][Var2],
           Var1 = -1*Var1)
  
  
  p <- ggplot(data = A.melt, aes(x=Var2, y=Var1, fill=value, label=round(value,1))) + 
    geom_tile() +
    coord_fixed() +
    labs(col=NULL, title=title, x=NULL, y=NULL) +
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      #legend.position="none",
      plot.title=element_text(hjust=0.5),
      plot.subtitle=element_text(hjust=0.5),
      #axis.text = element_blank(),
      axis.ticks = element_blank(),
      #axis.title = element_blank(),
      panel.background = element_rect(fill = "lightblue",
                                      colour = "lightblue")
    ) +
    scale_x_continuous(breaks = 1:20) +
    scale_y_continuous(breaks = -1:-20, labels = paste(1:20))
  if(black_white){
    p <- p + scale_fill_gradient(low = "white", high = "black")
  }else{
    p <- p + #scale_fill_gradient(low = "green", high = "red")
      #scale_fill_gradientn(colours = colorspace::diverge_hcl(7))
      scale_fill_gradient2(high = "red", low = "blue", breaks = seq(-.2,.2,.05),
                           #labels = c("-.9","-.6","-.3","0",".3",".6",".9")
                           labels = paste(seq(-.2,.2,.05)))
  }
  if (annotate){
    p <- p + geom_text(color="black")
  }
  
  return(p)
}

getHeatMaps_side_by_side <- function(mat_list, title="", black_white = T, facet_wrap_rows = 1, indices = NULL,
                                     annotate_signs_mat = NULL, scale_fill_breaks = seq(-.2,.2,.05)){
  
  
  mat_list_melt <- lapply(
    seq_along(mat_list),
    function(mat_index){
      A <- mat_list[[mat_index]]
      ids = colnames(A)
      dimnames(A) <- NULL
      
      if(is.null(indices)){
        indices = 1:nrow(A)
      }
      
      A.melt <-
        reshape2::melt(A[indices,indices]) %>%
        mutate(
          from = ids[indices][Var1],
          to = ids[indices][Var2],
          Var1 = -1*Var1,
          training_cutoff = names(mat_list)[mat_index]
        )
    }
  ) %>%
    bind_rows() %>%
    mutate(training_cutoff = factor(training_cutoff, levels = names(mat_list)))
  
  if(!is.null(annotate_signs_mat)){
    ids = colnames(annotate_signs_mat)
    dimnames(annotate_signs_mat) <- NULL
    
    if(is.null(indices)){
      indices = 1:nrow(annotate_signs_mat)
    }
    
    Annotate.melt <-
      reshape2::melt(annotate_signs_mat[indices,indices]) %>%
      mutate(
        from = ids[indices][Var1],
        to = ids[indices][Var2],
        Var1 = -1*Var1
      ) %>%
      as_tibble() %>%
      rowwise() %>%
      mutate(labelx = c("-", "", "+")[sign(value) + 2]) %>%
      ungroup() %>%
      select(-value)
    
    mat_list_melt <- as_tibble(mat_list_melt) %>%
      left_join(Annotate.melt, by = c("Var1", "Var2"))
  } else {
    mat_list_melt <- mutate(mat_list_melt, labelx = "")
  }
  
  
  p <- mat_list_melt %>%
    ggplot(mapping = aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() +
    coord_fixed() +
    labs(col=NULL, title=title, x=NULL, y=NULL) +
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      #legend.position="none",
      plot.title=element_text(hjust=0.5),
      plot.subtitle=element_text(hjust=0.5),
      #axis.text = element_blank(),
      axis.ticks = element_blank(),
      #axis.title = element_blank(),
      panel.background = element_rect(fill = "lightblue",
                                      colour = "lightblue")
    ) +
    scale_x_continuous(breaks = 1:20) +
    scale_y_continuous(breaks = -1:-20, labels = paste(1:20)) +
    facet_wrap(~training_cutoff, nrow = facet_wrap_rows, labeller = label_parsed)
  
  if(black_white){
    p <- p + scale_fill_gradient(low = "white", high = "black")
  }else{
    p <- p + #scale_fill_gradient(low = "green", high = "red")
      #scale_fill_gradientn(colours = colorspace::diverge_hcl(7))
      scale_fill_gradient2(high = "red", low = "blue", breaks = scale_fill_breaks,
                           #labels = c("-.9","-.6","-.3","0",".3",".6",".9")
                           labels = paste(scale_fill_breaks))
  }
  
  #if(!is.null(annotate_signs_mat)){
  #  p <- p + geom_text(mapping = aes(label = labelx), color="black", size = .5)
  #}
  
  return(p)
}


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
