###### ----------
# Plant spectrum functions
# All code was done with ❤️ by two special Piisamirotte
##### ---------------

# source files
library(readr)
library(stringr)
source_rmd <- function(file_path) {
  stopifnot(is.character(file_path) && length(file_path) == 1)
  .tmpfile <- tempfile(fileext = ".R")
  .con <- file(.tmpfile) 
  on.exit(close(.con))
  full_rmd <- read_file(file_path)
  codes <- str_match_all(string = full_rmd, pattern = "```(?s)\\{r[^{}]*\\}\\s*\\n(.*?)```")
  stopifnot(length(codes) == 1 && ncol(codes[[1]]) == 2)
  codes <- paste(codes[[1]][, 2], collapse = "\n")
  writeLines(codes, .con)
  flush(.con)
  cat(sprintf("R code extracted to tempfile: %s\nSourcing tempfile...", .tmpfile))
  source(.tmpfile)
}

library(tidyverse)
library(vegan) # to use envfit
library(StatMatch)

# Missing traits
missing_traits <- function(plants) {
  select(plants, -species) %>% 
  is.na() %>% 
  colSums() %>% 
  as.data.frame() %>%
  transmute(perc_na = ./nrow(plants) * 100) %>% 
  rownames_to_column(var = "trait") %>% 
  mutate(trait = fct_reorder(trait, perc_na))
}


# Plotting missing traits
plot_missing_traits <- function(miss_traits) {
  ggplot(miss_traits, aes(perc_na, trait)) +
    geom_col() +
    xlim(0, 100) +
    labs(x="Percentage of missing data", y="Trait")
}

# Scatter plot with linear regression
plot_scatter_lin_reg <- function(data, mapping){
  ggplot(data, mapping) +
    geom_smooth(formula = y ~ x, method = "lm",
                se = FALSE, colour = "tomato4", size=.5) +
    geom_point(size = .7, alpha = .7)
}
# Creating pair plot with ggplot
pair_plots <- function(dat_pca) {
  ggpairs(dat_pca,
          diag = list(continuous = wrap("barDiag", bins=18)), # n bins manually selected for current dataset
          lower=list(continuous = plot_scatter_lin_reg),
          progress=FALSE) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}


# Compares each column with each other to see the pairwise amount of data available
pairwise_data_availability <- function(dat, fraction=T, perc=T){
  # converting to presence of data (True/False)
  dat_pres <- dat %>% 
    transmute(across(-species, negate(is.na)))
  ncol <- length(dat_pres)
  amount_data <- map_dfr(seq_along(dat_pres), function(i){
    # each column is compared with every other column only with the diagonal
    row <- rep(NA, ncol) # empty column
    for (j in 1:i){
      present <- sum(dat_pres[i] & dat_pres[j])
      if (fraction){
        present <- (present / nrow(dat_pres) ) %>% 
          round(2)
      }
      row[j] <- present
    }
    
    names(row) <- names(dat_pres)
    if (perc) {row <- row * 100} # convert to percentage
    row
  }
  )
  amount_data %>%  
    mutate(trait = names(dat_pres)) %>% 
    relocate(trait)
}


# Plot pairwise data availability
plot_pariwise_data_availability <- function(data_avail) {
  trait_order <- names(data_avail)
  data_avail %>%   
    gather("trait2", "data", -trait) %>% 
    mutate(
      trait_1 = factor(trait, rev(trait_order)),
      trait_2 = factor(trait2, trait_order)
    ) %>% 
    drop_na() %>% 
    ggplot(aes(trait_1, trait_2, label=data, fill=data)) +
    geom_label() +
    scale_fill_gradient(low="#cc2d17", high="#ffff3e") +
    theme(axis.text.x = element_text(angle=30, hjust = 1)) +
    labs(subtitle = "Percentage of data available for each trait combination")
}


# PCA function to look like a PCoA using euclidian distance
# This is the Plant Spectrum (ps) PCA (can find better name here)
ps_pca <- function(plants) {
  pcoa_data <- plants %>% select(-species)
  dist <- pcoa_data %>% 
    as.data.frame() %>% 
    dist()
  pcoa <- cmdscale(dist, eig=TRUE)
  return(pcoa)
}


# PCoA function
pcoa_gower <- function(plants) {
 pcoa_data <- plants %>% select(-species)
  dist <- pcoa_data %>% 
    as.data.frame() %>% 
    gower.dist()
  pcoa <- cmdscale(dist, eig=TRUE)
  return(pcoa)
}



# PCoA var
get_pcoa_var <- function(pcoa) {
  round(pcoa$eig/sum(pcoa$eig)*100, 1)
}


# Plot PCoA (without categorical traits)
plot_pcoa <- function(pcoa, plants, label="PCoA") {
  plants <- plants %>% 
    select(-species) %>% 
    rename_with(function(col) {paste0(col, "__")}, where(is.character))
  envfit <- envfit(pcoa, plants)
  
  pcoa_var <- get_pcoa_var(pcoa)
  plants <- plants %>% 
    mutate(
      pcoa_ax1=pcoa$points[,1],
      pcoa_ax2=pcoa$points[,2]
    )
  
  max_ax1 <- range(pcoa$points[,1]) %>% 
    abs() %>% 
    max()
  max_ax2 <- range(pcoa$points[,2]) %>% 
    abs() %>% 
    max()
  scaling_factor <- max(max_ax1, max_ax2)
  
  ggplot() +
    add_vectors(envfit, arrow_scaling = scaling_factor) +
    geom_point(aes(pcoa_ax1, pcoa_ax2), alpha=.3, data=plants) +
    xlab(paste(label, " ax 1 - ", pcoa_var [1], "%", sep="")) +
    ylab(paste(label, " ax 2 - ", pcoa_var [2], "%", sep="")) 
}


# Plot PCoA (categorical traits)
plot_pcoa_cat <- function(pcoa, plants) {
  plants <- plants %>% 
    select(-species) %>% 
    rename_with(function(col) {paste0(col, "__")}, where(is.character))
  envfit <- envfit(pcoa, plants)
  
  pcoa_var <- get_pcoa_var(pcoa)
  plants <- plants %>% 
    mutate(
      pcoa_ax1=pcoa$points[,1],
      pcoa_ax2=pcoa$points[,2]
    )
  
  # to use the range points to estimate the best scale of the arrows
  max_ax1 <- range(pcoa$points[,1]) %>% 
    abs() %>% 
    max()
  max_ax2 <- range(pcoa$points[,2]) %>% 
    abs() %>% 
    max()
  scaling_factor <- max(max_ax1, max_ax2)
  
  
  ggplot() +
    add_vectors(envfit, arrow_scaling = scaling_factor) +
    geom_point(aes(pcoa_ax1, pcoa_ax2), alpha=.5, data=plants) +
    xlab(paste("PCoA ax 1 - ", pcoa_var [1], "%", sep="")) +
    ylab(paste("PCoA ax 2 - ", pcoa_var [2], "%", sep="")) +
    add_factors(envfit)
}


# Add continuous variable vectors
# arrow_scaling: to make sure that in the plot the arrows looks similar to the points range
add_vectors <- function(envfit, arrow_scaling = 1) { 
  vectors <- 
    envfit$vectors$arrows %>%
    as.data.frame() %>% 
    rownames_to_column("trait") %>% 
    mutate(
      Dim1 = Dim1 * arrow_scaling,
      Dim2 = Dim2 * arrow_scaling
    )

  add_origins <- function(dat){
    tibble(
      trait = dat$trait,
      Dim1 = 0,
      Dim2 = 0
    ) %>% 
      bind_rows(dat)
  }

  vectors_pos <- vectors %>% 
    filter(Dim1 > 0) %>% 
    add_origins()
  vectors_neg <- vectors %>% 
    filter(Dim1 < 0) %>% 
    add_origins()

  list(
    geom_line(aes(x= Dim1, y = Dim2, group=trait), data = vectors_pos, arrow=arrow()),
    geom_line(aes(x= Dim1, y = Dim2, group=trait), data = vectors_neg, arrow=arrow(ends = "first")),
    geom_label_repel(aes(x=Dim1, y=Dim2, label=trait), data=vectors)
  )
}


# Plot PCoA var
plot_pcoa_var <- function(pcoa_var) {
  pcoa_var <- head(pcoa_var, 5) 
  pcoa_var %>% 
    tibble(
      axis = paste("axis", 1:5),
      percentage = pcoa_var
    ) %>% 
  ggplot() +
    geom_col(aes(axis, percentage), fill = "tomato4")
}


# add categorical factors
add_factors <- function(envfit) {
  factors <- envfit$factors$centroids %>% 
    as.data.frame() %>% 
    rownames_to_column("trait_value")
  factors <- factors %>% 
    mutate(
      trait = str_extract(trait_value, ".*(?=__)"), # everything before __
      value = str_extract(trait_value, "(?<=__).*")
    )
  
  list(
    geom_point(aes(Dim1, Dim2, col=trait, shape=trait), data = factors, size=4),
      geom_text_repel(aes(Dim1, Dim2, label=value), data = factors, nudge_y = 0.02)
  )
}