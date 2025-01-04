# R/plot_rhythmic_gene.R

#' Plot Rhythmic Expression for a Specific Gene Across Conditions
#'
#' This function plots the rhythmic expression curves for a specified gene under two conditions using the fitted parameters.
#'
#' @param fitted_params A data frame containing fitted parameters from `perform_analysis`. It should include columns: Gene, Condition, Alpha, BetaCos, BetaSin.
#' @param gene_id The ID of the gene to plot. Should match the 'Gene' column in `fitted_params`.
#' @param period The rhythmic period (default is 24).
#' @param time_range A numeric vector specifying the range of time points for plotting (default is seq(0, 24, by = 0.1)).
#' @param colors A named vector specifying colors for conditions, e.g., c("1" = "blue", "2" = "red").
#' @param ... Additional arguments passed to `ggplot`.
#'
#' @return A ggplot object displaying the rhythmic expression curves.
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' PlotRhythmicGene(
#'   fitted_params = result$FittedParams,
#'   gene_id       = 166,
#'   period        = 24
#' )

PlotRhythmicGene <- function(
    fitted_params,
    gene_id,
    period     = 24,
    time_range = seq(0, 24, by = 0.1),
    colors     = c("Condition1" = "blue", "Condition2" = "red"),
    ...
) {

  # Load ggplot2 package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required but not installed.")
  }
  library(ggplot2)

  # Subset the data for the specified gene
  gene_data <- subset(fitted_params, Gene == gene_id)

  if (nrow(gene_data) < 2) {
    stop("Error: This function assumes 2 conditions. Check your data or ncond.")
  }

  # Extract gene name
  gene_name <- unique(gene_data$GeneName)
  if (length(gene_name) != 1) {
    stop("Error: GeneName should be unique for the specified gene_id.")
  }

  # Extract data for Condition = 1 and Condition = 2
  cond1_data <- subset(gene_data, Condition == 1)
  cond2_data <- subset(gene_data, Condition == 2)

  if (nrow(cond1_data) == 0 || nrow(cond2_data) == 0) {
    stop("Error: Both conditions must be present for the specified gene.")
  }

  # Extract parameters (Alpha, BetaCos, BetaSin)
  alpha1  <- cond1_data$Alpha
  bcos1   <- cond1_data$BetaCos
  bsin1   <- cond1_data$BetaSin

  alpha2  <- cond2_data$Alpha
  bcos2   <- cond2_data$BetaCos
  bsin2   <- cond2_data$BetaSin

  # formula = alpha + bcos * cos(2*pi*t/period) + bsin * sin(2*pi*t/period)
  pred_cond1 <- alpha1 + bcos1 * cos(2 * pi * time_range / period) +
    bsin1 * sin(2 * pi * time_range / period)

  pred_cond2 <- alpha2 + bcos2 * cos(2 * pi * time_range / period) +
    bsin2 * sin(2 * pi * time_range / period)

  # Create a data frame for plotting
  plot_df <- data.frame(
    Time = rep(time_range, 2),
    Expression = c(pred_cond1, pred_cond2),
    Condition = factor(rep(c("Condition1", "Condition2"), each = length(time_range)),
                       levels = c("Condition1", "Condition2"))
  )

  # Generate the plot
  p <- ggplot(plot_df, aes(x = Time, y = Expression, color = Condition)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = colors) +
    labs(
      title = paste("Gene", gene_id, ":", gene_name),
      x = "Time (h)",
      y = "Expression (fitted)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )

  return(p)
}
