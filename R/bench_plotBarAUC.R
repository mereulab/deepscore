#' Benchmarking: Plotting the output of the function `bench_calcAUC()`
#'
#' This function plots the output of the function `bench_calcAUC()`
#'
#' @param out The output of the function `bench_calcAUC()`
#' @param bar_cols Color vector for the bar plot layer
#' @param ct_mapping Option to fill bar plot layer by grouping cols, if not provided, by default will only take information from `out`
#' @param group_case Column specifying the benchmarking cases
#' @param group_cols Color vector for the line and point layers
#' @param alpha_vec Alpha vector for the bar plot layer
#'
#' @import dplyr
#' @import reshape2
#' @import ggnewscale
#'
#' @return A dataset of AUC per cell type of each case
#' @export
#'

bench_plotBarAUC<-function(out,bar_cols=NULL,ct_mapping=NULL,group_case="Case",group_cols=NULL,alpha_vec=NULL){
  out_melt<-melt(out)
  out_melt[,group_case]<-factor(out_melt[,group_case],levels = unique(out_melt[,group_case]))
  out_melt$variable<-as.factor(out_melt$variable)

  if(missing(bar_cols)){
    nb<-length(unique(out_melt$variable))
    bar_cols<-palette(rainbow(nb))
    names(bar_cols)<-unique(out_melt$variable)
  }
  if(missing(group_cols)){
    ng<-length(unique(out_melt[,group_case]))
    group_cols<-palette(rainbow(ng))
    names(group_cols)<-unique(out_melt[,group_case])
  }
  if(missing(alpha_vec)){
    an<-0.99/length(unique(out_melt[,group_case]))
    alpha_vec<-seq(0.3, 0.99, by=an)
    names(alpha_vec)<-names(group_cols)
  }
  if(missing(ct_mapping)){
    fill_by<-"variable"
  }
  else{
    out_melt$ct <- ct_mapping[as.character(out_melt$variable)]
    fill_by<-"ct"
  }

  plot<-ggplot(data = out_melt) +
    # Bar plot layer
    geom_bar(
      stat = "identity",
      position = position_dodge(),
      aes(x = variable, y = value, fill = .data[[fill_by]], group = .data[[group_case]], alpha = .data[[group_case]]),
      col = "black"
    ) +
    scale_fill_manual(values = bar_cols) +
    scale_alpha_manual(values = alpha_vec) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(fill = guide_legend(title = "Cell Type")) +
    # Line plot layer
    geom_line(
      aes(x = variable, y = value, col = .data[[group_case]], group = .data[[group_case]]),
      linewidth = 1
    ) +
    scale_color_manual(values = group_cols) +guides(color=guide_legend(title = paste0(group_case)))+
    new_scale_color() +

    # Point layer
    geom_point(
      aes(x = variable, y = value, col = .data[[group_case]]),
      size = 3
    ) +
    scale_color_manual(values = group_cols) +
    guides(
      alpha = guide_legend(title = paste0(group_case)),
      color = guide_legend(title = paste0(group_case))
    ) +
    scale_size_continuous(guide = "none")

  return(plot)
}
