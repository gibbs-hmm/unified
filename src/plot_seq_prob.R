#' Title plot_seq_prob - plot output of unified showing nucleotide and change point probability at each position.
#'
#' @param filename - _info_det file produced by unified
#'
#' @return a data frame of the contents of file name
#'          columns:
#'          Sequence, Sample_Count, A_prob', C_prob, G_prob, T_prob, 
#'          Position, Samples, Change_Point_probability 
#'          
plot_seq_prob <- function(filename) {
  require(tidyverse)
  require(ggpubr)
  
  df <- read.table(filename, 
                   col.names = c("Sequence", 'Sample_Count', 'A_prob', 'C_prob', 'G_prob', 'T_prob', 
                                 'Position', 'Samples', 'Change_Point_probability'))
  
  
  for(seq in unique(df$Sequence)) {
    df2 <- df %>% 
      filter(Sequence == seq)
    
    p1 <- ggplot(df2) +
      geom_line(aes(x = Position, y = A_prob, color = 'A')) +
      geom_line(aes(x = Position, y = C_prob, color = 'C')) +
      geom_line(aes(x = Position, y = G_prob, color = 'G')) +
      geom_line(aes(x = Position, y = T_prob, color = 'T')) +
      scale_color_manual(values = c('red', 'seagreen', 'blue', 'black'),
                         labels = c('A', 'C', 'G', 'T')) + 
      labs(y = 'Probability',
           color = 'Nucleotide') +
      ylim(c(0, 1))
  
    p2 <- ggplot(df2) +
      geom_linerange(aes(x = Position, ymin = 0, ymax = Change_Point_probability)) +
      ylab('Change Point Probability') +
      ylim(c(0,1))
    
    p3 <- ggarrange(plotlist = list(p1, p2),
                   ncol = 1,
                   nrow = 2,
                   common.legend = TRUE,
                   legend="bottom")
    
    p <- annotate_figure(p3, 
                         top = text_grob(paste('Sequence', seq), 
                                    color = "red",
                                    face = "bold", 
                                    size = 14))
    
    plot(p)
  }
  
  return(df)
}