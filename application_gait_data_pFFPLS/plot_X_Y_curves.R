library(tidyverse)
library(gridExtra)



plot_X_Y_curves <- function(X, Y, argvals_X, argvals_Y, out_folder_ref_plot, num_obs_to_plot = NULL) {
  
  colnames(X) <- argvals_X
  colnames(Y) <- argvals_Y
  
  if ( is.null(num_obs_to_plot) ) {
    num_obs_to_plot <- nrow(X)
  }
  
  obs_to_plot <- 1:num_obs_to_plot
  
  
  
  dfY <-Y[obs_to_plot, ] %>% as.data.frame() %>% 
    mutate(id = obs_to_plot) %>% 
    pivot_longer(!id, names_to = "q", values_to = "Yq", names_transform = list(q = as.double))
  
  dfX <- X[obs_to_plot, ] %>% as.data.frame() %>% 
    mutate(id = obs_to_plot) %>% 
    pivot_longer(!id, names_to = "p", values_to = "Xp", names_transform = list(p = as.double))
  
  
  pY <- ggplot(dfY, aes(x = q, y = Yq, group = id, color = as.factor(id)) )  +
    geom_line(alpha = 0.6) + 
    ylab("Y(q)") +
    xlab("q") +
    theme_bw() +
    theme(text = element_text(size = 20), legend.position =  "none")
  
  
  pX <- ggplot(dfX, aes(x = p, y = Xp, group = id, color = as.factor(id)) ) +
    geom_line(alpha = 0.6)+ 
    ylab("X(p)") +
    xlab("p") +
    theme_bw() +
    theme(text = element_text(size = 20), legend.position =  "none")
  
  
  p_both <- grid.arrange(pX, pY, nrow = 1)
  
  
  ggsave(p_both,
         filename = paste0(out_folder_ref_plot,
                           paste0("X_Y_plot.png")  ),
         width = 12, height = 6 )
  
  ggsave(p_both,
         filename = paste0(out_folder_ref_plot,
                           paste0("X_Y_plot.pdf")  ),
         width = 12, height = 6 )
  
  
  ggsave(pX,
         filename = paste0(out_folder_ref_plot,
                           paste0("X_plot.png")  ),
         width = 6, height = 6 )
  
  ggsave(pX,
         filename = paste0(out_folder_ref_plot,
                           paste0("X_plot.pdf")  ),
         width = 6, height = 6 )
  
  
  ggsave(pY,
         filename = paste0(out_folder_ref_plot,
                           paste0("Y_plot.png")  ),
         width = 6, height = 6 )
  
  ggsave(pY,
         filename = paste0(out_folder_ref_plot,
                           paste0("Y_plot.pdf")  ),
         width = 6, height = 6 )
  
}
