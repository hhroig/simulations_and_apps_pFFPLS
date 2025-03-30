source("compare_methods_fofr_with_ivanescu_ramsay_silverman.R")



list_subfolder <- c(
  "results/set1_rep60_pen100_K7L7/",
  "results/set1e_rep60_pen100_K7L7/",
  "results/set3_rep60_pen100_K40L40/",
  "results/set3e_rep60_pen100_K40L40/"
)



for (in_folds in list_subfolder) {
  print(in_folds)
  
  compare_methods_fun(in_folds, 
                      zoom_r2_lower = 0, 
                      do_rough_r2 = F,
                      theta = 30,   # Angle for viewing (rotation beta surface)
                      phi = 30    # Angle for viewing (tilt beta surface)
  )
  
}


# input_folder = "results_simulations/set3e_rep60_pen100_K40L40_/"
# 
# 
# source("compare_methods_fofr_with_ivanescu_ramsay_silverman.R")
# 
# compare_methods_fun(input_folder, 
#                                 zoom_r2_lower = 0, 
#                                 do_rough_r2 = False,
#                                 theta = 30,   # Angle for viewing (rotation beta surface)
#                                 phi = 30    # Angle for viewing (tilt beta surface)
# )
# 
# 
# # 
# # # -----------------------------------------------------------------------------
# # 
# # 
# zoom_r2_lower = 0
# do_rough_r2 = FALSE
# theta = 30   # Angle for viewing (rotation beta surface)
# phi = 30    # Angle for viewing (tilt beta surface)
