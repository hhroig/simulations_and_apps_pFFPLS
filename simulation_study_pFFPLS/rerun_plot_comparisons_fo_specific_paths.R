list_of_paths = 
  str_c("results_simulations/",
        c(
          "set1_rep30_pen100_K7L7",  
          "set2_rep30_pen100_K40L40",   
          "set3_rep30_pen100_K40L40",
          "set1e_rep30_pen100_K7L7", 
          "set2e_rep30_pen100_K40L40",  
          "set3e_rep30_pen100_K40L40"
        ), 
        "/")



for (pathito in list_of_paths) {
  
  compare_methods_fun(input_folder = pathito, 
                      zoom_r2_lower = 0.5, 
                      do_rough_r2 = TRUE)
}