library(tidyverse)
library(gridExtra)
library(viridis)
library(ggpubr)
library(scales)
library(plotly)
library(writexl)



beta_num_to_text <- function(beta_num_in) {
  beta_txt <- "est"
  return(beta_txt)
  
}




input_folder = "results_bike/reps100_pen400_K8L8/"



out_folder <- paste0(input_folder, "results_plots/")

if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}


# show_col(hue_pal()(6))

color_codes <- c(
  "pFFPLS" = hue_pal()(6)[1],
  "pFFR_I" = hue_pal()(6)[2],
  "pFFR_RS" = hue_pal()(6)[3],
  "FFPLS_OB" = hue_pal()(6)[4],
  "FFPLS" = hue_pal()(6)[5]
)

## Final Models Files  ---------------------------------------------------

# IMSE for Beta and validation Y:
all_final_res <- data.frame()

final_res_files <- list.files(path = input_folder, pattern = "final_models")

for (ind_file in final_res_files) {
  
  all_final_res <- rbind(
    all_final_res,
    readRDS(paste0(input_folder, ind_file))
  )
  
}

# CVEs:
all_cves <- data.frame()

all_cves_files <- list.files(path = input_folder, pattern = "cves_rep")

for (ind_file in all_cves_files) {
  
  all_cves <- rbind(
    all_cves,
    readRDS(paste0(input_folder, ind_file))
  )
  
}

# Computation time
all_computation_times <- data.frame()

all_computation_times_files <- list.files(path = input_folder, pattern = "computation_times")

for (ind_file in all_computation_times_files) {
  
  all_computation_times <- rbind(
    all_computation_times,
    readRDS(paste0(input_folder, ind_file))
  )
  
}


# Best number bases for FFPLS:

all_best_num_bases_FFPLS <- data.frame()

all_best_num_bases_FFPLS_files <- list.files(path = input_folder, pattern = "best_num_bases_FFPLS")

for (ind_file in all_best_num_bases_FFPLS_files) {
  
  all_best_num_bases_FFPLS <- rbind(
    all_best_num_bases_FFPLS,
    readRDS(paste0(input_folder, ind_file))
  )
  
}



# Best number bases for FFPLS:

all_best_num_bases_FFPLS <- data.frame()

all_best_num_bases_FFPLS_files <- list.files(path = input_folder, pattern = "best_num_bases_FFPLS")

for (ind_file in all_best_num_bases_FFPLS_files) {
  
  all_best_num_bases_FFPLS <- rbind(
    all_best_num_bases_FFPLS,
    readRDS(paste0(input_folder, ind_file))
  )
  
}


# Best Lambdas files:

all_best_lambdas <- data.frame()

all_best_lambdas_files <- list.files(path = input_folder, pattern = "best_lambdas")

for (ind_file in all_best_lambdas_files) {
  
  all_best_lambdas <- rbind(
    all_best_lambdas,
    readRDS(paste0(input_folder, ind_file))
  )
  
}

# All R^2 files:

all_best_r2 <- data.frame()

all_best_r2_files <- list.files(path = input_folder, pattern = "R2s_rep")

for (ind_file in all_best_r2_files) {
  
  all_best_r2 <- rbind(
    all_best_r2,
    readRDS(paste0(input_folder, ind_file))
  )
  
}

# Betas files:

all_betas <- data.frame()

all_betas_files <- list.files(path = input_folder, pattern = "betas_")

for (ind_file in all_betas_files) {
  
  all_betas <- rbind(
    all_betas,
    readRDS(paste0(input_folder, ind_file))
  )
  
}


# Limit the number of components to plot ----------------------------------


all_final_res <- all_final_res %>%
  filter(nComp <= cap_ncomp) %>% 
  mutate(nComp = as.factor(nComp))

all_cves <- all_cves %>%
  filter(nComp <= cap_ncomp) %>% 
  mutate(nComp = as.factor(nComp))

all_best_lambdas <- all_best_lambdas %>%
  filter(nComp <= cap_ncomp) %>% 
  mutate(nComp = as.factor(nComp)) 

all_computation_times <- all_computation_times %>%
  filter(nComp <= cap_ncomp) %>% 
  mutate(nComp = as.factor(nComp)) 

all_betas <- all_betas %>%
  filter(nComp <= cap_ncomp) %>% 
  mutate(nComp = as.factor(nComp))


all_best_num_bases_FFPLS <- all_best_num_bases_FFPLS %>%
  filter(nComp <= cap_ncomp) %>% 
  mutate(nComp = as.factor(nComp)) 


# Limits:

# cve_limits <- range(all_cves$CVE)
# imse_limits <- range(all_final_res$imse )
# mean_imse_Y_val_limits <- range(all_final_res$mean_imse_Y_val )



# IMSE + MSE --------------------------------------------------------------


out_folder_IMSE_CVEs <- paste0(out_folder, "IMSEs_CVEs_specific/")

if (!dir.exists(out_folder_IMSE_CVEs)) {
  dir.create(out_folder_IMSE_CVEs)
}




# Plot 1:
# IMSE predicción sobre la test para todos los métodos y con FFPLS (con 3 cp), 
# FFPLS_OB (con 5 cp) y pFFPLS (5 cp). Es decir, tal y como está, pero considerando 
# 3 componentes en el caso del FFPLS.


all_final_res_ffpls3 <- all_final_res %>% filter(method == "FFPLS", nComp == 3)

all_final_res_ncomp5_no_ffpls <-  all_final_res %>% 
  filter(method != "FFPLS") %>% 
  filter(nComp == 5)

all_final_res_1 <- rbind(all_final_res_ncomp5_no_ffpls, all_final_res_ffpls3)



p_imse_val_1 <- ggplot(all_final_res_1 ,
                     aes(x = method, y = mean_imse_Y_val, fill = method)) +
  geom_boxplot(position=position_dodge(0.8)) +
  ylab( expression( "Test: IMSE(Y)" )  ) +
  xlab("") +
  scale_fill_manual(values = color_codes) +
  theme_bw() +
  theme(legend.position="none", text = element_text(size = 20)) +
  labs(fill = "")

p_imse_val_1

ggsave(p_imse_val_1,
       filename = paste0(out_folder_IMSE_CVEs,
                         "test_imseY_fig1.png"  ),
       width = 8, height = 8 )

ggsave(p_imse_val_1,
       filename = paste0(out_folder_IMSE_CVEs,
                         "test_imseY_fig1.pdf" ),
       width = 8, height = 8 )



# Plot 2:
# IMSE predicción sobre la test para todos los métodos, excepto FFPLS_OB, y con 
# FFPLS (con 3 cp) y pFFPLS (5 cp).


all_final_res_2 <- all_final_res_1 %>% 
  filter(method != "FFPLS_OB") 



p_imse_val_2 <- ggplot(all_final_res_2 ,
                     aes(x = method, y = mean_imse_Y_val, fill = method)) +
  geom_boxplot(position=position_dodge(0.8)) +
  ylab( expression( "Test: IMSE(Y)" )  ) +
  xlab("") +
  scale_fill_manual(values = color_codes) +
  theme_bw() +
  theme(legend.position="none", text = element_text(size = 20)) +
  labs(fill = "")

p_imse_val_2

ggsave(p_imse_val_2,
       filename = paste0(out_folder_IMSE_CVEs,
                         "test_imseY_fig2.png"  ),
       width = 8, height = 8 )

ggsave(p_imse_val_2,
       filename = paste0(out_folder_IMSE_CVEs,
                         "test_imseY_fig2.pdf" ),
       width = 8, height = 8 )


# Plot 3
# CVE estimación sobre la training solo para los métodos FFPLS y pFFPLS y para 
# un número de componentes variando de 1 a 5.

all_cves_3 <- all_cves %>% 
  filter((method == "FFPLS") | (method == "pFFPLS"),  
         nComp <= 5)



p_cve_3 <- ggplot(all_cves_3,
                aes(x = nComp, y = CVE, fill = method)) +
  geom_boxplot(position=position_dodge(0.8))  +
  ylab("Training: CVE(Y)") +
  xlab("# of components") +
  scale_fill_manual(values = color_codes)+
  theme_bw()  +
  theme(legend.position="bottom", text = element_text(size = 20)) +
  labs(fill = "")


ggsave(p_cve_3,
       filename = paste0(out_folder_IMSE_CVEs,
                         "train_cveY_fig3.png"  ) ,
       width = 8, height = 8 )

ggsave(p_cve_3,
       filename = paste0(out_folder_IMSE_CVEs,
                         "train_cveY_fig3.pdf"  ) ,
       width = 8, height = 8 )