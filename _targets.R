# _targets.R
library(targets)
library(tarchetypes)
source("Rfunctions/Functions.R")
tar_option_set(packages = c("tidyverse","tibble","lubridate","mltools","stringi","data.table","mice","naniar","finalfit", "tidymodels","bnlearn","embed","visNetwork", "cowplot","ggplot2","brms","tidybayes","bayesplot"))
list(
  #record state of R environment
  tar_target(lockenv,renv::snapshot()),
  # set the data path
  tar_target(raw_data_path,"Data/seconddata/Final_data_v4.csv",format = "file"),
  # read the raw data
  tar_target(raw_data,read.csv(raw_data_path)),
  # Process the data
  ## Fix the date variables
  tar_target(raw_data_date, raw_data %>% mutate(across(starts_with(c("study", "date", "death")) , dmy))),
  ## get the OH for the drug data
  tar_target(fixed_data_oh, get_oh_drug(raw_data_date)),
  ## get the smoking data
  tar_target(fixed_smko_data, get_smok_data(raw_data_date)),
  # get the OSC vars
  tar_target(fixed_data_OSC, get_OSC(raw_data_date)),
  # get_burst_var
  tar_target(fixed_burst_data, get_burst_var(raw_data_date, fixed_data_OSC[[2]])),
  # get the Emergency data
  tar_target(fixed_em_data, get_Em_vars(raw_data_date)),
  # get the Hospitalization data
  tar_target(fixed_Hos_data, get_Hos_vars(raw_data_date)),
  # get num data
  tar_target(num_data, raw_data_date %>% select(where(is.numeric)) %>% select(-c(prednisolone, packages)) %>% distinct()),
  # get charechter data
  tar_target(char_data, raw_data_date %>% select(where(is.character) | id) %>%  distinct() %>% select(-c(starts_with(c("smok","drug")))) %>% mutate(across(everything(), ~ na_if(.,"")))),
  # bind data
  tar_target(ready_data, fixed_data_OSC[[1]] %>%  inner_join(fixed_burst_data, by=c("id")) %>% 
               inner_join(fixed_em_data, by=c("id", "date_surgery")) %>% 
               inner_join(fixed_Hos_data, by=c("id", "date_surgery")) %>% 
               inner_join(fixed_smko_data, by="id") %>% 
               inner_join(char_data, by="id") %>% 
               inner_join(num_data, by="id") %>% 
               inner_join(fixed_data_oh, by="id")),
  # save the data
  tar_target(ready_Data, write.csv(x = ready_data, file =  "/home/rstudio/ObesityAsthma/Intermediate/ProcessedData/ready_data_tar.csv", row.names = FALSE)),
  # save encrypted data for shiny
  # mydata_ready %>% select(-c( bleeding_amount, sleep_apnea, b_hba1c)) %>% encrypt(id)
  # Prepare the data for the BN analysis
  # tar_target(BN_processed_data, get_BN_processed_data(ready_data)),
  # tar_target(BN_processed_data_inc, get_BN_processed_data(ready_data, Not_inclusion = FALSE)),
  tar_target(BN_processed_data_inc_acc, get_BN_processed_data(BN_ready_data = ready_data, Is_acc_main = TRUE, Is_change = FALSE, Not_inclusion = FALSE)),
  # Define the bl and wl
  # tar_target(wl, matrix(c("wght_loss6w", "hos_Diff", "surgery", "wght_loss6w"), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))),
  # tar_target(wl_no_wl, matrix(c( "surgery", "wght_loss6w"), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))),
  tar_target(wl_bmi,matrix(c("surgery", "BMI_6w"), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))),  
  # tar_target(bl, get_bl(final_BN_data = BN_processed_data, bmi_var = TRUE)), # since we replaced the weight lost by bmi variable
  tar_target(bl_inc, get_bl(final_BN_data = BN_processed_data_inc_acc, bmi_var = TRUE)),
  # learn the dag and Aggregate
  # tar_target(BN_bootstrap, boot.strength(BN_processed_data, R = 100, algorithm = "tabu", algorithm.args = list(blacklist = bl, whitelist = wl))),
  # tar_target(BN_avg_dag, averaged.network(BN_bootstrap)),
  # tar_target(BN_bootstrap_inc, boot.strength(BN_processed_data_inc, R = 100, algorithm = "tabu", algorithm.args = list(blacklist = bl_inc, whitelist = wl))),
  # tar_target(BN_avg_dag_inc, averaged.network(BN_bootstrap_inc)),
  # tar_target(BN_bootstrap_inc_no_wl, boot.strength(BN_processed_data_inc, R = 100, algorithm = "tabu", algorithm.args = list(blacklist = bl_inc, whitelist = wl_no_wl))),
  # tar_target(BN_avg_dag_inc_no_wl, averaged.network(BN_bootstrap_inc_no_wl)),
  tar_target(BN_bootstrap_inc_acc, boot.strength(BN_processed_data_inc_acc, R = 300, algorithm = "tabu", algorithm.args = list(blacklist = bl_inc, whitelist = wl_bmi, score = "bde", iss = 10))),
  tar_target(BN_avg_dag_inc_acc, averaged.network(BN_bootstrap_inc_acc)),
  # Fit the BN
  # tar_target(BN_fitted_obj, bn.fit(cextend(BN_avg_dag), BN_processed_data)),
  # tar_target(BN_fitted_obj_inc, bn.fit(cextend(BN_avg_dag_inc), BN_processed_data_inc)),
  # tar_target(BN_fitted_obj_inc_no_wl, bn.fit(cextend(BN_avg_dag_inc_no_wl), BN_processed_data_inc)),
  tar_target(BN_fitted_obj_inc_acc, bn.fit(cextend(BN_avg_dag_inc_acc), BN_processed_data_inc_acc)),
  # impute the data
  tar_target(Imputed_data,  mice(data = ready_data, m = 10, maxit = 10)),
  # render the report for the missing data
  tar_render(report, "Intermediate/Reports/miss_data_validation_report.Rmd"),
  # render the report for the BN analysis
  # tar_render(report_BN, "Intermediate/Reports/BN_Analysis_Obesity.Rmd"),
  # 
  # tar_render(report_BN_par, "Intermediate/Reports/BN_Analysis_Parametrized.Rmd", params = list(BN_data = BN_processed_data_inc, BN_avg_dag_par = BN_avg_dag_inc, BN_fitted_par = BN_fitted_obj_inc )),
  # tar_render(report_BN_par_inc_acc, "Intermediate/Reports/BN_Analysis_Parametrized.Rmd", params = list(BN_data = BN_processed_data_inc_acc, BN_avg_dag_par = BN_avg_dag_inc_acc, BN_fitted_par = BN_fitted_obj_inc_acc), output_file = "BN_analysis_inc_acc"),
  #tar_render(rmarkdown::render("MyDocument.Rmd", params = list(region = "Asia")))
  # Prepare the data for Interrupted time series analysis
  tar_target(IST_data_init, get_time_df(ready_data = raw_data)),
  tar_target(ITS_data, get_ITS_data(timeVis_data = IST_data_init)),
  # set the prior distributions
  tar_target(prior_dis, c(set_prior("normal(0,1)", class = "b", coef = "elapsed_time"),
                                          set_prior("normal(0,1)", class = "b", coef = "step_change"),
                                          set_prior("normal(0,1)", class = "b", coef = "timeAfterSurg"),
                                          set_prior("normal(5,5)", class = "Intercept")
                                          # set_prior("uniform(-1,1)", class = "ar")
  )),
  tar_target(modelPoi_main, brms::brm(
    NumberPeople ~ elapsed_time + step_change + timeAfterSurg + deathMonthAccCount + ar(p = 1),
    # ppc_outcome ~ elapsed_time + step_change + timeAfterSurg ,
    data = ITS_data,
    family = poisson(),
    prior = prior_dis,
    chains = 4,
    iter = 3000,
    warmup = 1000,
    sample_prior = TRUE
  )),
  # final report
  tar_render(report_BN_par_inc_acc, "Intermediate/Reports/BN_Analysis_Parametrized.Rmd", params = list(BN_data = BN_processed_data_inc_acc, BN_avg_dag_par = BN_avg_dag_inc_acc, BN_fitted_par = BN_fitted_obj_inc_acc, ITS_data_par = ITS_data, Brms_model_par = modelPoi_main), output_file = "BN_analysis_inc_acc")
)
  



