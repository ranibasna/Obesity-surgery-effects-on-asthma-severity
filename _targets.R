# _targets.R
library(targets)
library(tarchetypes)
source("Rfunctions/Functions.R")
tar_option_set(packages = c("tidyverse","tibble","lubridate","mltools","stringi","data.table"))
list(
  #record state of R environment
  tar_target(lockenv,renv::snapshot()),
  # set the data path
  tar_target(raw_data_path,"Data/seconddata/Final_data_v4.csv",format = "file"),
  # read the raw dta data
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
  # get the Emergency data
  tar_target(fixed_Hos_data, get_Hos_vars(raw_data_date)),
  # get num data
  tar_target(num_data, raw_data_date %>% select(where(is.numeric)) %>% select(-c(prednisolone, packages)) %>% distinct()),
  # get charechter data
  tar_target(char_data, raw_data_date %>% select(where(is.character) | id) %>%  distinct() %>% select(-c(starts_with("smok"))) %>% mutate(across(everything(), ~ na_if(.,"")))),
  # bind data
  tar_target(ready_data, fixed_data_OSC[[1]] %>%  inner_join(fixed_burst_data, by=c("id", "date_surgery")) %>% 
               inner_join(fixed_em_data, by=c("id", "date_surgery")) %>% 
               inner_join(fixed_Hos_data, by=c("id", "date_surgery")) %>% 
               inner_join(fixed_smko_data, by="id") %>% 
               inner_join(char_data, by="id") %>% 
               inner_join(num_data, by="id") %>% 
               inner_join(fixed_data_oh, by="id")),
  # save the data
  tar_target(ready_Data, write.csv(x = ready_data, file =  "/home/rstudio/ObesityAsthma/Intermediate/ProcessedData/ready_data_tar.csv", row.names = FALSE))
)
  

