# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("sems1_ker.R")
sa_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric

set.seed(sa_id)
sim(sa_id) %>%
  write.csv(file = paste("res/res", sa_id, ".csv", sep = ""))




  


