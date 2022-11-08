library("smotefamily")
library("readr")
library("dplyr")
library("magrittr")
library("data.table")
library("caTools")
library("tidytable")
library("survival")
library("labelled")

##############
# Trae datos #
##############

# load(url('https://hbiostat.org/data/repo/support.sav'))
# # View(support)
# m <- length(names(support))
# V = c()
# Label = c()
# 
# for (v in 1:m){
#   V[v] = names(support)[v]
#   Label[v] = var_label(support)[names(support)[v]]
# }
# df <- data.frame(cbind(V,Label))

data <- read.csv("support2.csv") %>% as.data.table(.)

drop_col_1 <- c("sfdm2","ca","income","dzclass","edu","totmcst","prg2m","prg6m","pafi","alb","bili","ph","glucose","bun","urine","adlp","adls")
data_1 <- data %>% 
  filter(dzclass == "Cancer" , race != "") %>% 
  select(-all_of(drop_col_1)) %>% 
  na.omit() %>% 
  get_dummies() %>% 
  select(-all_of(c("sex","dzgroup","race","dnr")))



data_2 <- data %>%
  filter(dzclass == "ARF/MOSF" , race != "") %>% 
  select(-all_of(drop_col_1)) %>% 
  na.omit() %>% 
  get_dummies() %>% 
  select(-all_of(c("sex","dzgroup","race","dnr")))

data_3 <- data %>% 
  filter(dzclass == "COPD/CHF/Cirrhosis") %>% 
  select(-all_of(drop_col_1)) %>% 
  na.omit() %>% 
  get_dummies() %>% 
  select(-all_of(c("sex","dzgroup","race","dnr")))

data <- read.csv("METABRIC_RNA_Mutation.csv") %>% as.data.table(.)

data_4 <- data %>% 
  select(1:31) %>% 
  filter(death_from_cancer != "Died of Other Causes") %>% 
  mutate(
    time = overall_survival_months,
    status = case_when(overall_survival == 1 ~ 0, TRUE ~ 1)
    ) %>% 
  select(-all_of(c("overall_survival_months",
                   "overall_survival",
                   "patient_id","cancer_type",
                   "death_from_cancer",
                   "pam50_._claudin.low_subtype"))) %>% 
  dplyr::na_if("") %>% 
  na.omit() %>% 
  get_dummies() %>% 
  select(where(is.numeric))


nuestro_stepwise <- function(tiempo, estado, covariables){
  
  maximos <- c()
  minimos <- c()
  
  variables <- names(covariables)
  m <- length(variables)
  aux_c.index <- c()
  surv_obj <- Surv(tiempo, estado)
  for (id_variable in 1:m){
    surv_fit <- coxph(surv_obj ~ ., data = select(covariables, variables[id_variable]))
    c.index <- as.vector(surv_fit$concordance["concordance"])
    aux_c.index[id_variable] <- c.index
    names(aux_c.index)[id_variable] <- variables[id_variable]
  }
  min_c.index <- names(aux_c.index)[which.min(aux_c.index)]
  max_c.index <- names(aux_c.index)[which.max(aux_c.index)]
  
  minimos[1] <- min_c.index
  maximos[1] <- max_c.index
  
  while (m > 1){
    variables <- variables[!variables %in% c(minimos,maximos)]
    m <- length(variables) - 2
    aux_c.index <- c()
    for (id_variable in 1:m){
      covariables_aux <- select(covariables, c(variables[id_variable], all_of(maximos)))
      surv_fit <- coxph(surv_obj ~ ., data = covariables_aux)
      c.index <- as.vector(surv_fit$concordance["concordance"])
      aux_c.index[id_variable] <- c.index
      names(aux_c.index)[id_variable] <- variables[id_variable]
    }
    min_c.index <- names(aux_c.index)[which.min(aux_c.index)]
    max_c.index <- names(aux_c.index)[which.max(aux_c.index)]
    minimos <- c(minimos,min_c.index)
    maximos <- c(maximos,max_c.index)
  }
  out_covariables_aux  <- select(covariables, all_of(maximos))
  out_surv_fit <- coxph(surv_obj ~ ., data = out_covariables_aux)
  
  return(out_surv_fit)
}

fit_1 <- nuestro_stepwise(data_1$d.time, data_1$death, data_1[,-c(2,5)])
fit_2 <- nuestro_stepwise(data_2$d.time, data_2$death, data_2[,-c(2,5)])
fit_3 <- nuestro_stepwise(data_3$d.time, data_3$death, data_3[,-c(2,5)])

fit_4 <- nuestro_stepwise(data_4$time, data_4$status, data_4[,-c(12,13)])
