library("smotefamily")
library("readr")
library("dplyr")
library("magrittr")
library("data.table")
library("caTools")
library("tidytable")
library("survival")

##############
# Trae datos #
##############

data <- read.csv("support2.csv") %>% as.data.table(.)

drop_col_1 <- c("sfdm2","ca","income","dzclass","edu","totmcst","prg2m","prg6m","pafi","alb","bili","ph","glucose","bun","urine","adlp","adls")
data_1 <- data %>% 
  filter(dzclass == "Cancer" , race != "") %>% 
  select(-all_of(drop_col_1)) %>% 
  na.omit() %>% 
  get_dummies() %>% 
  select(-all_of(c("sex","dzgroup","race","dnr")))

# load(url('https://hbiostat.org/data/repo/support.sav'))
# View(support)

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
  na_if("") %>% 
  na.omit() %>% 
  get_dummies() %>% 
  select(where(is.numeric))


# ---------------------------
nuestro_stepwise <- function(tiempo, estado, datos_sin_t_ni_status){
  
  maximos <- c()
  minimos <- c()
  
  variables <- names(datos_sin_t_ni_status)
  m <- length(variables)
  aux_c.index <- c()
  surv_obj <- Surv(tiempo, estado)
  for (id_variable in 1:m){
    surv_fit <- coxph(surv_obj ~ ., data = select(datos_sin_t_ni_status, variables[id_variable]))
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
      datitos <- select(datos_sin_t_ni_status, c(variables[id_variable], all_of(maximos)))
      surv_fit <- coxph(surv_obj ~ ., data = datitos)
      c.index <- as.vector(surv_fit$concordance["concordance"])
      aux_c.index[id_variable] <- c.index
      names(aux_c.index)[id_variable] <- variables[id_variable]
    }
    min_c.index <- names(aux_c.index)[which.min(aux_c.index)]
    max_c.index <- names(aux_c.index)[which.max(aux_c.index)]
    minimos <- c(minimos,min_c.index)
    maximos <- c(maximos,max_c.index)
  }
  out_datitos  <- select(datos_sin_t_ni_status, all_of(maximos))
  out_surv_fit <- coxph(surv_obj ~ ., data = out_datitos)
  
  return(out_surv_fit)
}
# fitness <- nuestro_stepwise(data_4$time, data_4$status, data_4[,-c(12,13)])
