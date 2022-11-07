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

data_2 <- data %>% filter(dzclass == "ARF/MOSF")
data_3 <- data %>% filter(dzclass == "COPD/CHF/Cirrhosis")

data <- read.csv("METABRIC_RNA_Mutation.csv") %>% as.data.table(.)

data_4 <- data %>% 
  select(1:31) %>% 
  filter(death_from_cancer != "Died of Other Causes") %>% 
  mutate(status = case_when(death_from_cancer == "Died of Disease" ~ 1, TRUE ~ 0))

##############
# Train-Test #
##############

# sample <- sample.split(data_1_balanced$class, SplitRatio = 0.8)
# train  <- subset(data_1_balanced, sample == TRUE)
# test   <- subset(data_1_balanced, sample == FALSE)

# sample <- sample.split(data_1$death, SplitRatio = 0.8)
# train  <- subset(data_1, sample == TRUE)
# test   <- subset(data_1, sample == FALSE)
# train[,"dzgroup_Colon Cancer"] %>% table %>% prop.table
# test[,"dzgroup_Colon Cancer"]  %>% table %>% prop.table

fit <- coxph(Surv(data_1$d.time, data_1$death) ~ . , data = data_1) #Concordance= 0.765  (se = 0.007 )

pred.fit <- predict(fit)
concordance(Surv(data_1$d.time, data_1$death) ~ pred.fit, data=data_1)


# ---------------------------
tiempo <- data_1$d.time
estado <- data_1$death
datos_sin_t_ni_status <- data_1[,-c(2,5)]
nuestro_stepwise <- function(tiempo, estado, datos_sin_t_ni_status){
  variables <- names(datos_sin_t_ni_status)
  m <- length(variables)
  aux_c.index <- c()
  surv_obj <- Surv(tiempo, estado)
  for (id_variable in 1:m){
    surv_fit <- coxph(surv_obj ~ ., data = datos_sin_t_ni_status[,..id_variable])
    c.index <- as.vector(surv_fit$concordance["concordance"])
    aux_c.index[id_variable] <- c.index
    names(aux_c.index)[id_variable] <- variables[id_variable]
  }
  
}






steps <- 1000
c.index <- 0

while ( steps > 0 & c.index < 0.7){
  steps <- steps - 1
  surv_obj <- Surv(tiempo, estado)
  surv_fit <- coxph(surv_obj ~ . , data = datos)
  c.index <- as.vector(surv_fit$concordance["concordance"])
}
