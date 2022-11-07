library("smotefamily")
library("readr")
library("dplyr")
library("magrittr")
library("data.table")
library("caTools")

##############
# Trae datos #
##############

data <- read.csv("support2.csv") %>% as.data.table(.)

data_1 <- data %>% filter(dzclass == "Cancer")
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
A <- data_1[,-c(2)] %>% as.data.frame()
B <- data_1[,c(2)] %>% as.data.frame()
data_1_balanced <- SMOTE(A,B)

sample <- sample.split(data_1$death, SplitRatio = 0.8)
train  <- subset(data_1, sample == TRUE)
test   <- subset(data_1, sample == FALSE)

train$death %>% table %>% prop.table
test$death  %>% table %>% prop.table


SMOTed <- SMOTE(melanoma[,-c(3)],melanoma[,c(3)], K=1)
SMOTed <- SMOTed$data
DBSMOTed <- DBSMOTE(melanoma[,-c(3)],melanoma[,c(3)])
DBSMOTed <- DBSMOTed$data

# DBSMOTed %>% select("class") %>% table(.) %>% prop.table(.) %>% print(.)# DBSMOTed %>% select("class") %>% table(.) %>% prop.table(.) %>% print(.)


