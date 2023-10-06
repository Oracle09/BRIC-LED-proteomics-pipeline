# Import raw data into R
library(readxl)
mem_sp_hw = read_excel("C:\\Users\\romip\\Dropbox\\SharedWithAl\\ControlChecklists\\Mem_Hardware_ANOVA.xlsx", sheet = "HW")

#apply package tidyr to arrange data from rows into columns
library(tidyr)
mem_sp_hw1 = tidyr::gather(mem_sp_hw, key = Treatment, value = Measurement, 3:11)

mem_sp_hw1$Measurement = as.numeric(mem_sp_hw1$Measurement)
mem_sp_hw1$Treatment = as.factor(mem_sp_hw1$Treatment)

#Order according to sequence ID
mem_sp_hw2 = mem_sp_hw1[order(mem_sp_hw1$`Sequence Id`), ]

#Grouping the replicates according to the treatment done
library(stringr)
mem_sp_hw2$Treatment = str_trim(mem_sp_hw2$Treatment)

mem_sp_hw2$Group[mem_sp_hw2$Treatment == "OH1"] = "HW"
mem_sp_hw2$Group[mem_sp_hw2$Treatment == "OH2"] = "HW"
mem_sp_hw2$Group[mem_sp_hw2$Treatment == "OH3"] = "HW"

mem_sp_hw2$Group[mem_sp_hw2$Treatment == "S1"] = "Space"
mem_sp_hw2$Group[mem_sp_hw2$Treatment == "S2"] = "Space"
mem_sp_hw2$Group[mem_sp_hw2$Treatment == "S3"] = "Space"

mem_sp_hw2$Group[mem_sp_hw2$Treatment == "GC1"] = "Ground"
mem_sp_hw2$Group[mem_sp_hw2$Treatment == "GC2"] = "Ground"
mem_sp_hw2$Group[mem_sp_hw2$Treatment == "GC3"] = "Ground"

table(mem_sp_hw2$Group)

mem_sp_hw2$Measurement <- as.numeric(mem_sp_hw2$Measurement)
mem_sp_hw2$`Sequence Name` = as.character(mem_sp_hw2$`Sequence Name`)

#using package dplyr to run ANOVA on each protein and gene one by one
library(dplyr)
library(broom)
mem_sp_hw3 = data.frame(Id = mem_sp_hw2$`Sequence Id`, Measurement = mem_sp_hw2$Measurement, Group = factor(mem_sp_hw2$Group))
mem_sp_hw3$Id = as.character(mem_sp_hw3$Id)

#the result from the ANOVA is stored in obj1
obj1 <- mem_sp_hw3 %>% group_by(Id) %>% do(broom::tidy(model = anova(lm(Measurement ~ Group, data = mem_sp_hw3))))
obj1$model
write.csv(mem_sp_hw3, file = "C:\\Users\\romip\\Dropbox\\SharedWithAl\\ControlChecklists\\data_check_Mem_3.csv", row.names = FALSE)

#Running a post-hoc test
library(multcomp)
summary(glht(mod1, linfct=mcp(Treatment="Dunnett")))