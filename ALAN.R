# ============================================================
'R code for the paper "Disrupted darkness: 
The impact of artificial light at night on melatonin secretion of 
Hermodice carunculata (Polychaeta, Annelida)" by Keklikoglou et al.

Christina Pavloudi
christina.pavloudi@embrc.eu
https://cpavloud.github.io/mysite/

	Copyright (C) 2025 Christina Pavloudi
  
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
  
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.'

# =============================================================

################################################################################
############################ LOAD LIBRARIES ####################################
################################################################################

#list of CRAN packages needed
.packages = c("tidyverse", "ggforce", "rstatix", "report", "ggpubr")

#install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

#load packages into session 
lapply(.packages, require, character.only=TRUE)

packageVersion("tidyverse"); citation("tidyverse")
packageVersion("ggforce"); citation("ggforce")
packageVersion("rstatix"); citation("rstatix")
packageVersion("report"); citation("report")
packageVersion("ggpubr"); citation("ggpubr")

R.Version()
citation()

#define a default theme for ggplot graphics
theme_set(theme_bw()) 

#set seed for reproducibility
set.seed(12346)

################################################################################
######################### LOAD AND FORMAT DATA #################################
################################################################################

#import the tables with the concentrations
#each table contains three columns: "Time", "Concentration", "Treatment"
tryptamine <- read.csv("tryptamine.csv", sep = ",", header=TRUE)
serotonin <- read.csv("serotonin.csv", sep = ",", header=TRUE)
melatonin <- read.csv("melatonin.csv", sep = ",", header=TRUE)

#reorder based on Time
tryptamine$Time <- factor(tryptamine$Time , levels=c("8:00", "11:00", "14:00", "17:00", "20:00",
                                                     "23:00", "2:00" ,"5:00", "8:00 "))
serotonin$Time <- factor(serotonin$Time , levels=c("8:00", "11:00", "14:00", "17:00", "20:00",
                                                  "23:00", "2:00" ,"5:00", "8:00 "))
melatonin$Time <- factor(melatonin$Time , levels=c("8:00", "11:00", "14:00", "17:00", "20:00",
                                                   "23:00", "2:00" ,"5:00", "8:00 "))

#create subsets for each treatment
tryptamine_LD <- subset(tryptamine, Treatment == "L:D")
tryptamine_LL <- subset(tryptamine, Treatment == "L:L")

serotonin_LD <- subset(serotonin, Treatment == "L:D")
serotonin_LL <- subset(serotonin, Treatment == "L:L")

melatonin_LD <- subset(melatonin, Treatment == "L:D")
melatonin_LL <- subset(melatonin, Treatment == "L:L")

#create subsets for each 12 hour cycle
tryptamine_12H_8_20 <- subset(tryptamine, Time == "8:00" | Time == "11:00" | Time == "14:00"
                              | Time == "17:00" | Time == "20:00")
tryptamine_12H_20_8 <- subset(tryptamine, Time == "20:00" | Time == "23:00" | Time == "2:00"
                              | Time == "5:00" | Time == "8:00 ")
serotonin_12H_8_20 <- subset(serotonin, Time == "8:00" | Time == "11:00" | Time == "14:00"
                             | Time == "17:00" | Time == "20:00")
serotonin_12H_20_8 <- subset(serotonin, Time == "20:00" | Time == "23:00" | Time == "2:00"
                             | Time == "5:00" | Time == "8:00 ")
melatonin_12H_8_20 <- subset(melatonin, Time == "8:00" | Time == "11:00" | Time == "14:00"
                             | Time == "17:00" | Time == "20:00")
melatonin_12H_20_8 <- subset(melatonin, Time == "20:00" | Time == "23:00" | Time == "2:00"
                             | Time == "5:00" | Time == "8:00 ")

#split the melatonin_12H_8_20 according to Treatment
melatonin_12H_8_20_LD <- subset(melatonin_12H_8_20, Treatment == "L:D")
melatonin_12H_8_20_LL <- subset(melatonin_12H_8_20, Treatment == "L:L")

################################################################################
################################# TWO-WAY ANOVA ################################
################################################################################

################################################################################
################################### melatonin ##################################
################################################################################

#before checking the normality assumption, we first need to compute the ANOVA 
melatonin_two_way_anova <- stats::aov(Concentration ~ Treatment * Time, data = melatonin)
#check the summary of the ANOVA results
summary(melatonin_two_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(melatonin_two_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(melatonin_two_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(melatonin_two_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(melatonin_two_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(Concentration ~ Treatment * Time, data = melatonin)
#run the Levene’s test for homogeneity of variances
car::leveneTest(Concentration ~ Treatment * Time, data = melatonin)

#check if there are outliers
melatonin %>%
  group_by(Treatment, Time) %>%
  rstatix::identify_outliers(Concentration)

#interpretations of ANOVA results
report::report(melatonin_two_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(melatonin_two_way_anova)

#log transform the melatonin concentration for melatonin_12H_8_20
melatonin_12H_8_20_log <- melatonin_12H_8_20 %>%
  mutate(log_concentration = log2(Concentration))

#run anova on log transformed values of melatonin_12H_8_20
melatonin_12H_8_20_log_two_way_anova <- stats::aov(log_concentration ~ Treatment * Time, data = melatonin_12H_8_20_log)
#check the summary of the ANOVA results
summary(melatonin_12H_8_20_log_two_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(melatonin_12H_8_20_log_two_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(melatonin_12H_8_20_log_two_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(melatonin_12H_8_20_log_two_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(melatonin_12H_8_20_log_two_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(log_concentration ~ Treatment * Time, data = melatonin_12H_8_20_log)
#run the Levene’s test for homogeneity of variances
car::leveneTest(log_concentration ~ Treatment * Time, data = melatonin_12H_8_20_log)

#check if there are outliers
melatonin_12H_8_20_log %>%
  group_by(Treatment, Time) %>%
  rstatix::identify_outliers(log_concentration)

#interpretations of ANOVA results
report::report(melatonin_12H_8_20_log_two_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(melatonin_12H_8_20_log_two_way_anova)

################################################################################
################################### serotonin ##################################
################################################################################

#before checking the normality assumption, we first need to compute the ANOVA 
serotonin_two_way_anova <- stats::aov(Concentration ~ Treatment * Time, data = serotonin)
#check the summary of the ANOVA results
summary(serotonin_two_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(serotonin_two_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(serotonin_two_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(serotonin_two_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(serotonin_two_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(Concentration ~ Treatment * Time, data = serotonin)
#run the Levene’s test for homogeneity of variances
car::leveneTest(Concentration ~ Treatment * Time, data = serotonin)

#check if there are outliers
serotonin %>%
  group_by(Treatment, Time) %>%
  rstatix::identify_outliers(Concentration)

#interpretations of ANOVA results
report::report(serotonin_two_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(serotonin_two_way_anova)


################################################################################
################################### tryptamine #################################
################################################################################

#before checking the normality assumption, we first need to compute the ANOVA 
tryptamine_two_way_anova <- stats::aov(Concentration ~ Treatment * Time, data = tryptamine)
#check the summary of the ANOVA results
summary(tryptamine_two_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(tryptamine_two_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(tryptamine_two_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(tryptamine_two_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(tryptamine_two_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(Concentration ~ Treatment * Time, data = tryptamine)
#run the Levene’s test for homogeneity of variances
car::leveneTest(Concentration ~ Treatment * Time, data = tryptamine)

#check if there are outliers
tryptamine %>%
  group_by(Treatment, Time) %>%
  rstatix::identify_outliers(Concentration)

#interpretations of ANOVA results
report::report(tryptamine_two_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(tryptamine_two_way_anova)


#log transform the melatonin concentration for tryptamine
tryptamine_log <- tryptamine %>%
  mutate(log_concentration = log2(Concentration))

#run anova on log transformed values of tryptamine
tryptamine_log_two_way_anova <- stats::aov(log_concentration ~ Treatment * Time, data = tryptamine_log)
#check the summary of the ANOVA results
summary(tryptamine_log_two_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(tryptamine_log_two_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(tryptamine_log_two_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(tryptamine_log_two_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(tryptamine_log_two_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(log_concentration ~ Treatment * Time, data = tryptamine_log)
#run the Levene’s test for homogeneity of variances
car::leveneTest(log_concentration ~ Treatment * Time, data = tryptamine_log)

#check if there are outliers
tryptamine_log %>%
  group_by(Treatment, Time) %>%
  rstatix::identify_outliers(log_concentration)

#interpretations of ANOVA results
report::report(tryptamine_log_two_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(tryptamine_log_two_way_anova)


#log transform the melatonin concentration for tryptamine_12H_8_20
tryptamine_12H_8_20_log <- tryptamine_12H_8_20 %>%
  mutate(log_concentration = log2(Concentration))

#run anova on log transformed values of tryptamine_12H_8_20
tryptamine_12H_8_20_log_two_way_anova <- stats::aov(log_concentration ~ Treatment * Time, data = tryptamine_12H_8_20_log)
#check the summary of the ANOVA results
summary(tryptamine_12H_8_20_log_two_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(tryptamine_12H_8_20_log_two_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(tryptamine_12H_8_20_log_two_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(tryptamine_12H_8_20_log_two_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(tryptamine_12H_8_20_log_two_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(log_concentration ~ Treatment * Time, data = tryptamine_12H_8_20_log)
#run the Levene’s test for homogeneity of variances
car::leveneTest(log_concentration ~ Treatment * Time, data = tryptamine_12H_8_20_log)

#check if there are outliers
tryptamine_12H_8_20_log %>%
  group_by(Treatment, Time) %>%
  rstatix::identify_outliers(log_concentration)

#interpretations of ANOVA results
report::report(tryptamine_12H_8_20_log_two_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(tryptamine_12H_8_20_log_two_way_anova)

################################################################################
################################# ONE-WAY ANOVA ################################
################################################################################

################################################################################
################################### melatonin ##################################
################################################################################

#since we have established that Treatment and Time are significant 
#factors for melatonin_12H_8_20, we will use one-way ANOVA and 
#the respective post-hoc tests, to identify which are the time points
#significantly affecting each treatment separately

#before checking the normality assumption, we first need to compute the ANOVA 
melatonin_12H_8_20_LD_one_way_anova <- stats::aov(Concentration ~ Time, data = melatonin_12H_8_20_LD)
#check the summary of the ANOVA results
summary(melatonin_12H_8_20_LD_one_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(melatonin_12H_8_20_LD_one_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(melatonin_12H_8_20_LD_one_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(melatonin_12H_8_20_LD_one_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(melatonin_12H_8_20_LD_one_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(Concentration ~ Time, data = melatonin_12H_8_20_LD)
#run the Levene’s test for homogeneity of variances
car::leveneTest(Concentration ~ Time, data = melatonin_12H_8_20_LD)

#check if there are outliers
melatonin_12H_8_20_LD %>%
  group_by(Time) %>%
  rstatix::identify_outliers(Concentration)

#interpretations of ANOVA results
report::report(melatonin_12H_8_20_LD_one_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(melatonin_12H_8_20_LD_one_way_anova)

#log transform the melatonin concentration for melatonin_12H_8_20_LD
melatonin_12H_8_20_LD_log <- melatonin_12H_8_20_LD %>%
  mutate(log_concentration = log2(Concentration))

#run anova on log transformed values of melatonin_12H_8_20_LD
melatonin_12H_8_20_LD_log_one_way_anova <- stats::aov(log_concentration ~ Time, data = melatonin_12H_8_20_LD_log)
#check the summary of the ANOVA results
summary(melatonin_12H_8_20_LD_log_one_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(melatonin_12H_8_20_LD_log_one_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(melatonin_12H_8_20_LD_log_one_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(melatonin_12H_8_20_LD_log_one_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(melatonin_12H_8_20_LD_log_one_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(log_concentration ~ Time, data = melatonin_12H_8_20_LD_log)
#run the Levene’s test for homogeneity of variances
car::leveneTest(log_concentration ~ Time, data = melatonin_12H_8_20_LD_log)

#check if there are outliers
melatonin_12H_8_20_LD_log %>%
  group_by(Time) %>%
  rstatix::identify_outliers(log_concentration)

#interpretations of ANOVA results
report::report(melatonin_12H_8_20_LD_log_one_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(melatonin_12H_8_20_LD_log_one_way_anova)

#save the results of the post hoc tests in the appropriate format
#so that they can be added in the box plots later on
Tukey_melatonin_12H_8_20_LD_log_one_way_anova <- TukeyHSD(melatonin_12H_8_20_LD_log_one_way_anova, which = 'Time')
Tukey_melatonin_12H_8_20_LD_log_one_way_anova <- as.data.frame(Tukey_melatonin_12H_8_20_LD_log_one_way_anova$Time)
#add a group1 column
group1 <- c("8:00", "8:00", "8:00", "8:00", 
            "11:00", "11:00", "11:00", 
            "14:00", "14:00", 
            "17:00")
Tukey_melatonin_12H_8_20_LD_log_one_way_anova$group1 <- group1
#add a group2 column
group2 <- c("11:00", "14:00", "17:00", "20:00", 
            "14:00", "17:00", "20:00",
            "17:00", "20:00",
            "20:00")
Tukey_melatonin_12H_8_20_LD_log_one_way_anova$group2 <- group2
#add a significance column
Tukey_melatonin_12H_8_20_LD_log_one_way_anova <- Tukey_melatonin_12H_8_20_LD_log_one_way_anova %>%
  mutate(significance = case_when(`p adj` <= 0.01 ~ "**", 
                                  `p adj` <= 0.05 ~ "*"))

#before checking the normality assumption, we first need to compute the ANOVA 
melatonin_12H_8_20_LL_one_way_anova <- stats::aov(Concentration ~ Time, data = melatonin_12H_8_20_LL)
#check the summary of the ANOVA results
summary(melatonin_12H_8_20_LL_one_way_anova)

#run the Kolmogorov-Smirnov test for the normality of residuals
DescTools::LillieTest(melatonin_12H_8_20_LL_one_way_anova$residuals) 
#run the Shapiro-Wilk test for the normality of residuals
stats::shapiro.test(melatonin_12H_8_20_LL_one_way_anova$residuals)

#check normality visually using a histogram of the residuals
hist(melatonin_12H_8_20_LL_one_way_anova$residuals)
#check normality visually using a QQ-plot of the residuals
car::qqPlot(melatonin_12H_8_20_LL_one_way_anova$residuals,
            id = FALSE) # id = FALSE to remove point identification

#check the homogeneity of variances using a boxplot
boxplot(Concentration ~ Time, data = melatonin_12H_8_20_LL)
#run the Levene’s test for homogeneity of variances
car::leveneTest(Concentration ~ Time, data = melatonin_12H_8_20_LL)

#check if there are outliers
melatonin_12H_8_20_LL %>%
  group_by(Time) %>%
  rstatix::identify_outliers(Concentration)

#interpretations of ANOVA results
report::report(melatonin_12H_8_20_LL_one_way_anova)

#run post-hoc tests (Tukey HSD test)
stats::TukeyHSD(melatonin_12H_8_20_LL_one_way_anova)


#save the results of the post hoc tests in the appropriate format
#so that they can be added in the box plots later on
Tukey_melatonin_12H_8_20_LL_one_way_anova <- TukeyHSD(melatonin_12H_8_20_LL_one_way_anova, which = 'Time')
Tukey_melatonin_12H_8_20_LL_one_way_anova <- as.data.frame(Tukey_melatonin_12H_8_20_LL_one_way_anova$Time)
#add a group1 column
group1 <- c("8:00", "8:00", "8:00", "8:00", 
            "11:00", "11:00", "11:00", 
            "14:00", "14:00", 
            "17:00")
Tukey_melatonin_12H_8_20_LL_one_way_anova$group1 <- group1
#add a group2 column
group2 <- c("11:00", "14:00", "17:00", "20:00", 
            "14:00", "17:00", "20:00",
            "17:00", "20:00",
            "20:00")
Tukey_melatonin_12H_8_20_LL_one_way_anova$group2 <- group2
#add a significance column
Tukey_melatonin_12H_8_20_LL_one_way_anova <- Tukey_melatonin_12H_8_20_LL_one_way_anova %>%
  mutate(significance = case_when(`p adj` <= 0.01 ~ "**", 
                                  `p adj` <= 0.05 ~ "*"))


################################################################################
#################################### BOXPLOTS ##################################
################################################################################

################################################################################
################################### melatonin ##################################
################################################################################

#specify the y.position for each comparison
Tukey_melatonin_12H_8_20_LD_log_one_way_anova <- Tukey_melatonin_12H_8_20_LD_log_one_way_anova %>%
  mutate(y.position = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
Tukey_melatonin_12H_8_20_LL_one_way_anova <- Tukey_melatonin_12H_8_20_LL_one_way_anova %>%
  mutate(y.position = c(3, 4, 5, 6, 7, 8, 9, 10, 3, 12))

#plot the LD plot
melatonin_LD_plot <-
  ggboxplot(melatonin_LD, x = "Time", y = "Concentration", 
            fill = "Treatment", palette =c("powderblue"), 
            xlab = "Time",
            ylab = "Melatonin Concentration (pg mg-1)", 
            legend = "right")  +
  stat_summary(
    fun = median,
    geom = 'line',linewidth = 1,
    aes(group = Treatment, colour = Treatment),
    position = position_dodge(width = 0.9) #this has to be added
  )+
  scale_color_manual(values = c('L:D'='steelblue')) +
  stat_pvalue_manual(Tukey_melatonin_12H_8_20_LD_log_one_way_anova, label = "significance") +
  grids(linetype = "dashed") + 
  coord_cartesian(ylim=c(-0.6,3.5), clip="off") +
  annotate("segment", x = 1, xend = 4, y = -0.5, yend = -0.5, colour = "orange", linewidth = 2) +
  annotate("segment", x = 5, xend = 9, y = -0.5, yend = -0.5, colour = "lightgray", linewidth = 2)  +
  annotate("text", x = 2.5, y = -0.5, label = "Light", colour = "orangered4", size = 5) +
  annotate("text", x = 7, y = -0.5, label = "Dark", colour = "dimgrey", size = 5)  
melatonin_LD_plot
ggsave("melatonin_LD_plot_publication.png", width = 10, height = 6, dpi = 600)
ggsave("melatonin_LD_plot_publication.eps", width = 10, height = 6, dpi = 600)

#plot the LL plot
melatonin_LL_plot <-
  ggboxplot(melatonin_LL, x = "Time", y = "Concentration", 
            fill = "Treatment", palette =c("lightcoral"), 
            xlab = "Time",
            ylab = "Melatonin Concentration (pg mg-1)", 
            legend = "right")  +
  font("xlab", size = 12, color = "black", family="arial")+
  font("ylab", size = 12, color = "black", family="arial")+
  stat_summary(
    fun = median,
    geom = 'line',linewidth = 1,
    aes(group = Treatment, colour = Treatment),
    position = position_dodge(width = 0.9) #this has to be added
  )+
  scale_color_manual(values = c('L:L'='indianred4')) +
  stat_pvalue_manual(Tukey_melatonin_12H_8_20_LL_one_way_anova, label = "significance") +
  grids(linetype = "dashed") + 
  coord_cartesian(ylim=c(-0.6,3.5), clip="off") +
  annotate("segment", x = 1, xend = 9, y = -0.5, yend = -0.5, colour = "orange", linewidth = 2) +
  annotate("text", x = 4.5, y = -0.5, label = "Light", colour = "orangered4", size = 5)   
melatonin_LL_plot
ggsave("melatonin_LL_plot_publication.png", width = 10, height = 6, dpi = 600)
ggsave("melatonin_LL_plot_publication.eps", width = 10, height = 6, dpi = 600)

################################################################################
################################### serotonin ##################################
################################################################################

#plot the LD plot
serotonin_LD_plot <-
  ggboxplot(serotonin_LD, x = "Time", y = "Concentration", 
            fill = "Treatment", palette =c("powderblue"), 
            xlab = "Time",
            ylab = "Serotonin Concentration (ng mg-1)", 
            legend = "right")  +
  stat_summary(
    fun = median,
    geom = 'line',linewidth = 1,
    aes(group = Treatment, colour = Treatment),
    position = position_dodge(width = 0.9) #this has to be added
  )+
  scale_color_manual(values = c('L:D'='steelblue')) +
  grids(linetype = "dashed") + 
  coord_cartesian(ylim=c(200,1200), clip="off") +
  annotate("segment", x = 1, xend = 4, y = 250, yend = 250, colour = "orange", linewidth = 2) +
  annotate("segment", x = 5, xend = 9, y = 250, yend = 250, colour = "lightgray", linewidth = 2)  +
  annotate("text", x = 2.5, y = 250, label = "Light", colour = "orangered4", size = 5) +
  annotate("text", x = 7, y = 250, label = "Dark", colour = "dimgrey", size = 5)  
serotonin_LD_plot
ggsave("serotonin_LD_plot_publication.png", width = 10, height = 6, dpi = 600)
ggsave("serotonin_LD_plot_publication.eps", width = 10, height = 6, dpi = 600)

#plot the LL plot
serotonin_LL_plot <-
  ggboxplot(serotonin_LL, x = "Time", y = "Concentration", 
            fill = "Treatment", palette =c("lightcoral"), 
            xlab = "Time",
            ylab = "Serotonin Concentration (ng mg-1)", 
            legend = "right")  +
  stat_summary(
    fun = median,
    geom = 'line',linewidth = 1,
    aes(group = Treatment, colour = Treatment),
    position = position_dodge(width = 0.9) #this has to be added
  )+
  scale_color_manual(values = c('L:L'='indianred4')) +
  grids(linetype = "dashed") + 
  coord_cartesian(ylim=c(200,1200), clip="off") +
  annotate("segment", x = 1, xend = 9, y = 250, yend = 250, colour = "orange", linewidth = 2) +
  annotate("text", x = 4.5, y = 250, label = "Light", colour = "orangered4", size = 5)  
serotonin_LL_plot
ggsave("serotonin_LL_plot_publication.png", width = 10, height = 6, dpi = 600)
ggsave("serotonin_LL_plot_publication.eps", width = 10, height = 6, dpi = 600)


################################################################################
################################### tryptamine #################################
################################################################################

#plot the LD plot
tryptamine_LD_plot <-
  ggboxplot(tryptamine_LD, x = "Time", y = "Concentration", 
            fill = "Treatment", palette =c("powderblue"), 
            xlab = "Time",
            ylab = "Tryptamine Concentration (ng mg-1)", 
            legend = "right")  +
  stat_summary(
    fun = median,
    geom = 'line',linewidth = 1,
    aes(group = Treatment, colour = Treatment),
    position = position_dodge(width = 0.9) #this has to be added
  )+
  scale_color_manual(values = c('L:D'='steelblue')) +
  grids(linetype = "dashed") + 
  coord_cartesian(ylim=c(-4,180), clip="off") +
  annotate("segment", x = 1, xend = 4, y = -1, yend = -1, colour = "orange", linewidth = 2) +
  annotate("segment", x = 5, xend = 9, y = -1, yend = -1, colour = "lightgray", linewidth = 2)  +
  annotate("text", x = 2.5, y = -2, label = "Light", colour = "orangered4", size = 5) +
  annotate("text", x = 7, y = -2, label = "Dark", colour = "dimgrey", size = 5)  
tryptamine_LD_plot
ggsave("tryptamine_LD_plot_publication.png", width = 10, height = 6, dpi = 600)
ggsave("tryptamine_LD_plot_publication.eps", width = 10, height = 6, dpi = 600)

#plot the LL plot
tryptamine_LL_plot <-
  ggboxplot(tryptamine_LL, x = "Time", y = "Concentration", 
            fill = "Treatment", palette =c("lightcoral"), 
            xlab = "Time",
            ylab = "Tryptamine Concentration (ng mg-1)", 
            legend = "right")  +
  stat_summary(
    fun = median,
    geom = 'line',linewidth = 1,
    aes(group = Treatment, colour = Treatment),
    position = position_dodge(width = 0.9) #this has to be added
  )+
  scale_color_manual(values = c('L:L'='indianred4')) +
  grids(linetype = "dashed") + 
  coord_cartesian(ylim=c(-4,180), clip="off") +
  annotate("segment", x = 1, xend = 9, y = -1, yend = -1, colour = "orange", linewidth = 2) +
  annotate("text", x = 4.5, y = -2, label = "Light", colour = "orangered4", size = 5) 
tryptamine_LL_plot
ggsave("tryptamine_LL_plot_publication.png", width = 10, height = 6, dpi = 600)
ggsave("tryptamine_LL_plot_publication.eps", width = 10, height = 6, dpi = 600)

################################################################################
################################################################################
################################################################################

#save workspace
save.image("ALAN.RData")

################################################################################
################################################################################
################################################################################
