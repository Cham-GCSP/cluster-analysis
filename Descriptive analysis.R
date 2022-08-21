
setwd("D:/Graduation")

library(pacman)

pacman::p_load(rlang, dplyr,psych,plyr,sf,rgdal,sp,stringr,rgdal,ggplot2,tidyr,scales,RColorBrewer,
               gridExtra,readr,rgeos,mtcars,ggsn,AICmodavg,readxl,ggpubr)

#_______________________Description analysis__________________________

# Load shapefile

shap2010 <- st_read('Data/Shapefile_2000_2021/Shapefile_2000_2021/Sgg_2010.shp')

# Sigungu code change file

Sgg_change <- read.csv('Data/Sigungu_list_1990_2021/Sgg_change.csv')

Sgg_change <- subset(Sgg_change, select = -c(5,6))

code_change1 <- Sgg_change %>% 
  filter(from_year > 1998 & from_year < 2014) %>% 
  dplyr::select(fromcode, tocode)

#____________________________________________________________________
#Population

population <- read.csv('Data/Data/Pop_1999_2013.csv')


pop9913 <- population %>%   # calculate population by period
  mutate(Sgg_cd = plyr::mapvalues(Sigungu_cd, code_change1$fromcode, code_change1$tocode),
        year_agg = cut(Year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE)) %>%
  select(.,-c(Year, Sigungu_cd)) %>% 
  group_by(Sgg_cd, Gender, year_agg) %>% 
  summarise_all(sum,na.rm=T) %>%
  ungroup

pop9913$Sgg_cd <- as.character(pop9913$Sgg_cd)

#__________________________________________________________________
# Mortality
Mor <- read.csv('Data/Data/Mor_1999_2013.csv')

#create city name
Mor[Mor =='Malignant neoplasm of stomach (C16)'] <- 'stomach'
Mor[Mor =='Malignant neoplasms of colon, rectum and anus (C18-C21)'] <- 'colon'
Mor[Mor =='Malignant neoplasms of liver and intrahepatic bile ducts (C22)'] <- 'liver'
Mor[Mor =='Malignant neoplasms of trachea, bronchus and lung (C33-C34)'] <- 'lung'
Mor[Mor =='Malignant neoplasm of breast (C50)'] <- 'breast'
Mor[Mor =='Malignant neoplasm of uterus (C53-C55)'] <- 'uterus'
Mor[Mor =='Malignant neoplasm of prostate (C61)'] <- 'prostate'
Mor[Mor=='Deaths[Person]'] <- 'D'
Mor[Mor=='Age-standardized death rate[per 100 thousand person]'] <- 'ASMR'

Mor$sgg_cd <- as.character(Mor$sgg_cd)

Mor <- Mor %>% 
  mutate(Sgg_cd = plyr::mapvalues(sgg_cd, code_change1$fromcode, code_change1$tocode)) %>%
  pivot_longer(cols = Y1999:Y2013) %>% 
  mutate(Year = as.integer(str_sub(name,2,5))) %>% 
  mutate(value1 = replace(value,is.na(value),0),
        year_agg = cut(Year, breaks = c(1999,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE)) %>%
  pivot_wider(id_cols = c(Sgg_cd, cause,gender, year_agg), names_from = item, values_from =  value1, values_fn = sum) %>% 
  group_by(Sgg_cd,cause,gender, year_agg) %>% 
  summarise_all(sum) %>%
  right_join(.,pop9913, by = c('Sgg_cd'='Sgg_cd', 'gender' = 'Gender','year_agg' = 'year_agg')) 

Mor <- na.omit(Mor)
# calculate expected mortality adjusted for age

# Merge observed mortality data with expected mortality
mor_standard <- read_excel('D:/Graduation/Data/Data/Data_standard.xlsx', sheet = 'mor_standard')

pop_agespecific_mor <- read_excel('D:/Graduation/Data/Data/Data_standard.xlsx', sheet = 'pop_agespecific_mor')

#_________________ Calculate expected mortality______________
#1. merge standard mortality with study population

ex_mor <- merge(pop_agespecific_mor,mor_standard, by = c('gender', 'Age'))

#2. calculate expected mortality for each age-specific & Aggregate data into very 5 years

ex_mor <- pivot_longer(data = ex_mor, cols = 4:18, names_to = "year", values_to = "study_pop")

ex_mor1 <- ex_mor %>% 
  mutate(ex_mor = `mor_rate (per hundred thousands)`*study_pop/100000) %>% 
  mutate(Year = as.integer(str_sub(year,2,5)), 
         year_agg = cut(Year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE)) %>% 
  select(., -c(2,5,6,7,8,9,11)) %>% 
  group_by(Sgg_cd,gender, cause, year_agg) %>% 
  summarise_all(sum)


# Merge
Mor_data <- merge(ex_mor1, Mor, by = c('Sgg_cd','gender','cause','year_agg'))

write.csv(Mor_data,'Mor_data.csv')


#________________________________________________________________________

# Descriptive mortality

library(dplyr)
detach(package:plyr)

Mor$cause <- as.factor(Mor$cause)
Mor$gender <- as.factor(Mor$gender)
Mor$year_agg <- as.factor(Mor$year_agg)

Mor_summary <- Mor %>% 
  group_by(gender,cause, year_agg) %>% 
  summarize(D_sum  = round(sum(D,na.rm = T), digits = 2),
            D_mean = round(mean(D, na.rm = T), digits = 2),
            D_min  = round(min(D, na.rm = T), digits = 2),
            D_max  = round(max(D, na.rm = T), digits = 2),
            D_sd   = round(sd(D, na.rm = T), digits = 2))
Mor_summary

write.csv(Mor_summary, 'Mor_summary.csv')

# ANOVA test 
anova_func = function(c,g){
  data_anova <- Mor %>% 
    filter(cause == c & gender == g)
  onewaymor <- aov(D ~ year_agg, data = data_anova)

  summary(onewaymor)
  
}

anova_func("colon", "Female")
anova_func("colon", "Male")

anova_func("liver", "Female")
anova_func("liver", "Male")

anova_func("stomach", "Female")
anova_func("stomach", "Male")

anova_func("lung", "Female")
anova_func("lung", "Male")

anova_func("breast", "Female")
anova_func("uterus", "Female")
anova_func("prostate", "Male")


# Spatial distribution of cancer mortality over 3 periods, using shapefile 2010

# For 2 sexes

Mor_data$Sgg_cd <- as.character(Mor_data$Sgg_cd)

map_func1 <- function(c,y,g,t){
  qbrks = seq(0, 1, length.out = 5)
  
  dat1 <- shap2010 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,Mor_data, by = 'Sgg_cd') %>% 
    filter(cause == c & gender!= 'Total' & year_agg == y& gender == g) %>% 
    mutate(asmr_cut = cut(ASMR, breaks = quantile(ASMR, breaks = qbrks, na.rm = T),
                                     include.lowest = T ))
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = asmr_cut)) +
    scale_fill_brewer(name="ASMR", palette = "YlOrRd",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    facet_wrap(~ gender, ncol = 2) +
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t) +
    north(shap2010, symbol = 3)+
    theme(legend.position = c(.9, .3),
          legend.text=element_text(size=40),
          legend.title=element_text(size=40))
  
}


Lung_mor_male1 <- map_func1('lung','1999-2003','Male','Lung cancer-period 1')

Lung_mor_female1 <- map_func1('lung','1999-2003','Female','Lung cancer-period 1')

Liver_mor_male1 <- map_func1('liver','1999-2003','Male','Liver cancer-period 1')

Liver_mor_female1 <- map_func1('liver','1999-2003','Female','Liver cancer-period 1')

Colon_mor_male1 <- map_func1('colon','1999-2003','Male','Colorectal cancer-period 1')

Colon_mor_female1 <- map_func1('colon','1999-2003','Female','Colorectal cancer-period 1')

Stomach_mor_male1 <- map_func1('stomach','1999-2003','Male','Stomach cancer-period 1')

Stomach_mor_female1 <- map_func1('stomach','1999-2003','Female','Stomach cancer-period 1')

Lung_mor_male2 <- map_func1('lung','2004-2008','Male','Lung cancer-period 2')

Lung_mor_female2 <- map_func1('lung','2004-2008','Female','Lung cancer-period 2')

Liver_mor_male2 <- map_func1('liver','2004-2008','Male','Liver cancer-period 2')

Liver_mor_female2 <- map_func1('liver','2004-2008','Female','Liver cancer-period 2')

Colon_mor_male2 <- map_func1('colon','2004-2008','Male','Colorectal cancer-period 2')

Colon_mor_female2 <- map_func1('colon','2004-2008','Female','Colorectal cancer-period 2')

Stomach_mor_male2 <- map_func1('stomach','2004-2008','Male','Stomach cancer-period 2')

Stomach_mor_female2 <- map_func1('stomach','2004-2008','Female','Stomach cancer-period 2')

Lung_mor_male3 <- map_func1('lung','2009-2013','Male','Lung cancer-period 3')

Lung_mor_female3 <- map_func1('lung','2009-2013','Female','Lung cancer-period 3')

Liver_mor_male3 <- map_func1('liver','2009-2013','Male','Liver cancer-period 3')

Liver_mor_female3 <- map_func1('liver','2009-2013','Female','Liver cancer-period 3')

Colon_mor_male3 <- map_func1('colon','2009-2013','Male','Colorectal cancer-period 3')

Colon_mor_female3 <- map_func1('colon','2009-2013','Female','Colorectal cancer-period 3')

Stomach_mor_male3 <- map_func1('stomach','2009-2013','Male','Stomach cancer-period 3')

Stomach_mor_female3 <- map_func1('stomach','2009-2013','Female','Stomach cancer-period 3')

# mapping asmr for 1 sex cancer

map_func2 <- function(c,g,y,t){
 
  #quantile breaks
  qbrks = seq(0, 1, length.out = 5)
  
  dat1 <- shap2010 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,Mor_data, by = 'Sgg_cd') %>% 
    filter(cause == c & gender == g & year_agg == y) %>% 
    mutate(asmr_cut = cut(ASMR, breaks = quantile(ASMR, breaks = qbrks, na.rm = T),
                          include.lowest = T ))
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = asmr_cut)) +
    scale_fill_brewer(name="ASMR", palette = "YlOrRd",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    facet_wrap(~ year_agg, ncol = 3)+
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          legend.position = c(0.9, 0.2)) +
    labs( title = t)+
    north(shap2010, symbol = 3)
  
    

}

Breast_map1 <- map_func2('breast','Female', '1999-2003','Female breast cancer-Period 1')

Uterus_map1 <- map_func2('uterus','Female','1999-2003', 'Uterus cancer-period 1')

Prostate_map1 <- map_func2('prostate','Male', '1999-2003','Prostate cancer-period 1')

Breast_map2 <- map_func2('breast','Female', '2004-2008','Female breast cancer-Period 2')

Uterus_map2 <- map_func2('uterus','Female','2004-2008', 'Uterus cancer-Period 2')

Prostate_map2 <- map_func2('prostate','Male', '2004-2008','Prostate cancer-Period 2')

Breast_map3 <- map_func2('breast','Female', '2009-2013','Female breast cancer-Period 3')

Uterus_map3 <- map_func2('uterus','Female','2009-2013', 'Uterus cancer-Period 3')

Prostate_map3 <- map_func2('prostate','Male', '2009-2013','Prostate cancer-Period 3')

# plot map of spatial mortality
# male
png(filename = "ASMR_male.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)
ggarrange(Lung_mor_male1,
          Lung_mor_male2,
          Lung_mor_male3,
          Stomach_mor_male1,
          Stomach_mor_male2,
          Stomach_mor_male3,
          Colon_mor_male1,
          Colon_mor_male2,
          Colon_mor_male3,
          Liver_mor_male1,
          Liver_mor_male2,
          Liver_mor_male3,
          Prostate_map1,
          Prostate_map2,
          Prostate_map3,
          ncol = 3, nrow = 5)
dev.off()

png(filename = "ASMR_female.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)
ggarrange(Lung_mor_female1,
          Lung_mor_female2,
          Lung_mor_female3,
          Stomach_mor_female1,
          Stomach_mor_female2,
          Stomach_mor_female3,
          Colon_mor_female1,
          Colon_mor_female2,
          Colon_mor_female3,
          Liver_mor_female1,
          Liver_mor_female2,
          Liver_mor_female3,
          Breast_map1,
          Breast_map2,
          Breast_map3,
          Uterus_map1,
          Uterus_map2,
          Uterus_map3,
          ncol = 3, nrow = 6)
dev.off()


#____________Spatial distribution of population across 3 periods____

  qbrks = seq(0, 1, length.out = 5)
  
  pop_func = function(g,y,t){
    
    dat1 <- shap2010 %>% 
      mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
      left_join(.,Mor, by = 'Sgg_cd') %>% 
      filter(cause == 'lung' & gender == g & year_agg == y) %>% 
      mutate(pop1  = Population/100000,
             pop_cut = cut(pop1, breaks = quantile(pop1, breaks = qbrks, na.rm = T),
                           include.lowest = T ))

    ggplot(data = dat1) +
      geom_sf(data = dat1, aes(fill = pop_cut)) +
      scale_fill_brewer(name="Population", palette = "Blues",
                        guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
      theme(text = element_text(size = 30))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),
            legend.position = c(0.9, 0.2)) +
      labs( title = t,
            caption = 'Data from KOSIS website\n unit: hundred thousands') +
      north(shap2010, symbol = 3)
    
  }
 
  
png(filename = 'Male_pop1', 
      units = "in", 
      width = 10, 
      height = 10, 
      res=300)

pop_male <- pop_func('Male','1999-2003','Male - Period 1')

pop_male

dev.off()

png(filename = 'Male_pop2', 
    units = "in", 
    width = 10, 
    height = 10, 
    res=300)
pop_male2 <- pop_func('Male','2004-2008','Male- Period 2')

pop_male2

dev.off()
  
png(filename = 'Male_pop3', 
    units = "in", 
    width = 10, 
    height = 10, 
    res=300)
pop_male3 <- pop_func('Male','2009-2013','Male- Period 3')

pop_male3

dev.off()

png(filename = 'Female_pop1', 
    units = "in", 
    width = 10, 
    height = 10, 
    res=300)
pop_Female <- pop_func('Female','1999-2003','Female- Period 1')

pop_Female

dev.off()

png(filename = 'Female_pop2', 
    units = "in", 
    width = 10, 
    height = 10, 
    res=300)
pop_Female2 <- pop_func('Female','2004-2008','Female- Period 2')

pop_Female2

dev.off()

png(filename = 'Female_pop3', 
    units = "in", 
    width = 10, 
    height = 10, 
    res=300)
pop_Female3 <- pop_func('Female','2009-2013','Female- Period 3')

pop_Female3

dev.off()



# Boxplot
  
boxfunc1 = function(c,l,u,t){
  
  boxdata <- Mor_data %>% 
    filter(cause == c & gender != 'Total') 
  
  ggplot(data = boxdata,
          aes(x = year_agg, y = D, fill = gender))+
    geom_boxplot(outlier.shape = NA)+
    coord_cartesian(ylim =  c(l,u))+
    theme_bw()+ 
    theme(plot.title = element_text(size = 35,hjust = 0.5, face = "bold"),
          axis.text.x = element_text(hjust = 1, size = 30),
          axis.title.x = element_text(size = 30),
          axis.text.y  = element_text( size = 30),
          axis.title.y = element_text(size=30) )+
    theme(legend.position = c(.90, .90),
          legend.text=element_text(size=30),
          legend.title=element_text(size=30))+
    labs(title = t,
         x = '',
         y = 'Death (persons)') +
    theme( axis.ticks.length.y = unit(.30, "cm"),
           axis.ticks.length.x = unit(-.30, "cm"),
           axis.text.x = element_text(margin = margin(t = .3, unit = "cm") ) )+
    guides(fill= guide_legend("Gender")) 
   

}

# run function

Lung_count <- boxfunc1('lung',0,400,'Lung cancer death counts across three periods')

Stomach_count <- boxfunc1('stomach',0,300, 'Stomach cancer death counts across three periods')

Colon_count <- boxfunc1('colon',0,250,'Colorectal cancer death counts across three periods')

Liver_count<- boxfunc1('liver', 0,300,'Liver cancer death counts across three periods')


# 1 genders cancer
boxfunc2 <- function(c,g,l,u,t){
  boxdata <- Mor_data %>% 
    filter(cause == c & gender == g)  
  

  ggplot(data = boxdata,
        aes(x = year_agg, y = D, fill = gender))+
    geom_boxplot(outlier.shape = NA)+
    coord_cartesian(ylim =  c(l,u))+
    theme_bw()+
    theme(plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(hjust = 1, size = 30),
          axis.title.x=element_text(size = 30),
          axis.text.y = element_text( size = 30),
          axis.title.y=element_text(size = 30))+
    labs(title = t,
         x = '',
         y = 'Death (persons)') +
    theme( axis.ticks.length.y = unit(.40, "cm"),
           axis.ticks.length.x = unit(-.40, "cm"),
           axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
    theme(legend.position = "none")
}

Breast_count <- boxfunc2('breast', 'Female',0,75,'Female breast cancer death counts across three periods')

Uterus_count <- boxfunc2('uterus', 'Female',0,60,'Female uterus cancer death counts across three periods')

Prostate_count<- boxfunc2('prostate','Male', 0,60,'Male prostate cancer death counts across three periods')

# plot boxplots
png(filename = 'Box_mor', 
    units = "in", 
    width = 40, 
    height = 30, 
    res=300)
ggarrange(Lung_count,
          Stomach_count,
          Colon_count,
          Liver_count,
          Breast_count,
          Uterus_count,
          Prostate_count,
          ncol = 3,nrow = 3)
dev.off()


# _______________Cancer incidence descriptive___________________________

Inc <- read.csv('Data/Data/cancerInc_sgg.csv',encoding = 'UTF-8' )

colnames(Inc) <- c('gender','cause','sgg_cd','year_agg','item','value')

Inc$gender <- plyr::mapvalues(Inc$gender, from = c('???','??????','??????'), to = c('Total','Male','Female'))

Inc$cause <- plyr::mapvalues(Inc$cause, from = c('???(C16)','??????(C18-C20)','???(C22)','???(C33-C34)','??????(C50)','????????????(C53)','?????????(C61)',"?????????(C73)"),
                             to = c('stomach', 'colon', 'liver', 'lung', 'breast', 'uterus', 'prostate', 'thyroid'))

Inc$item <- plyr::mapvalues(Inc$item, from = c("????????????","????????????","????????????????????????"),
                            to = c('counts','rate','ASIR'))

Inc$sgg_cd <- as.character(Inc$sgg_cd)

Inc1 <- Inc %>% 
  mutate(Sgg_cd = plyr::mapvalues(sgg_cd, code_change1$fromcode, code_change1$tocode)) %>% 
  group_by(Sgg_cd, cause, gender, item, year_agg) %>% 
  select(.,-3) %>% 
  summarise_all(sum,na.rm=T) %>%
  data.frame() %>% 
  pivot_wider(names_from = 'item',
              values_from = 'value')

# compute expected incidence adjust for age

inc_standard <- read_excel('D:/Graduation/Data/Data/Data_standard.xlsx', sheet = 'inc_standard')

pop_agespecific_inc <- read_excel('D:/Graduation/Data/Data/Data_standard.xlsx', sheet = 'pop_agespecific_inc')

ex_inc <- merge(pop_agespecific_inc,inc_standard, by = c('gender', 'Age'))

#2. calculate expected mortality for each age-specific

ex_inc <- pivot_longer(data = ex_inc, cols = 4:18, names_to = "year", values_to = "study_pop")

ex_inc1 <- ex_inc %>% 
  mutate(ex_inc = `inc_rate (per hundred thousands)`*study_pop/100000) %>% 
  mutate(Year = as.integer(str_sub(year,2,5)), 
         year_agg = cut(Year, breaks = c(1998,2003,2008,2013), labels = c('1999-2003', '2004-2008', '2009-2013'), right = TRUE)) %>% 
  select(., -c(2,5,6,7,8,9,11)) %>% 
  group_by(Sgg_cd,gender,cause, year_agg) %>% 
  summarise_all(sum)

# merge with observe data

ex_inc1$Sgg_cd <- as.character(ex_inc1$Sgg_cd)

inc_data <- merge(Inc1, ex_inc1, by = c('Sgg_cd','cause','gender','year_agg'))

#merge with population

inc_data <- right_join(inc_data, pop9913, by = c('Sgg_cd'='Sgg_cd', 'gender' = 'Gender','year_agg' = 'year_agg'))


# Descriptive incidence

Inc_sumary <- inc_data %>% 
  group_by(gender,cause, year_agg) %>%
  summarise(I_sum = round(sum(counts, na.rm = T), digits = 2),
            I_mean = round(mean(counts, na.rm = T), digits = 2),
            I_min = round(min(counts, na.rm = T), digits = 2),
            I_max = round(max(counts, na.rm = T), digits = 2),
            I_sd  = round(sd(counts, na.rm = T), digits = 2))
           

write.csv(Inc_sumary, 'Incidence_summary.csv')

# anova 
anova_func_inc = function(c,g){
  data_anova <- inc_data %>% 
    filter(cause == c & gender == g)
  onewaymor <- aov(counts ~ year_agg, data = data_anova)
  
  summary(onewaymor)
  
}

anova_func_inc("colon", "Female")
anova_func_inc("colon", "Male")

anova_func_inc("liver", "Female")
anova_func_inc("liver", "Male")

anova_func_inc("stomach", "Female")
anova_func_inc("stomach", "Male")

anova_func_inc("lung", "Female")
anova_func_inc("lung", "Male")

anova_func_inc("breast", "Female")
anova_func_inc("uterus", "Female")
anova_func_inc("prostate", "Male")

anova_func_inc("thyroid", "Female")
anova_func_inc("thyroid", "Male")


# Mapping ASIR across 3 periods

map_func3 <- function(c,y,g,t){
  
  qbrks = seq(0, 1, length.out = 5)
  
  dat1 <- shap2010 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,inc_data, by = 'Sgg_cd') %>% 
    filter(cause == c & gender == g & year_agg == y) %>% 
    mutate(asir_cut = cut(ASIR, breaks = quantile(ASIR, breaks = qbrks, na.rm = T),
                          include.lowest = T ))
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = asir_cut)) +
    scale_fill_brewer(name="ASIR", palette = "Reds",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t) +
    north(shap2010, symbol = 3)+
    theme(legend.position = c(.9, .3),
          legend.text=element_text(size=40),
          legend.title=element_text(size=40))
  
}

Lung_inc_male1 <- map_func3('lung','1999-2003','Male','Lung cancer-period 1')

Lung_inc_female1 <- map_func3('lung','1999-2003','Female','Lung cancer-period 1')

Liver_inc_male1 <- map_func3('liver','1999-2003','Male','Liver cancer-period 1')

Liver_inc_female1 <- map_func3('liver','1999-2003','Female','Liver cancer-period 1')

Colon_inc_male1 <- map_func3('colon','1999-2003','Male','Colorectal cancer-period 1')

Colon_inc_female1 <- map_func3('colon','1999-2003','Female','Colorectal cancer-period 1')

Stomach_inc_male1 <- map_func3('stomach','1999-2003','Male','Stomach cancer-period 1')

Stomach_inc_female1 <- map_func3('stomach','1999-2003','Female','Stomach cancer-period 1')

Thyroid_inc_male1 <- map_func3('thyroid','1999-2003','Male','Thyroid cancer-period 1')

Thyroid_inc_female1 <- map_func3('thyroid','1999-2003','Female','Thyroid cancer-period 1')

Lung_inc_male2 <- map_func3('lung','2004-2008','Male','Lung cancer-period 2')

Lung_inc_female2 <- map_func3('lung','2004-2008','Female','Lung cancer-period 2')

Liver_inc_male2 <- map_func3('liver','2004-2008','Male','Liver cancer-period 2')

Liver_inc_female2 <- map_func3('liver','2004-2008','Female','Liver cancer-period 2')

Colon_inc_male2 <- map_func3('colon','2004-2008','Male','Colorectal cancer-period 2')

Colon_inc_female2 <- map_func3('colon','2004-2008','Female','Colorectal cancer-period 2')

Stomach_inc_male2 <- map_func3('stomach','2004-2008','Male','Stomach cancer-period 2')

Stomach_inc_female2 <- map_func3('stomach','2004-2008','Female','Stomach cancer-period 2')

Thyroid_inc_male2 <- map_func3('thyroid','2004-2008','Male','Thyroid cancer-period 2')

Thyroid_inc_female2 <- map_func3('thyroid','2004-2008','Female','Thyroid cancer-period 2')

Lung_inc_male3 <- map_func3('lung','2009-2013','Male','Lung cancer-period 3')

Lung_inc_female3 <- map_func3('lung','2009-2013','Female','Lung cancer-period 3')

Liver_inc_male3 <- map_func3('liver','2009-2013','Male','Liver cancer-period 3')

Liver_inc_female3 <- map_func3('liver','2009-2013','Female','Liver cancer-period 3')

Colon_inc_male3 <- map_func3("colon",'2009-2013','Male','Colorectal cancer-period 3')

Colon_inc_female3 <- map_func3('colon','2009-2013','Female','Colorectal cancer-period 3')

Stomach_inc_male3 <- map_func3('stomach','2009-2013','Male','Stomach cancer-period 3')

Stomach_inc_female3 <- map_func3('stomach','2009-2013','Female','Stomach cancer-period 3')

Thyroid_inc_male3 <- map_func3('thyroid','2009-2013','Male','Thyroid cancer-period 3')

Thyroid_inc_female3 <- map_func3('thyroid','2009-2013','Female','Thyroid cancer-period 3')


# mapping asir for 1 sex cancer

map_func4 <- function(c,g,y,t){
  
  #quantile breaks
  qbrks = seq(0, 1, length.out = 5)
  
  dat1 <- shap2010 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,inc_data, by = 'Sgg_cd') %>% 
    filter(cause == c, gender == g, year_agg == y) %>% 
    mutate(asir_cut = cut(ASIR, breaks = quantile(ASIR, breaks = qbrks, na.rm = T),
                          include.lowest = T ))
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = asir_cut)) +
    scale_fill_brewer(name="ASIR", palette = "Reds",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t)+
    north(shap2010, symbol = 3)+
    theme(legend.position = c(.9, .3),
          legend.text=element_text(size=40),
          legend.title=element_text(size=40))
  
}

Breast_inc_map1 <- map_func4('breast','Female', '1999-2003','Female breast cancer-Period 1')

Uterus_inc_map1 <- map_func4('uterus','Female','1999-2003', 'Uterus cancer-period 1')

Prostate_inc_map1 <- map_func4('prostate','Male', '1999-2003','Prostate cancer-period 1')

Breast_inc_map2 <- map_func4('breast','Female', '2004-2008','Female breast cancer-Period 2')

Uterus_inc_map2 <- map_func4('uterus','Female','2004-2008', 'Uterus cancer-Period 2')

Prostate_inc_map2 <- map_func4('prostate','Male', '2004-2008','Prostate cancer-Period 2')

Breast_inc_map3 <- map_func4('breast','Female', '2009-2013','Female breast cancer-Period 3')

Uterus_inc_map3 <- map_func4('uterus','Female','2009-2013', 'Uterus cancer-Period 3')

Prostate_inc_map3 <- map_func4('prostate','Male', '2009-2013','Prostate cancer-Period 3')
#plot
# male
png(filename = "ASIR_male.png", 
    units = "in", 
    width = 30, 
    height = 60, 
    res = 300)
ggarrange(Lung_inc_male1,
          Lung_inc_male2,
          Lung_inc_male3,
          Stomach_inc_male1,
          Stomach_inc_male2,
          Stomach_inc_male3,
          Colon_inc_male1,
          Colon_inc_male2,
          Colon_inc_male3,
          Liver_inc_male1,
          Liver_inc_male2,
          Liver_inc_male3,
          Prostate_inc_map1,
          Prostate_inc_map2,
          Prostate_inc_map3,
          Thyroid_inc_male1,
          Thyroid_inc_male2,
          Thyroid_inc_male3,
          ncol = 3, nrow = 6)
dev.off()

png(filename = "ASIR_female.png", 
    units = "in", 
    width = 30, 
    height = 70, 
    res = 300)
ggarrange(Lung_inc_female1,
          Lung_inc_female2,
          Lung_inc_female3,
          Stomach_inc_female1,
          Stomach_inc_female2,
          Stomach_inc_female3,
          Colon_inc_female1,
          Colon_inc_female2,
          Colon_inc_female3,
          Liver_inc_female1,
          Liver_inc_female2,
          Liver_inc_female3,
          Breast_inc_map1,
          Breast_inc_map2,
          Breast_inc_map3,
          Uterus_inc_map1,
          Uterus_inc_map2,
          Uterus_inc_map3,
          Thyroid_inc_female1,
          Thyroid_inc_female2,
          Thyroid_inc_female3,
          ncol = 3, nrow = 7)
dev.off()


#_____________Boxplot for cancer incidence__________________

#Incidence count

boxfunc_inc1 <- function(c,l,u,t){
  select <- inc_data %>% 
    filter(cause == c) %>% 
    ggplot(data = .,
         mapping = aes(x = year_agg, y = counts, fill = gender))+
    geom_boxplot(outlier.shape = NA)+
    coord_cartesian(ylim =  c(l,u))+
    theme_bw()+
    theme(plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
          axis.text.x = element_text( hjust = 1, size = 30),
          axis.title.x=element_text(size=30),
          axis.text.y = element_text( size = 30),
          axis.title.y=element_text(size=30))+
    labs(title = t,
         x = '',
         y = 'Incidence (cases)')+
    theme( axis.ticks.length.y = unit(.40, "cm"),
           axis.ticks.length.x = unit(-.40, "cm"),
           axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
    theme(legend.position = c(0.1,0.9),
          legend.text=element_text(size=30),
          legend.title=element_text(size=30))
}

Lung_inc <- boxfunc_inc1('lung', 0, 500, 'Lung cancer incidence across three periods')

Liver_inc <- boxfunc_inc1('liver', 0, 500, 'Liver cancer incidence across three periods')

Stomach_inc <- boxfunc_inc1('stomach', 0, 700,'Stomach cancer incidence across three periods')

Colon_inc <- boxfunc_inc1('colon', 0, 600, 'Colorectal cancer incidence across three periods')

Thyroid_inc <- boxfunc_inc1('thyroid', 0,1200, 'Thyroid cancer incidence across three periods')

#  for 1 sex cancer

boxfunc_inc2 <- function(c,g,l,u,t){
  select <- inc_data %>% 
    filter(cause == c & gender == g) %>% 
    ggplot(data =.,
         mapping = aes(x = year_agg, y = counts, fill = gender))+
    geom_boxplot(outlier.shape = NA, notch = T, alpha=0.8)+
    coord_cartesian(ylim =  c(l,u))+
    theme_bw()+
    theme(plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 30, hjust = 1),
          axis.text.y = element_text(size = 30),
          axis.title.x=element_text(size=30),
          axis.title.y=element_text(size=30))+
    labs(title = t,
         x = '',
         y = 'Incidence (cases)')+
    theme( axis.ticks.length.y = unit(.40, "cm"),
           axis.ticks.length.x = unit(-.40, "cm"),
           axis.text.x = element_text(margin = margin(t = .5, unit = "cm")))+
    theme(legend.position= 'none')
  
  
}

Breast_inc <- boxfunc_inc2('breast','Female', 0,700,'Female breast cancer incidence across three periods')

Uterus_inc <- boxfunc_inc2('uterus', 'Female', 0,200,'Female uterus cancer incidence across three periods')

Prostate_inc <- boxfunc_inc2('prostate', 'Male', 0, 300,'Male prostate cancer incidence across three periods')

png(filename = 'Box_inc', 
    units = "in", 
    width = 40, 
    height = 30, 
    res=300)
ggarrange(Lung_inc,
          Stomach_inc,
          Colon_inc,
          Liver_inc,
          Breast_inc,
          Uterus_inc,
          Prostate_inc,
          Thyroid_inc,
          ncol = 3,nrow = 3)
dev.off()


write.csv(inc_data,'inc_data.csv')


#________Spatial distribution of health risk factors___________________

riskfactor <- read.csv('Data/Data/riskfactor.csv')

screen <- read.csv('Data/Data/screening.csv')

ses <- read.csv('Data/Data/ses.csv', encoding = 'UTF-8')

sgg_cd <- read.csv('Data/Data/Sigungu_code.csv')

sgg_cd10kr <- read.csv('Data/Data/sgg10.csv', encoding = 'UTF-8')

riskfactor_cd <- merge(riskfactor, sgg_cd, by = 'sgg_nm', all.x = T)

ses <- merge(ses, sgg_cd10kr, by = 'sgg_nm')

ses <- subset(ses, select = -c(1,2,11,13,14))

screen <- merge(screen, sgg_cd, by = 'sgg_nm')

factors <- merge(screen, ses, by = 'sgg_cd', all = T)

factors <- merge(factors, riskfactor_cd, by = 'sgg_cd', all = T)

factors$sgg_cd <- as.character(factors$sgg_cd)

shap2010 <- st_read('Data/Data/Sgg_2010.shp')

shap2010_sub <- shap2010[-c(29,59,201,202,227,233,244),]

data_factor = shap2010_sub %>% 
  left_join(factors, by = c('SIGUNGU_CD' = 'sgg_cd'))


# mapping distribution of health risk factors
factors <- factors %>% 
  select(-19) %>% 
  pivot_longer(cols = 3:25, names_to = 'factor' ) 


map_risk <- function(c,n){
  
  qbrks = seq(0, 1, length.out = 5)
  
  dat <- shap2010 %>% 
    mutate(sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,factors, by = 'sgg_cd') %>%
    filter(factor == c) %>% 
    mutate(cut = cut(value, breaks = quantile(value, breaks = qbrks, na.rm = T),
                          include.lowest = T ))
  
  ggplot(data = dat) +
    geom_sf(data = dat, aes(fill = cut)) +
    scale_fill_brewer(name = n, palette = "YlOrRd",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    north(shap2010, symbol = 3)+
    theme(legend.position = c(.9, .3),
          legend.text=element_text(size=40),
          legend.title=element_text(size=40))
  
}

sto_scr <- map_risk('P_stomach_m', "Screening rate of stomach cancer")
colon_scr <- map_risk('P_colon_m', "Screening rate of colorectal cancer")
liver_scr <- map_risk('P_liver_m', "Screening rate of liver cancer")
bre_scr <- map_risk('P_breast', "Screening rate of breast cancer")
cer_scr <- map_risk('P_cervical', "Screening rate of cervical cancer")

png(filename = 'screening_map', 
    units = "in", 
    width = 50, 
    height = 10, 
    res=300)

ggarrange(sto_scr, colon_scr, liver_scr, bre_scr, cer_scr,
          ncol = 5,nrow = 1)
dev.off()

grdp_map <- map_risk('GRDP_2005', 'GRDP')
edu_map <- map_risk('edu_2005', "Higher education level")
elder_map <- map_risk('elderly_2005', "Elderly rate")

png(filename = 'ses_map', 
    units = "in", 
    width = 30, 
    height = 10, 
    res=300)

ggarrange(grdp_map, edu_map, elder_map,
          ncol = 3,nrow = 1)
dev.off()

drink_map <- map_risk('dinking', 'Drinking rate')
walking_map <- map_risk('walking', "Walking rate")
stress_map <- map_risk('stress', "Stress rate")
obesity_map <- map_risk('obesity', "Obesity rate")
salt_map <- map_risk('lowsaltuse', "Low salt use rate")
smoke_map <- map_risk('smoking', "Smoking rate")

png(filename = 'healthrisk_map', 
    units = "in", 
    width = 60, 
    height = 10, 
    res=300)

ggarrange(drink_map, walking_map, stress_map, obesity_map, salt_map, smoke_map,
          ncol = 6,nrow = 1)
dev.off()


# mapping crude incidence rate

map_func5 <- function(c,y,g,t){
  
  qbrks = seq(0, 1, length.out = 5)
  
  dat1 <- shap2010 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,inc_data, by = 'Sgg_cd') %>% 
    filter(cause == c & gender == g & year_agg == y) %>% 
    mutate(rate_crude = counts *100000/Population)%>%
    mutate(rate_crude_cut = cut(rate_crude, breaks = quantile(rate_crude, breaks = qbrks, na.rm = T),
                                include.lowest = T ))
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = rate_crude_cut)) +
    scale_fill_brewer(name="Incidence rate", palette = "Reds",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t) +
    north(shap2010, symbol = 3)+
    theme(legend.position = c(.9, .3),
          legend.text=element_text(size=40),
          legend.title=element_text(size=40))
  
}

Lung_inc_male1_r <- map_func5('lung','1999-2003','Male','Lung cancer-period 1')

Lung_inc_female1_r <- map_func5('lung','1999-2003','Female','Lung cancer-period 1')

Liver_inc_male1_r <- map_func5('liver','1999-2003','Male','Liver cancer-period 1')

Liver_inc_female1_r <- map_func5('liver','1999-2003','Female','Liver cancer-period 1')

Colon_inc_male1_r <- map_func5('colon','1999-2003','Male','Colorectal cancer-period 1')

Colon_inc_female1_r <- map_func5('colon','1999-2003','Female','Colorectal cancer-period 1')

Stomach_inc_male1_r <- map_func5('stomach','1999-2003','Male','Stomach cancer-period 1')

Stomach_inc_female1_r <- map_func5('stomach','1999-2003','Female','Stomach cancer-period 1')

Thyroid_inc_male1_r <- map_func5('thyroid','1999-2003','Male','Thyroid cancer-period 1')

Thyroid_inc_female1_r <- map_func5('thyroid','1999-2003','Female','Thyroid cancer-period 1')

Lung_inc_male2_r <- map_func5('lung','2004-2008','Male','Lung cancer-period 2')

Lung_inc_female2_r <- map_func5('lung','2004-2008','Female','Lung cancer-period 2')

Liver_inc_male2_r <- map_func5('liver','2004-2008','Male','Liver cancer-period 2')

Liver_inc_female2_r <- map_func5('liver','2004-2008','Female','Liver cancer-period 2')

Colon_inc_male2_r <- map_func5('colon','2004-2008','Male','Colorectal cancer-period 2')

Colon_inc_female2_r <- map_func5('colon','2004-2008','Female','Colorectal cancer-period 2')

Stomach_inc_male2_r <- map_func5('stomach','2004-2008','Male','Stomach cancer-period 2')

Stomach_inc_female2_r <- map_func5('stomach','2004-2008','Female','Stomach cancer-period 2')

Thyroid_inc_male2_r <- map_func5('thyroid','2004-2008','Male','Thyroid cancer-period 2')

Thyroid_inc_female2_r <- map_func5('thyroid','2004-2008','Female','Thyroid cancer-period 2')

Lung_inc_male3_r <- map_func5('lung','2009-2013','Male','Lung cancer-period 3')

Lung_inc_female3_r <- map_func5('lung','2009-2013','Female','Lung cancer-period 3')

Liver_inc_male3_r <- map_func5('liver','2009-2013','Male','Liver cancer-period 3')

Liver_inc_female3_r <- map_func5('liver','2009-2013','Female','Liver cancer-period 3')

Colon_inc_male3_r <- map_func5("colon",'2009-2013','Male','Colorectal cancer-period 3')

Colon_inc_female3_r <- map_func5('colon','2009-2013','Female','Colorectal cancer-period 3')

Stomach_inc_male3_r <- map_func5('stomach','2009-2013','Male','Stomach cancer-period 3')

Stomach_inc_female3_r <- map_func5('stomach','2009-2013','Female','Stomach cancer-period 3')

Thyroid_inc_male3_r <- map_func5('thyroid','2009-2013','Male','Thyroid cancer-period 3')

Thyroid_inc_female3_r <- map_func5('thyroid','2009-2013','Female','Thyroid cancer-period 3')


# incidence for 1 sex cancer

map_func6 <- function(c,g,y,t){
  
  #quantile breaks
  qbrks = seq(0, 1, length.out = 5)
  
  dat1 <- shap2010 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,inc_data, by = 'Sgg_cd') %>% 
    filter(cause == c, gender == g, year_agg == y) %>% 
    mutate(rate_crude = counts *100000/Population)%>%
    mutate(rate_crude_cut = cut(rate_crude, breaks = quantile(rate_crude, breaks = qbrks, na.rm = T),
                                include.lowest = T ))
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = rate_crude_cut)) +
    scale_fill_brewer(name="Incidence rate", palette = "Reds",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t)+
    north(shap2010, symbol = 3)+
    theme(legend.position = c(.9, .3),
          legend.text=element_text(size=40),
          legend.title=element_text(size=40))
  
}

Breast_inc_map1_r <- map_func6('breast','Female', '1999-2003','Female breast cancer-Period 1')

Uterus_inc_map1_r <- map_func6('uterus','Female','1999-2003', 'Uterus cancer-period 1')

Prostate_inc_map1_r <- map_func6('prostate','Male', '1999-2003','Prostate cancer-period 1')

Breast_inc_map2_r <- map_func6('breast','Female', '2004-2008','Female breast cancer-Period 2')

Uterus_inc_map2_r <- map_func6('uterus','Female','2004-2008', 'Uterus cancer-Period 2')

Prostate_inc_map2_r <- map_func6('prostate','Male', '2004-2008','Prostate cancer-Period 2')

Breast_inc_map3_r <- map_func6('breast','Female', '2009-2013','Female breast cancer-Period 3')

Uterus_inc_map3_r <- map_func6('uterus','Female','2009-2013', 'Uterus cancer-Period 3')

Prostate_inc_map3_r <- map_func6('prostate','Male', '2009-2013','Prostate cancer-Period 3')

# male

png(filename = "incidence_crude_rate_male.png", 
    units = "in", 
    width = 40, 
    height = 60, 
    res = 300)
ggarrange(Lung_inc_male1_r,
          Lung_inc_male2_r,
          Lung_inc_male3_r,
          Stomach_inc_male1_r,
          Stomach_inc_male2_r,
          Stomach_inc_male3_r,
          Colon_inc_male1_r,
          Colon_inc_male2_r,
          Colon_inc_male3_r,
          Liver_inc_male1_r,
          Liver_inc_male2_r,
          Liver_inc_male3_r,
          Prostate_inc_map1_r,
          Prostate_inc_map2_r,
          Prostate_inc_map3_r,
          Thyroid_inc_male1_r,
          Thyroid_inc_male2_r,
          Thyroid_inc_male3_r,
          ncol = 3, nrow = 6)
dev.off()



# incidence rate adjust for age

map_func7 <- function(c,y,g,t){
  
  qbrks = seq(0, 1, length.out = 5)
  
  dat1 <- shap2010 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,inc_data, by = 'Sgg_cd') %>% 
    filter(cause == c & gender == g & year_agg == y) %>% 
    mutate(rate_crude = ex_inc *100000/Population)%>%
    mutate(rate_crude_cut = cut(rate_crude, breaks = quantile(rate_crude, breaks = qbrks, na.rm = T),
                                include.lowest = T ))
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = rate_crude_cut)) +
    scale_fill_brewer(name="Incidence rate", palette = "Reds",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t) +
    north(shap2010, symbol = 3)+
    theme(legend.position = c(.9, .3),
          legend.text=element_text(size=40),
          legend.title=element_text(size=40))
  
}

Lung_inc_male1_e <- map_func7('lung','1999-2003','Male','Lung cancer-period 1')

Lung_inc_female1_e <- map_func7('lung','1999-2003','Female','Lung cancer-period 1')

Liver_inc_male1_e <- map_func7('liver','1999-2003','Male','Liver cancer-period 1')

Liver_inc_female1_e <- map_func7('liver','1999-2003','Female','Liver cancer-period 1')

Colon_inc_male1_e <- map_func7('colon','1999-2003','Male','Colorectal cancer-period 1')

Colon_inc_female1_e <- map_func7('colon','1999-2003','Female','Colorectal cancer-period 1')

Stomach_inc_male1_e <- map_func7('stomach','1999-2003','Male','Stomach cancer-period 1')

Stomach_inc_female1_e <- map_func7('stomach','1999-2003','Female','Stomach cancer-period 1')

Thyroid_inc_male1_e <- map_func7('thyroid','1999-2003','Male','Thyroid cancer-period 1')

Thyroid_inc_female1_e <- map_func7('thyroid','1999-2003','Female','Thyroid cancer-period 1')

Lung_inc_male2_e <- map_func7('lung','2004-2008','Male','Lung cancer-period 2')

Lung_inc_female2_e <- map_func7('lung','2004-2008','Female','Lung cancer-period 2')

Liver_inc_male2_e <- map_func7('liver','2004-2008','Male','Liver cancer-period 2')

Liver_inc_female2_e <- map_func7('liver','2004-2008','Female','Liver cancer-period 2')

Colon_inc_male2_e <- map_func7('colon','2004-2008','Male','Colorectal cancer-period 2')

Colon_inc_female2_e <- map_func7('colon','2004-2008','Female','Colorectal cancer-period 2')

Stomach_inc_male2_e <- map_func7('stomach','2004-2008','Male','Stomach cancer-period 2')

Stomach_inc_female2_e <- map_func7('stomach','2004-2008','Female','Stomach cancer-period 2')

Thyroid_inc_male2_e <- map_func7('thyroid','2004-2008','Male','Thyroid cancer-period 2')

Thyroid_inc_female2_e <- map_func7('thyroid','2004-2008','Female','Thyroid cancer-period 2')

Lung_inc_male3_e <- map_func7('lung','2009-2013','Male','Lung cancer-period 3')

Lung_inc_female3_e <- map_func7('lung','2009-2013','Female','Lung cancer-period 3')

Liver_inc_male3_e <- map_func7('liver','2009-2013','Male','Liver cancer-period 3')

Liver_inc_female3_e <- map_func7('liver','2009-2013','Female','Liver cancer-period 3')

Colon_inc_male3_e <- map_func7("colon",'2009-2013','Male','Colorectal cancer-period 3')

Colon_inc_female3_e <- map_func7('colon','2009-2013','Female','Colorectal cancer-period 3')

Stomach_inc_male3_e <- map_func7('stomach','2009-2013','Male','Stomach cancer-period 3')

Stomach_inc_female3_e <- map_func7('stomach','2009-2013','Female','Stomach cancer-period 3')

Thyroid_inc_male3_e <- map_func7('thyroid','2009-2013','Male','Thyroid cancer-period 3')

Thyroid_inc_female3_e <- map_func7('thyroid','2009-2013','Female','Thyroid cancer-period 3')


# ________for 1 sex_________________

map_func8 <- function(c,g,y,t){
  
  #quantile breaks
  qbrks = seq(0, 1, length.out = 5)
  
  dat1 <- shap2010 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,inc_data, by = 'Sgg_cd') %>% 
    filter(cause == c, gender == g, year_agg == y) %>% 
    mutate(rate_crude = ex_inc *100000/Population)%>%
    mutate(rate_crude_cut = cut(rate_crude, breaks = quantile(rate_crude, breaks = qbrks, na.rm = T),
                                include.lowest = T ))
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = rate_crude_cut)) +
    scale_fill_brewer(name="Incidence rate", palette = "Reds",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    theme(text = element_text(size = 40))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t)+
    north(shap2010, symbol = 3)+
    theme(legend.position = c(.9, .3),
          legend.text=element_text(size=40),
          legend.title=element_text(size=40))
  
}

Breast_inc_map1_e <- map_func8('breast','Female', '1999-2003','Female breast cancer-Period 1')

Uterus_inc_map1_e <- map_func8('uterus','Female','1999-2003', 'Uterus cancer-period 1')

Prostate_inc_map1_e <- map_func8('prostate','Male', '1999-2003','Prostate cancer-period 1')

Breast_inc_map2_e <- map_func8('breast','Female', '2004-2008','Female breast cancer-Period 2')

Uterus_inc_map2_e <- map_func8('uterus','Female','2004-2008', 'Uterus cancer-Period 2')

Prostate_inc_map2_e <- map_func8('prostate','Male', '2004-2008','Prostate cancer-Period 2')

Breast_inc_map3_e <- map_func8('breast','Female', '2009-2013','Female breast cancer-Period 3')

Uterus_inc_map3_e <- map_func8('uterus','Female','2009-2013', 'Uterus cancer-Period 3')

Prostate_inc_map3_e <- map_func8('prostate','Male', '2009-2013','Prostate cancer-Period 3')



png(filename = "incidence_ex_rate_male.png", 
    units = "in", 
    width = 40, 
    height = 60, 
    res = 300)
ggarrange(Lung_inc_male1_e,
          Lung_inc_male2_e,
          Lung_inc_male3_e,
          Stomach_inc_male1_e,
          Stomach_inc_male2_e,
          Stomach_inc_male3_e,
          Colon_inc_male1_e,
          Colon_inc_male2_e,
          Colon_inc_male3_e,
          Liver_inc_male1_e,
          Liver_inc_male2_e,
          Liver_inc_male3_e,
          Prostate_inc_map1_e,
          Prostate_inc_map2_e,
          Prostate_inc_map3_e,
          Thyroid_inc_male1_e,
          Thyroid_inc_male2_e,
          Thyroid_inc_male3_e,
          ncol = 3, nrow = 6)
dev.off()

