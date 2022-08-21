

memory.limit(size=10000000)
# Summary statistic and mapping change of screening attendance rate across year of South Korea 2010-2020

setwd(getwd())

# library
library(pacman)
library(rmarkdown)
rmarkdown::render
pacman::p_load(plyr,dplyr, tidyr, sf, data.table, psych, ggplot2, ggsn, tmap,ggpubr, tidyverse, purrr)

# Import dataset 

scr <- read.csv("D:/Graduation/Data/Data/screening_all.csv")

code <- read.csv("D:/Graduation/Data/Data/Sigungu_code.csv")

sgg_change <- read.csv('D:/Graduation/Data/Sigungu_list_1990_2021/Sgg_change.csv')

shap2016 <- st_read('D:/Graduation/Data/Shapefile_2000_2021/Shapefile_2000_2021/Sgg_2016.shp')


sgg_change <- subset(sgg_change, select = -c(5,6))

code_change1 <- sgg_change %>% 
  filter(from_year >= 2010) %>% 
  dplyr::select(fromcode, tocode)


colnames(scr) <- c("sigungu", "sgg_nm", "cancer", "item", "unit", "2010", '2011', "2012" , "2013", "2014", "2015", "2016","2017", "2018", "2019", "2020")

scr <- subset(scr, select = -c(sigungu, unit))

# changing value of columns 

scr <- scr %>% 
     mutate(sgg_nm = plyr::mapvalues (sgg_nm, from = c('Gyeonggi IisanDong-gu',"Gyeonggi lisanSeo-gu"), to =c('Gyeonggi Ilsandong-gu',"Gyeonggi Ilsanseo-gu")),
            cancer = plyr::mapvalues(cancer, from = c("Stomach Cancer", "Liver Cencer", "Colorectal Cancer", "Breast Cancer", "Cervical Cancer", "Lung Cancer"), to = c('stomach', 'liver', 'colon', 'breast', 'cervical','lung')),
            item = plyr::mapvalues(item, from = c("Number of Eligible Individuals[In person]", "Number of Actual Examinees[In person]"), to = c("eli", "act")))

scr$sgg_nm <- gsub("Chungbuk", "Chungcheongbuk", scr$sgg_nm)
scr$sgg_nm <- gsub("Chungnam", "Chungcheongnam", scr$sgg_nm)
scr$sgg_nm <- gsub("Jeonbuk", "Jeollabuk", scr$sgg_nm)
scr$sgg_nm <- gsub("Jeonnam", "Jeollanam", scr$sgg_nm)
scr$sgg_nm <- gsub("Gyeongbuk", "Gyeongsangbuk", scr$sgg_nm)
scr$sgg_nm <- gsub("Gyeongnam", "Gyeongsangnam", scr$sgg_nm)

# merge data with sgg code

scr1 <- merge (scr, code, by = "sgg_nm")

scr2 <- scr1 %>% 
  mutate(Sgg_cd = plyr::mapvalues(sgg_cd, code_change1$fromcode, code_change1$tocode)) %>% 
  select(.,-c(1,15)) %>% 
  group_by(cancer, item, Sgg_cd) %>% 
  summarise_all(sum,na.rm=T) %>%
  ungroup
  
scr3 <- scr2 %>% 
  pivot_longer(., cols = 4:14) %>% 
  pivot_wider( id_cols = c(Sgg_cd,cancer, name),
               names_from = item,
               values_from = value)

colnames(scr3)[which(names(scr3) == "name")] <- "year"

scr3 <- scr3 %>% 
  mutate(`Cancer type` = plyr::mapvalues(cancer, from = c("stomach", "colon", "liver", "breast", "cervical", "lung"),
                                         to = c("Stomach cancer", "Colorectal cancer", "Liver cancer", "Breast cancer", "Cervical cancer", "Lung cancer") ) )

# calculate rate of screening attendance

scr3$atd <- round(scr3$act*100/scr3$eli, 1)  #(number of actual screening * 100/number of eligible)

dat_sum <- scr3 %>%  # summarise annually
  select(.,-c(4,5)) %>% 
  group_by(cancer, year) %>% 
  summarise(across(where(is.numeric),
                   list(mean = ~ round(mean(.x, na.rm = TRUE), 1), 
                        count = ~ round(sum(is.na(.x)), 1),
                        sd = ~ round(sd(.x, na.rm = T),1),
                        min = ~ round(min(.x, na.rm = T),1),
                        max = ~ round(max(.x, na.rm = T),1) ))) 
write.csv(dat_sum, "dat_sum.csv")

# calculate sum of screening rate of all cancer types for each year (2010-2014-2018)

scr3_1 <- scr3 %>% 
  select(.,-c(2,7)) %>% 
  filter(year == '2010' | year == "2014" | year == "2018") %>% 
  group_by(Sgg_cd,year) %>% 
  summarise(across(where(is.numeric),
                   list(sum = ~sum(.x, na.rm = T) ) ) ) %>% 
  mutate(r_sum = round(act_sum*100/eli_sum))

# summarise rate for total cancer for each year

scr3_1_sum <- scr3_1 %>% 
  group_by(year) %>% 
  summarise(mean_sum = mean(r_sum, na.rm = T),
         median_sum = median(r_sum, na.rm = T),
         sd = sd(r_sum, na.rm = T),
         min = min(r_sum, na.rm = T),
         max = max(r_sum, na.rm = T) )
  

# How attendance rate change from 2010-2018

scr4 <- scr3 %>% 
  select(1,3,6,7) %>% 
  filter(year == "2010" | year =="2014" | year == "2018") 

# distribution of cancer screening attendance group

scr4 <- within(scr4, {
  d <- NA
  d[atd < 15] <- "0-15"
  d[atd < 30 & atd >= 15] <- "1. 15-30"
  d[atd >= 30 & atd < 45] <- "2. 30-45"
  d[atd >= 45 & atd < 60] <- "3. 45-60"
  d[atd >= 60] <- "4. >=60"
  
  d10 <- NA
  d10[atd < 52.4 & `Cancer type` == "Breast cancer"] <- "1. < 52.4"
  d10[atd >= 52.4 & `Cancer type` == "Breast cancer"] <- "2. >= 52.4"
  d10[atd < 38.75 & `Cancer type` == "Cervical cancer"] <- "1. < 38.75"
  d10[atd >= 38.75 & `Cancer type` == "Cervical cancer"] <- "2. >= 38.75"
  d10[atd < 32.9 & `Cancer type` == "Colorectal cancer"] <- "1. < 32.9"
  d10[atd >= 32.9 & `Cancer type` == "Colorectal cancer"] <- "2. >= 32.9.4"
  d10[atd < 45.65 & `Cancer type` == "Liver cancer"] <- "1. < 45.65"
  d10[atd >= 45.65 & `Cancer type` == "Liver cancer"] <- "2. >= 45.65"
  d10[atd < 46 & `Cancer type` == "Stomach cancer"] <- "1. < 46"
  d10[atd >= 46 & `Cancer type` == "Stomach cancer"] <- "2. >=46"
  
}) 

scr4$Sgg_cd <- as.character(scr4$Sgg_cd)

dat_sum1 <- scr4 %>% # summarize statistic by group of cancer type, year (Table 1)
  group_by(`Cancer type`, year) %>% 
  summarise(across(where(is.numeric),
                   list(mean = ~ round(mean(.x, na.rm = TRUE), 1), 
                        median = ~ median(.x, na.rm = T),
                        count = ~ round(sum(is.na(.x)), 1),
                        sd = ~ round(sd(.x, na.rm = T),1),
                        min = ~ round(min(.x, na.rm = T),1),
                        max = ~ round(max(.x, na.rm = T),1) )))

write.csv(dat_sum1, "dat_sum1.csv")

dat_sum2 <- scr4 %>%  # count number of sigungu by group of cancer type, year, and screening rate group (Table2)
  group_by(`Cancer type`, year, d) %>% 
  select(-3) %>% 
  summarise(count = n() )

write.csv(dat_sum2, "dat_sum2.csv")

dat_sum3 <- scr4 %>% # count number of sgg by group of median 2010
  group_by(`Cancer type`, year, d10) %>%
  select(.,-c(1,4,6)) %>% 
  summarise(count = n())

write.csv(dat_sum3, "dat_sum3.csv")
  

# Create group of screening rate for total cancer for each year

scr3_2 <- within(scr3_1, {
  ds <- NA
  ds[r_sum < 15] <- "0-15"
  ds[r_sum < 30 & r_sum >= 15] <- "1. 15-30"
  ds[r_sum >= 30 & r_sum < 45] <- "2. 30-45"
  ds[r_sum >= 45 & r_sum < 60] <- "3. 45-60"
  ds[r_sum >= 60] <- "4. >=60"
  
  ds10 <- NA
  ds10[r_sum < 42] <- "1. < 42"
  ds10[r_sum >= 42 ] <- "2. >= 42"
  
}) 


scr3_2_1 <- scr3_2 %>% # count number of sgg above and lower than total rate in 2010 for total cancer
  group_by(year, ds10) %>% 
  select(2,6) %>% 
  summarise(count = n())

# Illustrate the temporal trend of attendance rate 2010,2014,2018 (Figure 1.)

ggplot( 
  data = scr4,
  aes(x = year, y = atd,  fill = `Cancer type`))+
  geom_boxplot()+
  labs(title = "Fig 1. Screening attendance rate across years and cancer types",
       x = 'Year',
       y = 'Percentage (%)')+
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(hjust = 1, size = 15, angle = 45),
        axis.title.x=element_text(size = 15),
        axis.text.y = element_text( size = 15),
        axis.title.y=element_text(size = 15))


#________________________________________

# Create map using different scale for each cancer types for 2010,2014,2018 using different scale by cancer

# create function


map_func1 <- function(c,t) {
  
  dat1 <- shap2016 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,scr4, by = 'Sgg_cd') %>% 
    na.omit(Sgg_cd) %>% 
    filter(`Cancer type` == c & `Cancer type` != "Lung cancer")
  
  ggplot(data = dat1) +
    geom_sf(data = dat1, aes(fill = d)) +
    scale_fill_brewer(name ="Attendance rate (%)", palette = "YlOrRd",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    facet_wrap(.~ year) +
    theme(text = element_text(size = 60))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 60),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t) +
    theme(legend.position = "bottom",
          legend.text = element_text(size=50),
          legend.title = element_text(size=50))+
    north(shap2016, symbol = 3)
  
}

map_breast <- map_func1("Breast cancer", "Breast cancer")
map_cer <- map_func1("Cervical cancer", "Cervical cancer")
map_colon <- map_func1("Colorectal cancer", "Colorectal cancer")
map_liver <- map_func1("Liver cancer", "Liver cancer")
map_sto <- map_func1("Stomach cancer", "Stomach cancer")

# map for screening rate of total cancer types for each year

scr3_2$Sgg_cd <- as.character(scr3_2$Sgg_cd)

dat2 <- shap2016 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,scr3_2, by = 'Sgg_cd') %>% 
    na.omit(Sgg_cd) 

  s1 <- ggplot(data = dat2) +
    geom_sf(data = dat2, aes(fill = ds)) +
    scale_fill_brewer(name ="Attendance rate (%)", palette = "YlOrRd",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    facet_wrap(.~ year) +
    theme(text = element_text(size = 60))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 60),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = 'Total rate of all cancer types') +
    theme(legend.position = "bottom",
          legend.text = element_text(size=50),
          legend.title = element_text(size=50))+
    north(shap2016, symbol = 3)
  
# combine

png(filename = "map_diff_scale.png", 
    units = "in", 
    width = 30, 
    height = 60, 
    res = 300)

ggarrange(map_breast, map_cer, map_colon, map_liver,map_sto, s1, nrow = 6, ncol = 1)

dev.off()


# Illustrate district with lower and higher than median of rate in 2010 for each cancer types

map_func2 <- function(c, t){
  dat3 <- shap2016 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,scr4, by = 'Sgg_cd') %>% 
    na.omit(Sgg_cd) %>% 
    filter(`Cancer type` == c) 
  
  ggplot(data = dat3) +
    geom_sf(data = dat3, aes(fill = d10)) +
    scale_fill_brewer(name ="Attendance rate (%)", palette = "YlOrRd",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    facet_wrap(.~ year) +
    theme(text = element_text(size = 60))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 50),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t) +
    theme(legend.position = "bottom",
          legend.text = element_text(size=50),
          legend.title = element_text(size=50))+
    north(shap2016, symbol = 3)
  
}

map_breast_median <- map_func2("Breast cancer", "Breast cancer")
map_cer_median <- map_func2("Cervical cancer", "Cervical cancer")
map_colon_median <- map_func2("Colorectal cancer", "Colorectal cancer")
map_liver_median <- map_func2("Liver cancer", "Liver cancer")
map_sto_median <- map_func2("Stomach cancer", "Stomach cancer")
  

s2 <- ggplot(data = dat2) +
  geom_sf(data = dat2, aes(fill = ds10)) +
  scale_fill_brewer(name ="Attendance rate (%)", palette = "YlOrRd",
                    guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
  facet_wrap(.~ year) +
  theme(text = element_text(size = 60))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 60),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", color = NA)) +
  labs( title = 'Total rate of all cancer types') +
  theme(legend.position = "bottom",
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))+
  north(shap2016, symbol = 3)

  
png(filename = "map_diff_scale10.png", 
    units = "in", 
    width = 30, 
    height = 60, 
    res = 300)

ggarrange(map_breast_median, map_cer_median, map_colon_median, map_liver_median,map_sto_median, s2, nrow = 6, ncol = 1)

dev.off()

#____________________________________________________
# How screening rate change over 2 period 2010-2014 and 2010-218

scr4_change <- scr4 %>% 
  select(1,2,3,4) %>% 
  pivot_wider(id_cols = c(Sgg_cd,`Cancer type`), values_from = atd, names_from = year) %>% 
  mutate(d10_14 = `2014` - `2010`, # calculate change 
         d10_18 = `2018` - `2010`) 


# Create group of % change

scr4_change <- within(scr4_change, {
  `2010-2014` <- NA
  `2010-2014`[d10_14 < 0 ] <- "1.Decrease"
  `2010-2014`[d10_14 >= 0 & d10_14 < 5 ] <- "2.0-5"
  `2010-2014`[d10_14 >= 5 & d10_14 < 10 ] <- "3.5-10"
  `2010-2014`[d10_14 >= 10 & d10_14 < 15 ] <- "4.10-15"
  `2010-2014`[d10_14 >= 15 & d10_14 < 20 ] <- "5.15-20"
  `2010-2014`[d10_14 >= 20 & d10_14 <25 ] <- "6.20-25"
  `2010-2014`[d10_14 >= 25 & d10_14 < 30 ] <- "7.25-30"
  `2010-2014`[d10_14 >= 30 ] <- "8.>=30"
  
  `2010-2018` <- NA
  `2010-2018`[d10_18 < 0 ] <- "1.Decrease"
  `2010-2018`[d10_18 >= 0 & d10_18 < 5 ] <- "2.0-5"
  `2010-2018`[d10_18 >= 5 & d10_18 < 10 ] <- "3.5-10"
  `2010-2018`[d10_18 >= 10 & d10_18 < 15 ] <- "4.10-15"
  `2010-2018`[d10_18 >= 15 & d10_18 < 20 ] <- "5.15-20"
  `2010-2018`[d10_18 >= 20 & d10_18 <25 ] <- "6.20-25"
  `2010-2018`[d10_18 >= 25 & d10_18 < 30 ] <- "7.25-30"
  `2010-2018`[d10_18 >= 30 ] <- "8.>=30"
  
}) 

# count value of d1, d2

scr4_change1 <- scr4_change %>% 
  select(1,2,9) %>% 
  group_by(`Cancer type`) %>% 
  count(`2010-2014`) 

scr4_change2 <- scr4_change %>% 
  select(1,2,8) %>% 
  group_by(`Cancer type`) %>% 
  count(`2010-2018`)

write.csv(scr4_change1, "scr4_change1.csv")
write.csv(scr4_change2, "scr4_change2.csv")

# Illustrate the distribution of rate change across year # Fig 5

map_func3 <- function(c, t){
  dat3 <- shap2016 %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(.,scr4_change, by = 'Sgg_cd') %>% 
    na.omit(Sgg_cd) %>% 
    select(1,2,4,5,11,12,13) %>% 
    pivot_longer(cols = c(`2010-2014`, `2010-2018`), names_to = 'group') %>% 
    filter(`Cancer type` == c) 
  
  ggplot(data = dat3) +
    geom_sf(data = dat3, aes(fill = value)) +
    scale_fill_brewer(name ="Attendance rate change (%)", palette = "YlGnBu",
                      guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
    facet_wrap(.~ group) +
    theme(text = element_text(size = 60))+
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 50),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA)) +
    labs( title = t) +
    theme(legend.position = "bottom",
          legend.text = element_text(size=50),
          legend.title = element_text(size=50))+
    north(shap2016, symbol = 3)
  
}

map_breast_change <- map_func3("Breast cancer", "Breast cancer")
map_cer_change <- map_func3("Cervical cancer", "Cervical cancer")
map_colon_change <- map_func3("Colorectal cancer", "Colorectal cancer")
map_liver_change <- map_func3("Liver cancer", "Liver cancer")
map_sto_change <- map_func3("Stomach cancer", "Stomach cancer")

# overall rate change

scr3_1_change <- scr3_1 %>% 
  pivot_wider(id_cols = Sgg_cd, names_from = year, values_from = r_sum) %>% 
  mutate("2010-2014" = `2014`- `2010`,
         "2010-2018" = `2018` - `2010`) %>% 
  pivot_longer(cols = c(`2010-2014`, `2010-2018`), names_to = "period")

scr3_1_change <- within(scr3_1_change, {
  r_sum_change <- NA
  r_sum_change [value < 0 ] <- "1.Decrease"
  r_sum_change [value >= 0 & value < 5 ] <- "2.0-5"
  r_sum_change [value >= 5 & value < 10 ] <- "3.5-10"
  r_sum_change [value >= 10 & value < 15 ] <- "4.10-15"
  r_sum_change [value >= 15 & value < 20 ] <- "5.15-20"
  r_sum_change [value >= 20 & value <25 ] <- "6.20-25"
  r_sum_change [value >= 25 & value < 30 ] <- "7.25-30"
  r_sum_change [value >= 30 ] <- "8.>=30"
  
})

 scr3_1_change$Sgg_cd <- as.character(scr3_1_change$Sgg_cd)
 
dat4 <- shap2016 %>% 
  mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
  left_join(.,scr3_1_change, by = 'Sgg_cd') %>% 
  na.omit(Sgg_cd) 

s3 <- ggplot(data = dat4) +
  geom_sf(data = dat4, aes(fill = r_sum_change)) +
  scale_fill_brewer(name =" Total attendance rate change (%)", palette = "YlGnBu",
                    guide = guide_legend(override.aes = list(linetype = "blank", shape = NA)))+
  facet_wrap(.~ period) +
  theme(text = element_text(size = 60))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 50),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", color = NA)) +
  labs( title = "All cancer rate change") +
  theme(legend.position = "bottom",
        legend.text = element_text(size=50),
        legend.title = element_text(size=50))+
  north(shap2016, symbol = 3)

png(filename = "map_change.png", 
    units = "in", 
    width = 25, 
    height = 60, 
    res = 300)

ggpubr::ggarrange(map_breast_change, map_cer_change, map_colon_change, map_liver_change,map_sto_change, s3, nrow = 6)

dev.off()



