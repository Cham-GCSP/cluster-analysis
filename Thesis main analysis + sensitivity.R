
setwd(getwd())
#_________________ Mortality Cluster using flexible scan________________
library(pacman)
pacman::p_load(sp,data.table,plyr,readr,sf,rflexscan,dplyr,ggplot2,tidyverse,spData,spdep,
               RColorBrewer,smerc)

# Load shapefile

shap2010 <- st_read('Data/Data/Sgg_2010.shp')

shap2010_sp <- as(shap2010, 'Spatial')
#__________________________________________________________________

#Check neighbor

nb   <- poly2nb(shap2010_sp, queen = T)

nb1  <- poly2nb(shap2010_sp, queen = F)

isTRUE(all.equal(nb, nb1, check.attributes=T))
isTRUE(all.equal(nb, nb1, check.attributes=F))


oopar <- par(mfrow=c(1,2), mar=c(3,3,1,1)+0.1)

par(mar=c(1, 1, 1, 1))

plot(shap2010_sp, border="grey60")

plot(nb, coordinates(shap2010_sp), add=TRUE, pch=19, cex=0.6)

text(bbox(shap2010_sp)[1,1], bbox(shap2010_sp)[2,2], labels="a)", cex=0.8)

plot(shap2010_sp, border="grey60")

plot(nb, coordinates(shap2010_sp), add=TRUE, pch=19, cex=0.6)

plot(diffnb(nb, nb1, verbose=FALSE), coordinates(shap2010_sp),
     add=TRUE, pch=".", cex=0.6, lwd=4,col="red")

text(bbox(shap2010_sp)[1,1], bbox(shap2010_sp)[2,2], labels="b)", cex=0.8)

#__________________________________________________________________

#Import csv file

mor_cluster <- read.csv('Data/Data/mor_data.csv')

mor_cluster$Sgg_cd <- as.character(mor_cluster$Sgg_cd)

Sgg_change <- read.csv('Data/Data/Sgg_change.csv')

Sgg_change <- subset(Sgg_change, select = -c(5,6))

shap2010_sub <- shap2010[-c(29,59,201,202,227,233,244),]

code_change1 <- Sgg_change %>% 
  filter(from_year > 1999 & from_year < 2014) %>% 
  dplyr::select(fromcode, tocode)

#_______________________________________________________________________

# Separate data for each period
coord_matrix <- as.matrix(st_coordinates(st_centroid(shap2010_sub)))

coord <- as.data.frame(coord_matrix)

write.csv(coord, 'coord.csv')



mor1 <- mor_cluster %>%
  mutate(Sgg_cd1 = plyr::mapvalues(Sgg_cd, code_change1$fromcode, code_change1$tocode)) %>% 
  filter(year_agg == "1999-2003") %>% 
  select(.,-c(1,2,6,8)) %>% 
  group_by(gender,cause,Sgg_cd1) %>% 
  summarize(D = sum(D),Population = sum(Population)) %>% 
  pivot_wider(.,id_cols = c('Sgg_cd1'), names_from = c('gender','cause'), values_from = c('D','Population'), names_sep = '-' )


mor1_merge <- shap2010_sub %>% 
  mutate(Sgg_cd1 = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
  left_join(., mor1, by = 'Sgg_cd1')

mor1 <- st_as_sf(mor1_merge, coords  = c("x","y"))

st_write(mor1, 'mor_scan1.csv')


mor2 <- mor_cluster %>%
  mutate(Sgg_cd1 = plyr::mapvalues(Sgg_cd, code_change1$fromcode, code_change1$tocode)) %>% 
  filter(year_agg == "2004-2008") %>% 
  select(.,-c(1,2,6,8)) %>% 
  group_by(gender,cause,Sgg_cd1) %>% 
  summarize(D = sum(D),Population = sum(Population)) %>% 
  pivot_wider(.,id_cols = c('Sgg_cd1'), names_from = c('gender','cause'), values_from = c('D','Population'), names_sep = '-' )

mor2_merge <- shap2010_sub %>% 
  mutate(Sgg_cd1 = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
  left_join(., mor2, by = 'Sgg_cd1')

mor2 <- st_as_sf(mor2_merge, coords  = c("x","y"))

st_write(mor2, 'mor_scan2.csv')



mor3 <- mor_cluster %>%
  mutate(Sgg_cd1 = plyr::mapvalues(Sgg_cd, code_change1$fromcode, code_change1$tocode)) %>% 
  filter(year_agg == "2009-2013") %>% 
  select(.,-c(1,2,6,8)) %>% 
  group_by(gender,cause,Sgg_cd1) %>% 
  summarize(D = sum(D),Population = sum(Population)) %>% 
  pivot_wider(.,id_cols = c('Sgg_cd1'), names_from = c('gender','cause'), values_from = c('D','Population'), names_sep = '-' )


mor3_merge <- shap2010_sub %>% 
  mutate(Sgg_cd1 = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
  left_join(., mor3, by = 'Sgg_cd1')

mor3 <- st_as_sf(mor3_merge, coords  = c("x","y"))

st_write(mor3, 'mor_scan3.csv')

# flexible scan, remove areas without connection (7 district)

flexrun = function(g,c,y,k){
  
  data = shap2010_sub %>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(mor_cluster, by = 'Sgg_cd') %>% 
    filter(gender == g & cause == c & year_agg == y ) 
  
  names_vec    <- data$Sgg_cd
  
  adj_matrix   <- nb2mat(poly2nb(data), style = 'B', zero.policy = TRUE)
  
  coord_matrix <- as.matrix(st_coordinates(st_centroid(data)))
  
  case_max  <- data %>% st_drop_geometry 
  
  case_matrix <- as.matrix(case_max %>%.[,'D'])
  
  pmat <- data$Population
    
  expected <- sum(as.vector(case_matrix))/sum(pmat)

  rflexscan(x       = coord_matrix[,1],
            y       = coord_matrix[,2],
            name    = names_vec,
            observed = as.vector(case_matrix),
            expected = expected*pmat,
            nb         = adj_matrix,
            stattype   = 'RESTRICTED',
            scanmethod = 'FLEXIBLE',
            ralpha     = 0.05, 
            simcount   = 999,
            rantype    = 'POISSON',
            comments   = "", verbose = FALSE, secondary = NULL,
            clustersize  = k)
  
}
  

# Main analysis 

# n= 37

Male_lung_mor1_37  <- flexrun('Male', 'lung','1999-2003', 37)

Female_lung_mor1_37 <- flexrun('Female', 'lung' , '1999-2003', 37)


# Stomach

Male_stomach_mor1_37  <- flexrun('Male', 'stomach', '1999-2003',37)

Female_stomach_mor1_37 <- flexrun( 'Female', 'stomach', '1999-2003',37)


# Uterus
Female_uterus_mor1_37  <- flexrun('Female', 'uterus', '1999-2003',37)


# Breast
Female_breast_mor1_37  <- flexrun('Female', 'breast', '1999-2003',37)


# Colon

Male_colon_mor1_37   <- flexrun( 'Male', 'colon', '1999-2003',37)

Female_colon_mor1_37 <- flexrun( 'Female', 'colon', '1999-2003',37)

Male_liver_mor1_37  <- flexrun('Male', 'liver', '1999-2003',37)

Female_liver_mor1_37 <- flexrun('Female', 'liver', '1999-2003',37)


# Prostate
Male_prostate_mor1_37  <- flexrun('Male', 'prostate', '1999-2003',37)


#2004-2008

#Lung cancer

Male_lung_mor2_37  <- flexrun('Male', 'lung', '2004-2008',37)

Female_lung_mor2_37 <- flexrun('Female', 'lung', '2004-2008',37)


# Stomach

Male_stomach_mor2_37  <- flexrun('Male', 'stomach', '2004-2008',37)

Female_stomach_mor2_37 <- flexrun( 'Female', 'stomach', '2004-2008',37)


# Uterus
Female_uterus_mor2_37  <- flexrun('Female', 'uterus', '2004-2008',37)


# Breast
Female_breast_mor2_37  <- flexrun('Female', 'breast', '2004-2008',37)

# Colon

Male_colon_mor2_37   <- flexrun( 'Male', 'colon', '2004-2008',37)

Female_colon_mor2_37 <- flexrun( 'Female', 'colon', '2004-2008',37)

#liver cancer

Male_liver_mor2_37   <- flexrun('Male', 'liver', '2004-2008',37)

Female_liver_mor2_37 <- flexrun('Female', 'liver', '2004-2008',37)


# Prostate
Male_prostate_mor2_37  <- flexrun('Male', 'prostate', '2004-2008',37)

#_______________________________________________________________
#2009-2013

#Lung cancer

Male_lung_mor3_37  <- flexrun('Male', 'lung', '2009-2013',37)

Female_lung_mor3_37 <- flexrun('Female', 'lung', '2009-2013',37)


# Stomach

Male_stomach_mor3_37  <- flexrun('Male', 'stomach', '2009-2013',37)

Female_stomach_mor3_37 <- flexrun( 'Female', 'stomach', '2009-2013',37)

# Uterus
Female_uterus_mor3_37  <- flexrun('Female', 'uterus', '2009-2013',37)


# Breast
Female_breast_mor3_37  <- flexrun('Female', 'breast', '2009-2013',37)

# Colon

Male_colon_mor3_37   <- flexrun( 'Male', 'colon', '2009-2013',37)

Female_colon_mor3_37 <- flexrun( 'Female', 'colon', '2009-2013',37)

#liver cancer

Male_liver_mor3_37   <- flexrun('Male', 'liver', '2009-2013',37)

Female_liver_mor3_37 <- flexrun('Female', 'liver', '2009-2013',37)

# Prostate
Male_prostate_mor3_37  <- flexrun('Male', 'prostate', '2009-2013',37)

#Plot cluster
fun_color_range <- colorRampPalette(c("#1b98e0", "red"))


# Incidence

cluster_inc <- read.csv('Data/Data/inc_data.csv')

cluster_inc$Sgg_cd <- as.character(cluster_inc$Sgg_cd)

# export incidence for Satscan

inc1 <- cluster_inc %>%
  mutate(Sgg_cd1 = plyr::mapvalues(Sgg_cd, code_change1$fromcode, code_change1$tocode)) %>% 
  filter(year_agg == "1999-2003") %>% 
  select(.,-c(1,2,5,6,8,9)) %>% 
  group_by(gender,cause,Sgg_cd1) %>% 
  summarize(counts = sum(counts),Population = sum(Population)) %>% 
  pivot_wider(.,id_cols = c('Sgg_cd1'), names_from = c('gender','cause'), values_from = c('counts','Population'), names_sep = '-' )


inc1_merge <- shap2010_sub %>% 
  mutate(Sgg_cd1 = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
  left_join(., inc1, by = 'Sgg_cd1')

inc1 <- st_as_sf(inc1_merge, coords  = c("x","y"))

st_write(inc1, 'inc_scan1.csv')


inc2 <- cluster_inc %>%
  mutate(Sgg_cd1 = plyr::mapvalues(Sgg_cd, code_change1$fromcode, code_change1$tocode)) %>% 
  filter(year_agg == "2004-2008") %>% 
  select(.,-c(1,2,6,8,9)) %>% 
  group_by(gender,cause,Sgg_cd1) %>% 
  summarize(counts = sum(counts),Population = sum(Population)) %>% 
  pivot_wider(.,id_cols = c('Sgg_cd1'), names_from = c('gender','cause'), values_from = c('counts','Population'), names_sep = '-' )


inc2_merge <- shap2010_sub %>% 
  mutate(Sgg_cd1 = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
  left_join(., inc2, by = 'Sgg_cd1')

inc2 <- st_as_sf(inc2_merge, coords  = c("x","y"))

st_write(inc2, 'inc_scan2.csv')


inc3 <- cluster_inc %>%
  mutate(Sgg_cd1 = plyr::mapvalues(Sgg_cd, code_change1$fromcode, code_change1$tocode)) %>% 
  filter(year_agg == "2009-2013") %>% 
  select(.,-c(1,2,6,8,9)) %>% 
  group_by(gender,cause,Sgg_cd1) %>% 
  summarize(counts = sum(counts),Population = sum(Population)) %>% 
  pivot_wider(.,id_cols = c('Sgg_cd1'), names_from = c('gender','cause'), values_from = c('counts','Population'), names_sep = '-' )


inc3_merge <- shap2010_sub %>% 
  mutate(Sgg_cd1 = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
  left_join(., inc3, by = 'Sgg_cd1')

inc3 <- st_as_sf(inc3_merge, coords  = c("x","y"))

st_write(inc3, 'inc_scan3.csv')

# function

flexrun_inc0 = function(g,c,y,k){
  
  data <- shap2010_sub%>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(., cluster_inc, by = 'Sgg_cd') %>% 
    filter(gender == g & cause == c & year_agg == y) 
  
  names_vec    <- data$Sgg_cd
  
  adj_matrix   <- nb2mat(poly2nb(data), style = 'B', zero.policy = TRUE)
  
  coord_matrix <- as.matrix(st_coordinates(st_centroid(data)))
  
  case_max  <- data %>% st_drop_geometry
  
  case_matrix <- as.matrix(case_max %>%.[,'counts'])
  
  pmat <- data$Population
  
  expected <- sum(as.vector(case_matrix)/sum(pmat))
  
  rflexscan(x       = coord_matrix[,1],
            y       = coord_matrix[,2],
            name    = names_vec,
            observed = as.vector(case_matrix),
            expected = expected*pmat,
            nb         = adj_matrix,
            stattype   = 'RESTRICTED',
            scanmethod = 'FLEXIBLE',
            ralpha     = 0.05, 
            simcount   = 999,
            rantype    = 'POISSON',
            comments   = "", verbose = FALSE, secondary = NULL,
            clustersize = k)
  
}

# ________________________________________________-
#Main analysis for size = 15%

# run function

Male_lung_inc1_37  <- flexrun_inc0('Male', 'lung', '1999-2003',37)

Female_lung_inc1_37   <- flexrun_inc0('Female', 'lung', '1999-2003',37)


# Stomach

Male_stomach_inc1_37    <- flexrun_inc0('Male', 'stomach', '1999-2003',37)

Female_stomach_inc1_37   <- flexrun_inc0( 'Female', 'stomach', '1999-2003',37)


# Uterus
Female_uterus_inc1_37    <- flexrun_inc0('Female', 'uterus', '1999-2003',37)


# Breast
Female_breast_inc1_37    <- flexrun_inc0('Female', 'breast', '1999-2003',37)


# Colon

Male_colon_inc1_37     <- flexrun_inc0( 'Male', 'colon', '1999-2003',37)

Female_colon_inc1_37   <- flexrun_inc0( 'Female', 'colon', '1999-2003',37)


#liver cancer


Male_liver_inc1_37     <- flexrun_inc0('Male', 'liver', '1999-2003',37)

Female_liver_inc1_37   <- flexrun_inc0('Female', 'liver', '1999-2003',37)


# Prostate
Male_prostate_inc1_37    <- flexrun_inc0('Male', 'prostate', '1999-2003',37)
# thyroid

Male_thyroid_inc1_37     <- flexrun_inc0('Male', 'thyroid', '1999-2003',37)

Female_thyroid_inc1_37   <- flexrun_inc0('Female', 'thyroid', '1999-2003',37)


#----------------------------------------------------------------------

#2004-2008

#Lung cancer

Male_lung_inc2_37    <- flexrun_inc0('Male', 'lung', '2004-2008',37)

Female_lung_inc2_37   <- flexrun_inc0('Female', 'lung', '2004-2008',37)


# Stomach

Male_stomach_inc2_37    <- flexrun_inc0('Male', 'stomach', '2004-2008',37)

Female_stomach_inc2_37   <- flexrun_inc0( 'Female', 'stomach', '2004-2008',37)


# Uterus
Female_uterus_inc2_37    <- flexrun_inc0('Female', 'uterus', '2004-2008',37)


# Breast
Female_breast_inc2_37    <- flexrun_inc0('Female', 'breast', '2004-2008',37)


# Colon

Male_colon_inc2_37     <- flexrun_inc0( 'Male', 'colon', '2004-2008',37)

Female_colon_inc2_37   <- flexrun_inc0( 'Female', 'colon', '2004-2008',37)


#liver cancer

Male_liver_inc2_37     <- flexrun_inc0('Male', 'liver', '2004-2008',37)

Female_liver_inc2_37   <- flexrun_inc0('Female', 'liver', '2004-2008',37)


# Prostate
Male_prostate_inc2_37    <- flexrun_inc0('Male', 'prostate', '2004-2008',37)

# thyroid

Male_thyroid_inc2_37     <- flexrun_inc0('Male', 'thyroid', '2004-2008',37)

Female_thyroid_inc2_37   <- flexrun_inc0('Female', 'thyroid', '2004-2008',37)
#_______________________________________________________________
#2009-2013

#Lung cancer

Male_lung_inc3_37    <- flexrun_inc0('Male', 'lung', '2009-2013',37)

Female_lung_inc3_37   <- flexrun_inc0('Female', 'lung', '2009-2013',37)


# Stomach

Male_stomach_inc3_37    <- flexrun_inc0('Male', 'stomach', '2009-2013',37)

Female_stomach_inc3_37   <- flexrun_inc0( 'Female', 'stomach', '2009-2013',37)


# Uterus
Female_uterus_inc3_37    <- flexrun_inc0('Female', 'uterus', '2009-2013',37)


# Breast
Female_breast_inc3_37    <- flexrun_inc0('Female', 'breast', '2009-2013',37)


# Colon

Male_colon_inc3_37     <- flexrun_inc0( 'Male', 'colon', '2009-2013',37)

Female_colon_inc3_37   <- flexrun_inc0( 'Female', 'colon', '2009-2013',37)


#liver cancer


Male_liver_inc3_37     <- flexrun_inc0('Male', 'liver', '2009-2013',37)

Female_liver_inc3_37  <- flexrun_inc0('Female', 'liver', '2009-2013',37)


# Prostate
Male_prostate_inc3_37    <- flexrun_inc0('Male', 'prostate', '2009-2013',37)

# thyroid 

Male_thyroid_inc3_37     <- flexrun_inc0('Male', 'thyroid', '2009-2013',37)

Female_thyroid_inc3_37  <- flexrun_inc0('Female', 'thyroid', '2009-2013',37)


# color function

cols<- brewer.pal(3, "PuBuGn")

pal <- colorRampPalette(cols)

cols1<- brewer.pal(3, "Blues")

pal1 <- colorRampPalette(cols1)

fun_color_range1 <- colorRampPalette(c("#8DD3C7", "orange"))

fun_color_range <- colorRampPalette(c("#1b98e0", "red"))

#__________________________ plot main results_______________________--

# plot main results different type of cancer incidence & mortality

png(filename = "Inc_male_cluster.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(6,3))

choropleth(shap2010_sub$geometry, Male_lung_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_lung_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_lung_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_lung_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_stomach_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_stomach_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_stomach_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Male_colon_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_colon_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_colon_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_colon_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_liver_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_liver_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_liver_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_prostate_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Prostate cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_prostate_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Prostate cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_prostate_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Prostate cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_thyroid_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_thyroid_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_thyroid_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 3', cex.sub = 7)
dev.off()

# female incidence


png(filename = "Inc_female_cluster.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(7,3))

choropleth(shap2010_sub$geometry, Female_lung_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_lung_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_lung_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_lung_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_stomach_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_stomach_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_stomach_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Female_colon_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_colon_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_colon_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_colon_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_liver_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_liver_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_liver_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_breast_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Breast cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_breast_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Breast cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_breast_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Breast cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_uterus_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Uterus cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_uterus_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Uterus cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_uterus_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Uterus cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_thyroid_inc1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_thyroid_inc2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_thyroid_inc3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 3', cex.sub = 7)
dev.off()


# mortality 

png(filename = "Mor_male_cluster.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(5,3))


choropleth(shap2010_sub$geometry, Male_lung_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_lung_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_lung_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_lung_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Male_stomach_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_stomach_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_stomach_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_stomach_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_colon_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_colon_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_colon_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_liver_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_liver_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_liver_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_prostate_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Prostate cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_prostate_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Prostate cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_prostate_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Prostate cancer - Period 3', cex.sub = 7)
dev.off()


# female mortality

png(filename = "Mor_female_cluster.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(6,3))

choropleth(shap2010_sub$geometry, Female_lung_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_lung_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_lung_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_lung_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_stomach_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_stomach_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_stomach_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Female_colon_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_colon_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_colon_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_colon_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_liver_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_liver_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_liver_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Female_breast_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_breast_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Breast cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_breast_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Breast cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_breast_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Breast cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_uterus_mor1_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Uterus cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_uterus_mor2_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Uterus cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_uterus_mor3_37$cluster[[1]],
           pval = 0.05)
title(sub = 'Uterus cancer - Period 3', cex.sub = 7)

dev.off()

# check overlap 

png(filename =  " overlap_male_per1.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(5,2))

choropleth(shap2010_sub$geometry, Male_lung_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_lung_inc1_37$cluster[[1]],
           pval = 0.05)
title('Lung cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_lung_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_lung_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Lung cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_stomach_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_stomach_inc1_37$cluster[[1]],
           pval = 0.05)
title('Stomach cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_stomach_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_stomach_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Stomach cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_colon_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_colon_inc1_37$cluster[[1]],
           pval = 0.05)
title('Colorectal cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_colon_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_colon_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Colorectal cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_liver_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_liver_inc1_37$cluster[[1]],
           pval = 0.05)
title('Liver cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_liver_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_liver_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Liver cancer - Period 1', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_prostate_inc1_37$cluster[[1]],
           pval = 0.05)
title('Prostate cancer - Period 1', line = -5, cex.main = 7)
dev.off()
choropleth(shap2010_sub$geometry, Male_prostate_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_prostate_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Prostate cancer - Period 1', line = -5, cex.main = 7)
dev.off()

# period 2

png(filename =  " overlap_male_per2.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(5,2))

choropleth(shap2010_sub$geometry, Male_lung_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_lung_inc2_37$cluster[[1]],
           pval = 0.05)
title('Lung cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_lung_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_lung_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Lung cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_stomach_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_stomach_inc2_37$cluster[[1]],
           pval = 0.05)
title('Stomach cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_stomach_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_stomach_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Stomach cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_colon_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_colon_inc2_37$cluster[[1]],
           pval = 0.05)
title('Colorectal cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_colon_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_colon_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Colorectal cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_liver_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_liver_inc2_37$cluster[[1]],
           pval = 0.05)
title('Liver cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_liver_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_liver_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Liver cancer - Period 2', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_prostate_inc2_37$cluster[[1]],
           pval = 0.05)
title('Prostate cancer - Period 2', line = -5, cex.main = 7)
dev.off()
choropleth(shap2010_sub$geometry, Male_prostate_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_prostate_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Prostate cancer - Period 2', line = -5, cex.main = 7)
dev.off()

# period 3

png(filename =  " overlap_male_per3.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(5,2))

choropleth(shap2010_sub$geometry, Male_lung_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_lung_inc3_37$cluster[[1]],
           pval = 0.05)
title('Lung cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_lung_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_lung_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Lung cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_stomach_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_stomach_inc3_37$cluster[[1]],
           pval = 0.05)
title('Stomach cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_stomach_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_stomach_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Stomach cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_colon_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_colon_inc3_37$cluster[[1]],
           pval = 0.05)
title('Colorectal cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_colon_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_colon_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Colorectal cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_liver_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_liver_inc3_37$cluster[[1]],
           pval = 0.05)
title('Liver cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Male_liver_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_liver_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Liver cancer - Period 3', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Male_prostate_inc3_37$cluster[[1]],
           pval = 0.05)
title('Prostate cancer - Period 3', line = -5, cex.main = 7)
dev.off()
choropleth(shap2010_sub$geometry, Male_prostate_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Male_prostate_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Prostate cancer - Period 3', line = -5, cex.main = 7)
dev.off()

# overlap in female
# Female

png(filename =  " overlap_Female_per1.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(6,2))
choropleth(shap2010_sub$geometry, Female_lung_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_lung_inc1_37$cluster[[1]],
           pval = 0.05)
title('Lung cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_lung_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_lung_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Lung cancer - Period 1', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_stomach_inc1_37$cluster[[1]],
           pval = 0.05)
title('Stomach cancer - Period 1', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_stomach_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Stomach cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_colon_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_colon_inc1_37$cluster[[1]],
           pval = 0.05)
title('Colorectal cancer - Period 1', line = -5, cex.main = 7)


choropleth(shap2010_sub$geometry, Female_colon_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_colon_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Colorectal cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_liver_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_liver_inc1_37$cluster[[1]],
           pval = 0.05)
title('Liver cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_liver_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_liver_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Liver cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_breast_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_breast_inc1_37$cluster[[1]],
           pval = 0.05)
title('Breast cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_breast_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_breast_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Breast cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_uterus_inc1_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_uterus_inc1_37$cluster[[1]],
           pval = 0.05)
title('Uterus cancer - Period 1', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_uterus_mor1_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_uterus_mor1_37$cluster[[1]],
           pval = 0.05)
title( 'Uterus cancer - Period 1', line = -5, cex.main = 7)
dev.off()

# period 2

png(filename =  " overlap_Female_per2.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(6,2))
choropleth(shap2010_sub$geometry, Female_lung_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_lung_inc2_37$cluster[[1]],
           pval = 0.05)
title('Lung cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_lung_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_lung_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Lung cancer - Period 2', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_stomach_inc2_37$cluster[[1]],
           pval = 0.05)
title('Stomach cancer - Period 2', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_stomach_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Stomach cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_colon_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_colon_inc2_37$cluster[[1]],
           pval = 0.05)
title('Colorectal cancer - Period 2', line = -5, cex.main = 7)


choropleth(shap2010_sub$geometry, Female_colon_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_colon_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Colorectal cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_liver_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_liver_inc2_37$cluster[[1]],
           pval = 0.05)
title('Liver cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_liver_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_liver_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Liver cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_breast_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_breast_inc2_37$cluster[[1]],
           pval = 0.05)
title('Breast cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_breast_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_breast_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Breast cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_uterus_inc2_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_uterus_inc2_37$cluster[[1]],
           pval = 0.05)
title('Uterus cancer - Period 2', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_uterus_mor2_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_uterus_mor2_37$cluster[[1]],
           pval = 0.05)
title( 'Uterus cancer - Period 2', line = -5, cex.main = 7)
dev.off()

# period 3

png(filename =  " overlap_Female_per2.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(6,2))
choropleth(shap2010_sub$geometry, Female_lung_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_lung_inc3_37$cluster[[1]],
           pval = 0.05)
title('Lung cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_lung_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_lung_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Lung cancer - Period 3', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_stomach_inc3_37$cluster[[1]],
           pval = 0.05)
title('Stomach cancer - Period 3', line = -5, cex.main = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_stomach_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Stomach cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_colon_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_colon_inc3_37$cluster[[1]],
           pval = 0.05)
title('Colorectal cancer - Period 3', line = -5, cex.main = 7)


choropleth(shap2010_sub$geometry, Female_colon_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_colon_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Colorectal cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_liver_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_liver_inc3_37$cluster[[1]],
           pval = 0.05)
title('Liver cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_liver_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_liver_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Liver cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_breast_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_breast_inc3_37$cluster[[1]],
           pval = 0.05)
title('Breast cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_breast_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_breast_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Breast cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_uterus_inc3_37,
           col = fun_color_range1(1),
           region_color = 'white',
           rank = Female_uterus_inc3_37$cluster[[1]],
           pval = 0.05)
title('Uterus cancer - Period 3', line = -5, cex.main = 7)
choropleth(shap2010_sub$geometry, Female_uterus_mor3_37,
           col = fun_color_range(1),
           region_color = 'white',
           rank = Female_uterus_mor3_37$cluster[[1]],
           pval = 0.05)
title( 'Uterus cancer - Period 3', line = -5, cex.main = 7)
dev.off()

#summary statistic of cluster

Male_colon_mor1_37$cluster[[1]]
summary(Male_colon_mor1_37)

Female_colon_mor1_37$cluster[[1]]
summary(Female_colon_mor1_37)


Male_liver_mor1_37$cluster[[1]]
summary(Male_liver_mor1_37)

Female_liver_mor1_37$cluster[[1]]
summary(Female_liver_mor1_37)


Male_stomach_mor1_37$cluster[[1]]
summary(Male_stomach_mor1_37)

Female_stomach_mor1_37$cluster[[1]]
summary(Female_stomach_mor1_37)


Male_lung_mor1_37$cluster[[1]]
summary(Male_lung_mor1_37)

Female_lung_mor1_37$cluster[[1]]
summary(Female_lung_mor1_37)

Female_breast_mor1_37$cluster[[1]]
summary(Female_breast_mor1_37)

Female_uterus_mor1_37$cluster[[1]]
summary(Female_uterus_mor1_37)

Male_prostate_mor1_37$cluster[[1]]
summary(Male_prostate_mor1_37)

# period 2
Male_colon_mor2_37$cluster[[1]]
summary(Male_colon_mor2_37)

Female_colon_mor2_37$cluster[[1]]
summary(Female_colon_mor2_37)


Male_liver_mor2_37$cluster[[1]]
summary(Male_liver_mor2_37)

Female_liver_mor2_37$cluster[[1]]
summary(Female_liver_mor2_37)


Male_stomach_mor2_37$cluster[[1]]
summary(Male_stomach_mor2_37)

Female_stomach_mor2_37$cluster[[1]]
summary(Female_stomach_mor2_37)


Male_lung_mor2_37$cluster[[1]]
summary(Male_lung_mor2_37)

Female_lung_mor2_37$cluster[[1]]
summary(Female_lung_mor2_37)

Female_breast_mor2_37$cluster[[1]]
summary(Female_breast_mor2_37)

Female_uterus_mor2_37$cluster[[1]]
summary(Female_uterus_mor2_37)

Male_prostate_mor2_37$cluster[[1]]
summary(Male_prostate_mor2_37)

# period 3

Male_colon_mor3_37$cluster[[1]]
summary(Male_colon_mor3_37)

Female_colon_mor3_37$cluster[[1]]
summary(Female_colon_mor3_37)


Male_liver_mor3_37$cluster[[1]]
summary(Male_liver_mor3_37)

Female_liver_mor3_37$cluster[[1]]
summary(Female_liver_mor3_37)


Male_stomach_mor3_37$cluster[[1]]
summary(Male_stomach_mor3_37)

Female_stomach_mor3_37$cluster[[1]]
summary(Female_stomach_mor3_37)


Male_lung_mor3_37$cluster[[1]]
summary(Male_lung_mor3_37)

Female_lung_mor3_37$cluster[[1]]
summary(Female_lung_mor3_37)

Female_breast_mor3_37$cluster[[1]]
summary(Female_breast_mor3_37)

Female_uterus_mor3_37$cluster[[1]]
summary(Female_uterus_mor3_37)

Male_prostate_mor3_37$cluster[[1]]
summary(Male_prostate_mor3_37)

# incidence

Male_colon_inc1_37$cluster[[1]]
summary(Male_colon_inc1_37)

Female_colon_inc1_37$cluster[[1]]
summary(Female_colon_inc1_37)


Male_liver_inc1_37$cluster[[1]]
summary(Male_liver_inc1_37)

Female_liver_inc1_37$cluster[[1]]
summary(Female_liver_inc1_37)


Male_stomach_inc1_37$cluster[[1]]
summary(Male_stomach_inc1_37)

Female_stomach_inc1_37$cluster[[1]]
summary(Female_stomach_inc1_37)


Male_lung_inc1_37$cluster[[1]]
summary(Male_lung_inc1_37)

Female_lung_inc1_37$cluster[[1]]
summary(Female_lung_inc1_37)

Female_breast_inc1_37$cluster[[1]]
summary(Female_breast_inc1_37)

Female_uterus_inc1_37$cluster[[1]]
summary(Female_uterus_inc1_37)

Male_prostate_inc1_37$cluster[[1]]
summary(Male_prostate_inc1_37)

Male_thyroid_inc1_37$cluster[[1]]
summary(Male_thyroid_inc1_37)

Female_thyroid_inc1_37$cluster[[1]]
summary(Female_thyroid_inc1_37)

# period 2
Male_colon_inc2_37$cluster[[1]]
summary(Male_colon_inc2_37)

Female_colon_inc2_37$cluster[[1]]
summary(Female_colon_inc2_37)


Male_liver_inc2_37$cluster[[1]]
summary(Male_liver_inc2_37)

Female_liver_inc2_37$cluster[[1]]
summary(Female_liver_inc2_37)


Male_stomach_inc2_37$cluster[[1]]
summary(Male_stomach_inc2_37)

Female_stomach_inc2_37$cluster[[1]]
summary(Female_stomach_inc2_37)


Male_lung_inc2_37$cluster[[1]]
summary(Male_lung_inc2_37)

Female_lung_inc2_37$cluster[[1]]
summary(Female_lung_inc2_37)

Female_breast_inc2_37$cluster[[1]]
summary(Female_breast_inc2_37)

Female_uterus_inc2_37$cluster[[1]]
summary(Female_uterus_inc2_37)

Male_prostate_inc2_37$cluster[[1]]
summary(Male_prostate_inc2_37)
Male_thyroid_inc2_37$cluster[[1]]
summary(Male_thyroid_inc2_37)

Female_thyroid_inc2_37$cluster[[1]]
summary(Female_thyroid_inc2_37)


# period 3

Male_colon_inc3_37$cluster[[1]]
summary(Male_colon_inc3_37)

Female_colon_inc3_37$cluster[[1]]
summary(Female_colon_inc3_37)


Male_liver_inc3_37$cluster[[1]]
summary(Male_liver_inc3_37)

Female_liver_inc3_37$cluster[[1]]
summary(Female_liver_inc3_37)


Male_stomach_inc3_37$cluster[[1]]
summary(Male_stomach_inc3_37)

Female_stomach_inc3_37$cluster[[1]]
summary(Female_stomach_inc3_37)


Male_lung_inc3_37$cluster[[1]]
summary(Male_lung_inc3_37)

Female_lung_inc3_37$cluster[[1]]
summary(Female_lung_inc3_37)

Female_breast_inc3_37$cluster[[1]]
summary(Female_breast_inc3_37)

Female_uterus_inc3_37$cluster[[1]]
summary(Female_uterus_inc3_37)

Male_prostate_inc3_37$cluster[[1]]
summary(Male_prostate_inc3_37)

Male_thyroid_inc3_37$cluster[[1]]
summary(Male_thyroid_inc3_37)

Female_thyroid_inc3_37$cluster[[1]]
summary(Female_thyroid_inc3_37)


#____________plot all cluster_________________

png(filename = "Inc_male_cluster_all.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(6,3))

choropleth(shap2010_sub$geometry, Male_lung_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Male_colon_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_thyroid_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_thyroid_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_thyroid_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 3', cex.sub = 7)
dev.off()

# female incidence


png(filename = "Inc_female_cluster_all.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(7,3))

choropleth(shap2010_sub$geometry, Female_lung_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Female_colon_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc1_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_thyroid_inc1_37$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc2_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_thyroid_inc2_37$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc3_37,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_thyroid_inc3_37$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 3', cex.sub = 7)
dev.off()


# mortality 

png(filename = "Mor_male_cluster_all.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(5,3))


choropleth(shap2010_sub$geometry, Male_lung_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Male_stomach_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 3', cex.sub = 7)
dev.off()


# female mortality

png(filename = "Mor_female_cluster_all.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(6,3))

choropleth(shap2010_sub$geometry, Female_lung_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Female_colon_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Female_breast_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_mor1_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_mor1_37$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_mor2_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_mor2_37$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_mor3_37,
           col = rev(pal(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_mor3_37$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 3', cex.sub = 7)

dev.off()

#_______________t-test for difference of characteristics of MLC with other areas____________

factors <- read.csv('Data/Data/factors.csv')
factors$sgg_cd <- as.character(factors$sgg_cd)

data_factor = shap2010_sub %>% 
  left_join(factors, by = c('SIGUNGU_CD' = 'sgg_cd')) %>% 
  st_drop_geometry() %>% 
  select(.,-c(3,4)) %>% 
  mutate(m_lung_m = if_else(SIGUNGU_CD %in% c('36370', '36380', '36390', '36400', '36410', '36420', '36430', '36440', '36450', '36480', '38060', '38340', '38360', '38370', '38380'),'1','0'),
         m_lung_f = if_else(SIGUNGU_CD  %in%  c('33340', '35330', '37020', '37030', '37060', '37070', '37080', '37090', '37310', '37320', '37330', '37350', 
                                            '37360', '37370', '37380', '37400', '37410', '38080', '38330', '38390', '38400'),'1','0'),
         m_stomach_m = if_else(SIGUNGU_CD  %in%  c('33320', '33330', '33340', '34310', '35320', '35330', '35340', '37030', '37070', '37080', '37310', '37320', 
                                               '37360', '37370', '37380', '38080', '38310', '38320', '38330', '38370', '38380', '38390', '38400'),'1','0'),
         m_stomach_f = if_else(SIGUNGU_CD  %in%  c('33330', '33340', '34310', '35320', '35330', '35340', '37030', '37070', '37080','37310', '37320', '37360', 
                                               '37370', '37380', '38080', '38310', '38320', '38330', '38370', '38380', '38390 ','38400'),'1','0'),
         m_colon_m = if_else(SIGUNGU_CD  %in%  c('33320', '33330', '33340', '34020', '34030', '34060', '34070', '34310', '34320', '34330', '34340', '34350', '35040', '35050', '35060', '35310', '35330', '35340' , '35350',  '35360',  '35370',  '35380', '36310',  '36450', '38380', '38390'),'1','0'),
         m_colon_f = if_else(SIGUNGU_CD  %in%  c('32040','32070','32320', '32330', '32340', '33020', '33310', '33320', '33360', '33370', '37030', '37040', 
                                             '37060', '37070', '37080', '37320', '37330', '37340', '37350', '37400','37410','37420'),'1','0'),
         m_liver_m = if_else(SIGUNGU_CD  %in%  c('35350', '35360', '35370', '36040', '36310', '36320', '36330', '36350', '36360', '36370', '36380', '36390', 
                                             '36400', '36410', '36420', '36430', '36440', '36450', '36480', '38030', '38050', '38060', '38340', '38360', 
                                             '38370', '38380'),'1','0'),
         m_liver_f = if_else(SIGUNGU_CD  %in%  c('35320', '35330', '35340', '35350', '35360', '36320', '36330', '37030', '37370', '37380', '38050', '38060', 
                                             '38080', '38310', '38320','38330', '38340', '38360', '38370', '38380', '38390', '38400'),'1','0'),
         m_uterus = if_else(SIGUNGU_CD  %in%  c('37030', '37080', '37310', '37320', '38390', '38400'),'1','0'),
         m_prostate = if_else(SIGUNGU_CD  %in%  c('37012', '37020', '37040', '37070', '37080', '37310', '37320', '37360', '38080'),'1','0'),
         i_lung_m = if_else(SIGUNGU_CD  %in%  c('35040', '35050', '35340', '35350', '35360', '35370', '35380', '36040', '36310', '36320', '36330', '36350', 
                                            '36360', '36370', '36380', '36390', '36400', '36410', '36420', '36430', '36440', '36450', '36480', '38360', 
                                            '38370', '38380'),'1','0'),
         i_lung_f = if_else(SIGUNGU_CD  %in%  c('33340', '35330', '37020', '37030', '37070', '37080', '37090', '37310', '37320', '37330', '37340', '37350', 
                                            '37360', '37370', '37380', '37400', '37410', '38080', '38330', '38390', '38400'),'1','0'),
         i_stomach_m = if_else(SIGUNGU_CD  %in%  c('33320', '33330', '33340', '34020', '34060', '34310', '34330', '34340', '35040', '35050', '35060', '35310', 
                                               '35320', '35330', '35340', '35350', '35360', '35380', '36310', '36320', '36330', '37030', '37380', '38370', 
                                               '38380', '38390', '38400'),'1','0'),
         i_stomach_f = if_else(SIGUNGU_CD  %in%  c('35050', '35340' ,'35350', '35360','35370', '36040', '36310', '36320', '36330', '36350', '36360', '36370' ,
                                               '36380' ,'36390', '36400', '36410', '36420', '36430', '36440', '36450' ,'36480','38060', '38340', '38360' ,
                                               '38370' ,'38380'),'1','0'),
         i_colon_m = if_else(SIGUNGU_CD  %in%  c('25020', '33330', '34020', '34030', '34060', '34310', '34320', '34330', '34340', '34350', '34360', '34370' ,
                                             '34390', '35040', '35060', '35310', '35320', '35340', '35350', '35360', '35370', '35380'),'1','0'),
         i_colon_f = if_else(SIGUNGU_CD  %in%  c('33320', '33330', '33340', '33360', '37020', '37030', '37040', '37060', '37070', '37080', '37090', '37310', 
                                             '37320', '37330', '37340', '37350', '37360', '37370', '37380', '37400'),'1','0'),
         i_liver_m = if_else(SIGUNGU_CD  %in%  c('35040' ,'35050', '35350', '35360', '35370' ,'36040', '36310', '36320', '36330', '36350', '36360', '36370', 
                                             '36380' ,'36390' ,'36400', '36410' ,'36420', '36430' ,'36440' ,'36450', '36480', '38030', '38050', '38060', 
                                             '38340', '38360', '38370' ,'38380'),'1','0'),
         i_liver_f = if_else(SIGUNGU_CD %in% c('21010', '21020', '21100', '21120' ,'35340', '36320', '36330', '36350' ,'36360', '36370' ,'37370', '38050' ,
                                             '38060', '38070', '38080', '38310', '38320', '38330', '38340' ,'38360' ,'38370' ,'38380' ,'38390' ,'38400'),'1','0'),
         i_breast = if_else(SIGUNGU_CD %in% c('11010', '11020', '11030' ,'11040' ,'11050' ,'11060', '11070', '11080', '11110', '11120', '11130', '11140', 
                                            '11150', '11160' ,'11180' ,'11190' ,'11200', '11210', '11220', '11230', '11240', '11250' ,'31101', '31102', 
                                            '31102' ,'31110' ,'31260'),'1','0'),
         m_breast = if_else(SIGUNGU_CD %in% c('11010', '11020', '11030', '11120', '11140', '31101', '31260', '31350'),'1','0'),
         i_uterus = if_else(SIGUNGU_CD %in% c('31260', '31270'),'1','0'),
         i_prostate = if_else(SIGUNGU_CD %in% c('11010', '11020', '11030', '11040', '11080', '11090', '11130', '11140', '11220', '11230', '11240', '31101', 
                                              '31102', '31102', '31200', '31260', '31270', '31370', '31380', '32010'),'1','0'),
         i_thyroid_m = if_else(SIGUNGU_CD %in% c('24010', '24020', '24030', '24040', '36020', '36030', '36060', '36360', '36370', '36380', '36430', '36440' ),'1','0'),
         i_thyroid_f = if_else(SIGUNGU_CD %in% c('24010', '24020', '24030', '24040', '24050', '36010', '36020', '36030', '36060', '36320', '36350', '36360', 
                                                 '36380', '36410', '36420', '36430', '36440', '36450'),'1','0'))
# var test -> ttest

# male mortality

test_factor = function(x,y) {
  a = var.test(x~y)
  if (a$p.value <= 0.05){
    b = t.test(x~y, alternative = "less", var.equal = F)}
  else{
    b = t.test(x~y, alternative = "less", var.equal = T)}
  print(b$p.value)
  print(b$estimate)

  }
  

# mortality

# ses

test_factor(data_factor$elderly_2005, data_factor$m_lung_m)

test_factor(data_factor$elderly_2005, data_factor$m_stomach_m)

test_factor(data_factor$elderly_2005, data_factor$m_colon_m)

test_factor(data_factor$elderly_2005, data_factor$m_liver_m)

test_factor(data_factor$elderly_2005, data_factor$m_prostate)


test_factor(data_factor$elderly_2005, data_factor$m_lung_f)

test_factor(data_factor$elderly_2005, data_factor$m_stomach_f)

test_factor(data_factor$elderly_2005, data_factor$m_colon_f)

test_factor(data_factor$elderly_2005, data_factor$m_liver_f)

test_factor(data_factor$elderly_2005, data_factor$m_uterus)

test_factor(data_factor$elderly_2005, data_factor$m_breast)


#edu

test_factor(data_factor$edu_2005, data_factor$m_lung_m)

test_factor(data_factor$edu_2005, data_factor$m_stomach_m)

test_factor(data_factor$edu_2005, data_factor$m_colon_m)

test_factor(data_factor$edu_2005, data_factor$m_liver_m)

test_factor(data_factor$edu_2005, data_factor$m_prostate)


test_factor(data_factor$edu_2005, data_factor$m_lung_f)

test_factor(data_factor$edu_2005, data_factor$m_stomach_f)

test_factor(data_factor$edu_2005, data_factor$m_colon_f)

test_factor(data_factor$edu_2005, data_factor$m_liver_f)

test_factor(data_factor$edu_2005, data_factor$m_uterus)

test_factor(data_factor$edu_2005, data_factor$m_breast)


#population

test_factor(data_factor$pop_2005, data_factor$m_lung_m)

test_factor(data_factor$pop_2005, data_factor$m_stomach_m)

test_factor(data_factor$pop_2005, data_factor$m_colon_m)

test_factor(data_factor$pop_2005, data_factor$m_liver_m)

test_factor(data_factor$pop_2005, data_factor$m_prostate) # 


test_factor(data_factor$pop_2005, data_factor$m_lung_f)

test_factor(data_factor$pop_2005, data_factor$m_stomach_f)

test_factor(data_factor$pop_2005, data_factor$m_colon_f)

test_factor(data_factor$pop_2005, data_factor$m_liver_f)

test_factor(data_factor$pop_2005, data_factor$m_uterus)

test_factor(data_factor$pop_2005, data_factor$m_breast)


#grdp
test_factor(data_factor$GRDP_2005, data_factor$m_lung_m)

test_factor(data_factor$GRDP_2005, data_factor$m_stomach_m)

test_factor(data_factor$GRDP_2005, data_factor$m_colon_m)

test_factor(data_factor$GRDP_2005, data_factor$m_liver_m)

test_factor(data_factor$GRDP_2005, data_factor$m_prostate)


test_factor(data_factor$GRDP_2005, data_factor$m_lung_f)

test_factor(data_factor$GRDP_2005, data_factor$m_stomach_f)

test_factor(data_factor$GRDP_2005, data_factor$m_colon_f)

test_factor(data_factor$GRDP_2005, data_factor$m_liver_f)

test_factor(data_factor$GRDP_2005, data_factor$m_uterus)

test_factor(data_factor$GRDP_2005, data_factor$m_breast)



#health behaviors

test_factor(data_factor$walking, data_factor$m_lung_m) #

test_factor(data_factor$walking, data_factor$m_stomach_m) #

test_factor(data_factor$walking, data_factor$m_colon_m) #

test_factor(data_factor$walking, data_factor$m_liver_m) #

test_factor(data_factor$walking, data_factor$m_prostate)


test_factor(data_factor$walking, data_factor$m_lung_f)

test_factor(data_factor$walking, data_factor$m_stomach_f)

test_factor(data_factor$walking, data_factor$m_colon_f) #

test_factor(data_factor$walking, data_factor$m_liver_f) #

test_factor(data_factor$walking, data_factor$m_uterus)

test_factor(data_factor$walking, data_factor$m_breast)


# drinking

test_factor(data_factor$dinking, data_factor$m_lung_m)

test_factor(data_factor$dinking, data_factor$m_stomach_m)

test_factor(data_factor$dinking, data_factor$m_colon_m)

test_factor(data_factor$dinking, data_factor$m_liver_m)

test_factor(data_factor$dinking, data_factor$m_prostate) 


test_factor(data_factor$dinking, data_factor$m_lung_f)

test_factor(data_factor$dinking, data_factor$m_stomach_f)

test_factor(data_factor$dinking, data_factor$m_colon_f)

test_factor(data_factor$dinking, data_factor$m_liver_f)

test_factor(data_factor$dinking, data_factor$m_uterus)

test_factor(data_factor$dinking, data_factor$m_breast)


# antismoking

test_factor(data_factor$antismoking, data_factor$m_lung_m)

test_factor(data_factor$antismoking, data_factor$m_stomach_m)

test_factor(data_factor$antismoking, data_factor$m_colon_m)

test_factor(data_factor$antismoking, data_factor$m_liver_m)  #

test_factor(data_factor$antismoking, data_factor$m_prostate) #


test_factor(data_factor$antismoking, data_factor$m_lung_f)

test_factor(data_factor$antismoking, data_factor$m_stomach_f)

test_factor(data_factor$antismoking, data_factor$m_colon_f) #

test_factor(data_factor$antismoking, data_factor$m_liver_f) #

test_factor(data_factor$antismoking, data_factor$m_uterus) #

test_factor(data_factor$antismoking, data_factor$m_breast) #


# smoking
test_factor(data_factor$smoking, data_factor$m_lung_m)

test_factor(data_factor$smoking, data_factor$m_stomach_m) #

test_factor(data_factor$smoking, data_factor$m_colon_m)

test_factor(data_factor$smoking, data_factor$m_liver_m)

test_factor(data_factor$smoking, data_factor$m_prostate)  #
 

test_factor(data_factor$smoking, data_factor$m_lung_f) #

test_factor(data_factor$smoking, data_factor$m_stomach_f) #

test_factor(data_factor$smoking, data_factor$m_colon_f) #

test_factor(data_factor$smoking, data_factor$m_liver_f) 

test_factor(data_factor$smoking, data_factor$m_uterus) #

test_factor(data_factor$smoking, data_factor$m_breast) #


# obesity
test_factor(data_factor$obesity, data_factor$m_lung_m)

test_factor(data_factor$obesity, data_factor$m_stomach_m)

test_factor(data_factor$obesity, data_factor$m_colon_m)# 

test_factor(data_factor$obesity, data_factor$m_liver_m)

test_factor(data_factor$obesity, data_factor$m_prostate) #


test_factor(data_factor$obesity, data_factor$m_lung_f) #

test_factor(data_factor$obesity, data_factor$m_stomach_f)

test_factor(data_factor$obesity, data_factor$m_colon_f) #

test_factor(data_factor$obesity, data_factor$m_liver_f)

test_factor(data_factor$obesity, data_factor$m_uterus) #

test_factor(data_factor$obesity, data_factor$m_breast) #


# stress
test_factor(data_factor$stress, data_factor$m_lung_m)

test_factor(data_factor$stress, data_factor$m_stomach_m)

test_factor(data_factor$stress, data_factor$m_colon_m) #

test_factor(data_factor$stress, data_factor$m_liver_m)

test_factor(data_factor$stress, data_factor$m_prostate)


test_factor(data_factor$stress, data_factor$m_lung_f)

test_factor(data_factor$stress, data_factor$m_stomach_f)

test_factor(data_factor$stress, data_factor$m_colon_f) #

test_factor(data_factor$stress, data_factor$m_liver_f) #

test_factor(data_factor$stress, data_factor$m_uterus) #

test_factor(data_factor$stress, data_factor$m_breast) #


# low salt use

test_factor(data_factor$lowsaltuse, data_factor$m_lung_m) #
 
test_factor(data_factor$lowsaltuse, data_factor$m_stomach_m) #

test_factor(data_factor$lowsaltuse, data_factor$m_colon_m)  #

test_factor(data_factor$lowsaltuse, data_factor$m_liver_m) #

test_factor(data_factor$lowsaltuse, data_factor$m_prostate) #


test_factor(data_factor$lowsaltuse, data_factor$m_lung_f) #

test_factor(data_factor$lowsaltuse, data_factor$m_stomach_f)#

test_factor(data_factor$lowsaltuse, data_factor$m_colon_f)#

test_factor(data_factor$lowsaltuse, data_factor$m_liver_f)#

test_factor(data_factor$lowsaltuse, data_factor$m_uterus)#

test_factor(data_factor$lowsaltuse, data_factor$m_breast)


# screening

test_factor(data_factor$P_stomach_m, data_factor$i_stomach_m)

test_factor(data_factor$P_stomach_f, data_factor$i_stomach_f)

test_factor(data_factor$P_colon_m, data_factor$i_colon_m) #

test_factor(data_factor$P_colon_f, data_factor$i_colon_f) #

test_factor(data_factor$P_liver_m, data_factor$i_liver_m)

test_factor(data_factor$P_liver_f, data_factor$i_liver_f)

test_factor(data_factor$P_breast, data_factor$i_breast) 

test_factor(data_factor$P_cervical, data_factor$i_uterus) #


# test for incidence

# ses
test_factor(data_factor$elderly_2005, data_factor$i_lung_m)

test_factor(data_factor$elderly_2005, data_factor$i_stomach_m)

test_factor(data_factor$elderly_2005, data_factor$i_colon_m)

test_factor(data_factor$elderly_2005, data_factor$i_liver_m)

test_factor(data_factor$elderly_2005, data_factor$i_prostate)

test_factor(data_factor$elderly_2005, data_factor$i_thyroid_m)


test_factor(data_factor$elderly_2005, data_factor$i_lung_f)

test_factor(data_factor$elderly_2005, data_factor$i_stomach_f)

test_factor(data_factor$elderly_2005, data_factor$i_colon_f)

test_factor(data_factor$elderly_2005, data_factor$i_liver_f)

test_factor(data_factor$elderly_2005, data_factor$i_uterus)

test_factor(data_factor$elderly_2005, data_factor$i_breast)

test_factor(data_factor$elderly_2005, data_factor$i_thyroid_f)



#edu

test_factor(data_factor$edu_2005, data_factor$i_lung_m)

test_factor(data_factor$edu_2005, data_factor$i_stomach_m)

test_factor(data_factor$edu_2005, data_factor$i_colon_m)

test_factor(data_factor$edu_2005, data_factor$i_liver_m)

test_factor(data_factor$edu_2005, data_factor$i_prostate)

test_factor(data_factor$edu_2005, data_factor$i_thyroid_m)


test_factor(data_factor$edu_2005, data_factor$i_lung_f)

test_factor(data_factor$edu_2005, data_factor$i_stomach_f)

test_factor(data_factor$edu_2005, data_factor$i_colon_f)

test_factor(data_factor$edu_2005, data_factor$i_liver_f)

test_factor(data_factor$edu_2005, data_factor$i_uterus) #

test_factor(data_factor$edu_2005, data_factor$i_breast) #

test_factor(data_factor$edu_2005, data_factor$i_thyroid_f) #



#GRPD

test_factor(data_factor$GRDP_2005, data_factor$i_lung_m)

test_factor(data_factor$GRDP_2005, data_factor$i_stomach_m)

test_factor(data_factor$GRDP_2005, data_factor$i_colon_m)

test_factor(data_factor$GRDP_2005, data_factor$i_liver_m)

test_factor(data_factor$GRDP_2005, data_factor$i_prostate)

test_factor(data_factor$GRDP_2005, data_factor$i_thyroid_m)



test_factor(data_factor$GRDP_2005, data_factor$i_lung_f)

test_factor(data_factor$GRDP_2005, data_factor$i_stomach_f)

test_factor(data_factor$GRDP_2005, data_factor$i_colon_f)

test_factor(data_factor$GRDP_2005, data_factor$i_liver_f)

test_factor(data_factor$GRDP_2005, data_factor$i_uterus) #

test_factor(data_factor$GRDP_2005, data_factor$i_breast) #

test_factor(data_factor$GRDP_2005, data_factor$i_thyroid_f) #



#health behaviors

test_factor(data_factor$walking, data_factor$i_lung_m) #

test_factor(data_factor$walking, data_factor$i_stomach_m)#

test_factor(data_factor$walking, data_factor$i_colon_m)#

test_factor(data_factor$walking, data_factor$i_liver_m)#

test_factor(data_factor$walking, data_factor$i_prostate)#

test_factor(data_factor$walking, data_factor$i_thyroid_m)#


test_factor(data_factor$walking, data_factor$i_lung_f)

test_factor(data_factor$walking, data_factor$i_stomach_f)#

test_factor(data_factor$walking, data_factor$i_colon_f)#

test_factor(data_factor$walking, data_factor$i_liver_f)#

test_factor(data_factor$walking, data_factor$i_uterus)#

test_factor(data_factor$walking, data_factor$i_breast)#

test_factor(data_factor$walking, data_factor$i_thyroid_f)#




# drinking

test_factor(data_factor$dinking, data_factor$i_lung_m)

test_factor(data_factor$dinking, data_factor$i_stomach_m)

test_factor(data_factor$dinking, data_factor$i_colon_m)

test_factor(data_factor$dinking, data_factor$i_liver_m)

test_factor(data_factor$dinking, data_factor$i_prostate)

test_factor(data_factor$dinking, data_factor$i_thyroid_m)


test_factor(data_factor$dinking, data_factor$i_lung_f)

test_factor(data_factor$dinking, data_factor$i_stomach_f)

test_factor(data_factor$dinking, data_factor$i_colon_f)

test_factor(data_factor$dinking, data_factor$i_liver_f)

test_factor(data_factor$dinking, data_factor$i_uterus)

test_factor(data_factor$dinking, data_factor$i_breast)

test_factor(data_factor$dinking, data_factor$i_thyroid_f)


# antismoking

test_factor(data_factor$antismoking, data_factor$i_lung_m)

test_factor(data_factor$antismoking, data_factor$i_stomach_m)

test_factor(data_factor$antismoking, data_factor$i_colon_m)

test_factor(data_factor$antismoking, data_factor$i_liver_m)

test_factor(data_factor$antismoking, data_factor$i_prostate)

test_factor(data_factor$antismoking, data_factor$i_thyroid_m)



test_factor(data_factor$antismoking, data_factor$i_lung_f)#

test_factor(data_factor$antismoking, data_factor$i_stomach_f)

test_factor(data_factor$antismoking, data_factor$i_colon_f)#

test_factor(data_factor$antismoking, data_factor$i_liver_f)#

test_factor(data_factor$antismoking, data_factor$i_uterus)#

test_factor(data_factor$antismoking, data_factor$i_breast)#

test_factor(data_factor$antismoking, data_factor$i_thyroid_f)


# smoking
test_factor(data_factor$smoking, data_factor$i_lung_m)

test_factor(data_factor$smoking, data_factor$i_stomach_m)

test_factor(data_factor$smoking, data_factor$i_colon_m)

test_factor(data_factor$smoking, data_factor$i_liver_m)

test_factor(data_factor$smoking, data_factor$i_prostate)#

test_factor(data_factor$smoking, data_factor$i_thyroid_m)


test_factor(data_factor$smoking, data_factor$i_lung_f)#

test_factor(data_factor$smoking, data_factor$i_stomach_f)#

test_factor(data_factor$smoking, data_factor$i_colon_f)

test_factor(data_factor$smoking, data_factor$i_liver_f)#

test_factor(data_factor$smoking, data_factor$i_uterus)#

test_factor(data_factor$smoking, data_factor$i_breast)

test_factor(data_factor$smoking, data_factor$i_thyroid_f)


# obesity
test_factor(data_factor$obesity, data_factor$i_lung_m)#

test_factor(data_factor$obesity, data_factor$i_stomach_m)#

test_factor(data_factor$obesity, data_factor$i_colon_m)#

test_factor(data_factor$obesity, data_factor$i_liver_m)

test_factor(data_factor$obesity, data_factor$i_prostate)#

test_factor(data_factor$obesity, data_factor$i_thyroid_m)#



test_factor(data_factor$obesity, data_factor$i_lung_f)#

test_factor(data_factor$obesity, data_factor$i_stomach_f)

test_factor(data_factor$obesity, data_factor$i_colon_f)#

test_factor(data_factor$obesity, data_factor$i_liver_f)

test_factor(data_factor$obesity, data_factor$i_uterus)#

test_factor(data_factor$obesity, data_factor$i_breast)#

test_factor(data_factor$obesity, data_factor$i_thyroid_f)


# stress
test_factor(data_factor$stress, data_factor$i_lung_m)

test_factor(data_factor$stress, data_factor$i_stomach_m)#

test_factor(data_factor$stress, data_factor$i_colon_m)#

test_factor(data_factor$stress, data_factor$i_liver_m)

test_factor(data_factor$stress, data_factor$i_prostate)#

test_factor(data_factor$stress, data_factor$i_thyroid_m)


test_factor(data_factor$stress, data_factor$i_lung_f)

test_factor(data_factor$stress, data_factor$i_stomach_f)

test_factor(data_factor$stress, data_factor$i_colon_f)

test_factor(data_factor$stress, data_factor$i_liver_f)#

test_factor(data_factor$stress, data_factor$i_uterus)#

test_factor(data_factor$stress, data_factor$i_breast)#

test_factor(data_factor$stress, data_factor$i_thyroid_f)


# low salt use

test_factor(data_factor$lowsaltuse, data_factor$i_lung_m)#

test_factor(data_factor$lowsaltuse, data_factor$i_stomach_m)#

test_factor(data_factor$lowsaltuse, data_factor$i_colon_m)#

test_factor(data_factor$lowsaltuse, data_factor$i_liver_m)

test_factor(data_factor$lowsaltuse, data_factor$i_prostate)

test_factor(data_factor$lowsaltuse, data_factor$i_thyroid_m)


test_factor(data_factor$lowsaltuse, data_factor$i_lung_f)

test_factor(data_factor$lowsaltuse, data_factor$i_stomach_f)

test_factor(data_factor$lowsaltuse, data_factor$i_colon_f)

test_factor(data_factor$lowsaltuse, data_factor$i_liver_f)#

test_factor(data_factor$lowsaltuse, data_factor$i_uterus)

test_factor(data_factor$lowsaltuse, data_factor$i_breast)

test_factor(data_factor$lowsaltuse, data_factor$i_thyroid_f)





# cluster of incidence after standardize
flexrun_inc1 = function(g,c,y,k){
  
  data <- shap2010_sub%>% 
    mutate(Sgg_cd = plyr::mapvalues(SIGUNGU_CD, code_change1$fromcode, code_change1$tocode)) %>% 
    left_join(., cluster_inc, by = 'Sgg_cd') %>% 
    filter(gender == g & cause == c & year_agg == y) 
  
  names_vec    <- data$Sgg_cd
  
  adj_matrix   <- nb2mat(poly2nb(data), style = 'B', zero.policy = TRUE)
  
  coord_matrix <- as.matrix(st_coordinates(st_centroid(data)))
  
  case_max  <- data %>% st_drop_geometry
  
  case_matrix <- as.matrix(case_max %>%.[,'counts'])
  
  pmat <- data$Population
  
  expected <- data$ex_inc
  
  rflexscan(x       = coord_matrix[,1],
            y       = coord_matrix[,2],
            name    = names_vec,
            observed = as.vector(case_matrix),
            expected = expected,
            nb         = adj_matrix,
            stattype   = 'RESTRICTED',
            scanmethod = 'FLEXIBLE',
            ralpha     = 0.05, 
            simcount   = 999,
            rantype    = 'POISSON',
            comments   = "", verbose = FALSE, secondary = NULL,
            clustersize = k)
  
}

# ________________________________________________-
#Main analysis for size = 15%

# run function

Male_lung_inc1_37_ex  <- flexrun_inc1('Male', 'lung', '1999-2003',37)


# Stomach

Male_stomach_inc1_37_ex    <- flexrun_inc1('Male', 'stomach', '1999-2003',37)

# Colon

Male_colon_inc1_37_ex     <- flexrun_inc1( 'Male', 'colon', '1999-2003',37)


#liver cancer

Male_liver_inc1_37_ex    <- flexrun_inc1('Male', 'liver', '1999-2003',37)


# Prostate
Male_prostate_inc1_37_ex    <- flexrun_inc1('Male', 'prostate', '1999-2003',37)
# thyroid

Male_thyroid_inc1_37_ex    <- flexrun_inc1('Male', 'thyroid', '1999-2003',37)


#----------------------------------------------------------------------

#2004-2008

#Lung cancer

Male_lung_inc2_37_ex    <- flexrun_inc1('Male', 'lung', '2004-2008',37)



# Stomach

Male_stomach_inc2_37_ex    <- flexrun_inc1('Male', 'stomach', '2004-2008',37)


# Colon

Male_colon_inc2_37_ex     <- flexrun_inc1( 'Male', 'colon', '2004-2008',37)


#liver cancer

Male_liver_inc2_37_ex     <- flexrun_inc1('Male', 'liver', '2004-2008',37)


# Prostate
Male_prostate_inc2_37_ex    <- flexrun_inc1('Male', 'prostate', '2004-2008',37)

# thyroid

Male_thyroid_inc2_37_ex     <- flexrun_inc1('Male', 'thyroid', '2004-2008',37)

#_______________________________________________________________
#2009-2013

#Lung cancer

Male_lung_inc3_37_ex    <- flexrun_inc1('Male', 'lung', '2009-2013',37)

# Stomach

Male_stomach_inc3_37_ex    <- flexrun_inc1('Male', 'stomach', '2009-2013',37)


# Colon

Male_colon_inc3_37_ex     <- flexrun_inc1( 'Male', 'colon', '2009-2013',37)


#liver cancer


Male_liver_inc3_37_ex     <- flexrun_inc1('Male', 'liver', '2009-2013',37)



# Prostate
Male_prostate_inc3_37_ex    <- flexrun_inc1('Male', 'prostate', '2009-2013',37)

# thyroid 

Male_thyroid_inc3_37_ex     <- flexrun_inc1('Male', 'thyroid', '2009-2013',37)

# female
# run function

Female_lung_inc1_37_ex  <- flexrun_inc1('Female', 'lung', '1999-2003',37)


# Stomach

Female_stomach_inc1_37_ex    <- flexrun_inc1('Female', 'stomach', '1999-2003',37)

# Colon

Female_colon_inc1_37_ex     <- flexrun_inc1( 'Female', 'colon', '1999-2003',37)


#liver cancer

Female_liver_inc1_37_ex    <- flexrun_inc1('Female', 'liver', '1999-2003',37)


# Breast
Female_breast_inc1_37_ex    <- flexrun_inc1('Female', 'breast', '1999-2003',37)
# uterus
Female_uterus_inc1_37_ex    <- flexrun_inc1('Female', 'uterus', '1999-2003',37)

# thyroid

Female_thyroid_inc1_37_ex    <- flexrun_inc1('Female', 'thyroid', '1999-2003',37)


#----------------------------------------------------------------------

#2004-2008

#Lung cancer

Female_lung_inc2_37_ex    <- flexrun_inc1('Female', 'lung', '2004-2008',37)



# Stomach

Female_stomach_inc2_37_ex    <- flexrun_inc1('Female', 'stomach', '2004-2008',37)


# Colon

Female_colon_inc2_37_ex     <- flexrun_inc1( 'Female', 'colon', '2004-2008',37)


#liver cancer

Female_liver_inc2_37_ex     <- flexrun_inc1('Female', 'liver', '2004-2008',37)


# Breast
Female_breast_inc2_37_ex    <- flexrun_inc1('Female', 'breast', '2004-2008',37)
# uterus
Female_uterus_inc2_37_ex    <- flexrun_inc1('Female', 'uterus', '2004-2008',37)

# thyroid

Female_thyroid_inc2_37_ex     <- flexrun_inc1('Female', 'thyroid', '2004-2008',37)

#_______________________________________________________________
#2009-2013

#Lung cancer

Female_lung_inc3_37_ex    <- flexrun_inc1('Female', 'lung', '2009-2013',37)

# Stomach

Female_stomach_inc3_37_ex    <- flexrun_inc1('Female', 'stomach', '2009-2013',37)


# Colon

Female_colon_inc3_37_ex     <- flexrun_inc1( 'Female', 'colon', '2009-2013',37)


#liver cancer


Female_liver_inc3_37_ex     <- flexrun_inc1('Female', 'liver', '2009-2013',37)


# Breast
Female_breast_inc3_37_ex    <- flexrun_inc1('Female', 'breast', '2009-2013',37)
# uterus
Female_uterus_inc3_37_ex    <- flexrun_inc1('Female', 'uterus', '2009-2013',37)

# thyroid 

Female_thyroid_inc3_37_ex     <- flexrun_inc1('Female', 'thyroid', '2009-2013',37)

# plot
png(filename = "Inc_male_cluster_all_ex.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(6,3))

choropleth(shap2010_sub$geometry, Male_lung_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_lung_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_lung_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_stomach_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_stomach_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Male_colon_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_colon_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_colon_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_liver_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_liver_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_prostate_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_prostate_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Prostate cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_thyroid_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_thyroid_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Male_thyroid_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Male_thyroid_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 3', cex.sub = 7)
dev.off()

# 
png(filename = "Inc_female_cluster_all_ex.png", 
    units = "in", 
    width = 30, 
    height = 40, 
    res = 300)

par(mfrow=c(7,3))

choropleth(shap2010_sub$geometry, Female_lung_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_lung_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_lung_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Lung cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_stomach_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_stomach_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Stomach cancer - Period 3', cex.sub = 7)


choropleth(shap2010_sub$geometry, Female_colon_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_colon_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_colon_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Colorectal cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_liver_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_liver_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Liver cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_breast_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_breast_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Breast cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_uterus_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_uterus_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Uterus cancer - Period 3', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc1_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_thyroid_inc1_37_ex$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 1', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc2_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_thyroid_inc2_37_ex$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 2', cex.sub = 7)

choropleth(shap2010_sub$geometry, Female_thyroid_inc3_37_ex,
           col = rev(pal1(30)),
           region_color = 'white',
           rank = 1:length(Female_thyroid_inc3_37_ex$cluster),
           pval = 0.05)
title(sub = 'Thyroid cancer - Period 3', cex.sub = 7)
dev.off()

