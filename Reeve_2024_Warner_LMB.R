#### Packages ####
.libPaths("C:/Users/NEW USER/Documents/R/win-library/4.0")
library(data.table)
library(readr)
library(dplyr)
library(lubridate)
library(lme4)
library(nlme)
library(lmerTest)
library(car)
library(emmeans)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggpubr) # for ggqqplots
library(rstatix) # for shapiro normality test and kruskal test
library(multcomp)
library(multcompView)
library(MASS)
library(scales)
library(ggThemeAssist)
library(MuMIn)

####Accel data ####
dat<-fread("D:\\2021 Recovered Tags - Manual Download\\Connor Reeve accell project\\Data\\Max Minute Data\\alltagsmin.csv")

#removing duplicates and erroneous datapoints

dat<-dat %>% 
  distinct() %>% 
  mutate(depth_m = if_else(depth_m == 0, NA_real_, depth_m)) %>% 
  mutate(temp_m = if_else(temp_m < 0, NA_real_, temp_m)) %>% 
  mutate(temp_m = if_else(temp_m > 30, NA_real_, temp_m))


#importing morphometric data (includes release and recapture data also)
dat_morphs<-fread("D:\\2021 Recovered Tags - Manual Download\\Connor Reeve accell project\\Data\\Max Minute Data\\morphometrics.csv")

#converting ODBA from gs to m/s2 - necessary for Reeve et al. 2024 calibrations

dat$ODBA_ms2<-dat$ODBA_m*9.80665
dat$ODBA_ms2_sd<-dat$ODBA_sd*9.80665
dat$ODBA_ms2_max<-dat$ODBA_max*9.80665
dat$ODBA_ms2_min<-dat$ODBA_min*9.80665
dat$ODBA_ms2_sum<-dat$ODBA_sum*9.80665

#calculating speed as per Reeve et al. 2024
dat$speed_bls<-10^(1.3156+0.0277*dat$temp_m+0.3423*log10(dat$ODBA_ms2)-0.6505*log10(dat$length))

dat<-dat %>% 
  mutate(speed_bls = if_else(speed_bls < 0, 0, speed_bls)) %>% 
  mutate(speed_bls = if_else(ODBA_ms2 < 0.1, 0, speed_bls))

# averaging / totaling by day

glimpse(dat)

df_hr_m<-dat %>% 
  group_by(id,hour) %>%
  summarise(temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            speed_bls_m = mean(speed_bls),
            ODBA_ms2 = mean(ODBA_ms2),
            ODBA_ms2_sum = sum(ODBA_ms2_sum),
            ODBA_ms2_max = max(ODBA_ms2_max),
            ODBA_ms2_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(hour,id,fill=list(temp_m=NA,depth_m=NA,ODBA_ms2=NA, ODBA_ms2_sum=NA, ODBA_ms2_max=NA, ODBA_ms2_min=NA))


df_hr_m<-dat %>% 
  group_by(id,hour) %>%
  summarise(temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            speed_bls = mean(speed_bls),
            ODBA_ms2 = mean(ODBA_ms2),
            ODBA_ms2_sum = sum(ODBA_ms2_sum),
            ODBA_ms2_max = max(ODBA_ms2_max),
            ODBA_ms2_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(hour,id,fill=list(temp_m=NA,depth_m=NA,speed_bls = NA, ODBA_ms2=NA, ODBA_ms2_sum=NA, ODBA_ms2_max=NA, ODBA_ms2_min=NA))

df_day_m<-dat %>% 
  group_by(id,day) %>%
  summarise(temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            speed_bls = mean(speed_bls),
            ODBA_ms2 = mean(ODBA_ms2),
            ODBA_ms2_sum = sum(ODBA_ms2_sum),
            ODBA_ms2_max = max(ODBA_ms2_max),
            ODBA_ms2_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(day,id,fill=list(temp_m=NA,depth_m=NA,speed_bls = NA, ODBA_ms2=NA, ODBA_ms2_sum=NA, ODBA_ms2_max=NA, ODBA_ms2_min=NA))

df_day_sd<-dat %>% 
  group_by(id,day) %>%
  summarise(timestamp = sd(timestamp),
            temp_m = sd(temp_m, na.rm = TRUE),
            speed_bls = sd(speed_bls),
            depth_m = sd(depth_m, na.rm = TRUE),
            ODBA_m = sd(ODBA_ms2)) %>% 
  ungroup() %>%
  complete(day,id,fill=list(temp_m=NA,depth_m=NA,speed_bls = NA, ODBA_ms2=NA))

# averaging / totaling by month_character_cat (i.e. month as a categorical unit)

df_month_character_m<-dat %>% 
  group_by(id,month_character) %>%
  summarise(timestamp = mean(timestamp),
            temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            speed_bls = mean(speed_bls),
            ODBA_ms2 = mean(ODBA_ms2),
            ODBA_ms2_sum = sum(ODBA_ms2_sum),
            ODBA_ms2_max = max(ODBA_ms2_max),
            ODBA_ms2_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(month_character,id,fill=list(temp_m=NA,depth_m=NA,speed_bls = NA, ODBA_ms2=NA, ODBA_ms2_sum=NA, ODBA_ms2_max=NA, ODBA_ms2_min=NA))

df_month_character_sd<-dat %>% 
  group_by(id,month_character) %>%
  summarise(temp_sd = sd(temp_m, na.rm = TRUE),
            depth_sd = sd(depth_m, na.rm = TRUE),
            speed_bls = sd(speed_bls),
            ODBA_sd = sd(ODBA_ms2)) %>%
  ungroup() %>%
  complete(month_character,id,fill=list(temp_m=NA,depth_m=NA,speed_bls = NA, ODBA_ms2=NA))

glimpse(dat)

dat_averages<-dat %>% 
  group_by(timestamp) %>% 
  summarise(depth_mean = mean(depth_m,na.rm = TRUE),
            temp_mean = mean(temp_m, na.rm = TRUE))

glimpse(dat)

#inferring temperature for some the few loggers where the temp. sensor stopped working

dat_filled<-dat %>% 
  dplyr::select(id,timestamp,temp_m,hour,depth_m,ODBA_ms2, speed_bls) %>% 
  left_join(dat_averages, by = "timestamp") %>% 
  mutate(depth_m = if_else(is.na(depth_m), depth_mean, depth_m)) %>% 
  mutate(temp_m = if_else(is.na(temp_m), temp_mean, temp_m)) %>% 
  group_by(id,hour) %>% 
  summarise(mean_depth = mean(depth_m),
            mean_temp = mean(temp_m),
            mean_ODBA = mean(ODBA_ms2),
            mean_speed_bls = mean(speed_bls)) %>% 
  dplyr::select(id,timestamp = hour,temp_m = mean_temp, depth_m = mean_depth, ODBA_ms2 = mean_ODBA, speed_bls= mean_speed_bls)

weight_df<-dat %>% 
  dplyr::select(id,weight,length) %>% 
  distinct()

dat_filled<-dat_filled %>% 
  left_join(weight_df, by = "id")

dat_filled$speed_bls<-10^(1.3156+0.0277*dat_filled$temp_m+0.3423*log10(dat_filled$ODBA_ms2)-0.6505*log10(dat_filled$length))

dat_filled<-dat_filled %>% 
  mutate(speed_bls = if_else(speed_bls < 0, 0, speed_bls)) %>% 
  mutate(speed_bls = if_else(ODBA_ms2 < 0.1, 0, speed_bls))

####HR data####

datHR<-fread("D:\\2021 Recovered Tags - Manual Download\\WinterHR_2020_2022.csv")

glimpse(datHR)

datHR$timestamp = as.POSIXct(datHR$timestamp, format="%Y-%m-%d %H:%M:%S")

glimpse(datHR)

datHR<-datHR %>% 
  dplyr::select(-RMSSD,-SDNN)

#removing measurements prior to deployment

datHR<-datHR %>% 
  filter(!(tag == "HR0853" & timestamp < "2021-02-17 17:00:00")) %>% 
  filter(!(tag == "HR0854" & timestamp < "2021-01-13 14:50:00")) %>% 
  filter(!(tag == "HR0856" & timestamp < "2021-02-02 16:11:00")) %>%
  filter(!(tag == "HR0858" & timestamp < "2021-02-21 17:06:00")) %>%
  filter(!(tag == "HR0860" & timestamp < "2021-02-21 18:03:00")) %>%
  filter(!(tag == "HR0861" & timestamp < "2021-01-13 16:08:00")) %>%
  filter(!(tag == "HR0862" & timestamp < "2021-02-06 17:32:00")) %>%
  filter(!(tag == "HR0864" & timestamp < "2020-11-13 18:43:00")) %>%
  filter(!(tag == "HR0865" & timestamp < "2020-11-13 16:08:00")) %>%
  filter(!(tag == "HR0866" & timestamp < "2020-11-06 17:05:00")) %>%
  filter(!(tag == "HR0868" & timestamp < "2021-01-23 18:40:00")) %>%
  filter(!(tag == "HR0869" & timestamp < "2020-11-06 16:04:00")) %>%
  filter(!(tag == "HR0871" & timestamp < "2020-11-06 20:00:00")) %>%
  filter(!(tag == "HR0872" & timestamp < "2020-11-07 20:14:00"))  


datHR<-datHR %>% 
  filter(!(tag == "HR0853" & timestamp > "2021-06-15 10:00:00")) %>% 
  filter(!(tag == "HR0854" & timestamp > "2021-06-06 12:20:00")) %>% 
  filter(!(tag == "HR0856" & timestamp > "2021-06-06 13:40:00")) %>%
  filter(!(tag == "HR0858" & timestamp > "2021-06-02 09:00:00")) %>%
  filter(!(tag == "HR0860" & timestamp > "2021-06-11 10:30:00")) %>%
  filter(!(tag == "HR0861" & timestamp > "2021-06-03 11:30:00")) %>%
  filter(!(tag == "HR0862" & timestamp > "2021-06-21 09:15:00")) %>%
  filter(!(tag == "HR0864" & timestamp > "2021-06-06 11:49:00")) %>%
  filter(!(tag == "HR0865" & timestamp > "2021-06-12 09:36:00")) %>%
  filter(!(tag == "HR0866" & timestamp > "2021-06-22 07:48:00")) %>%
  filter(!(tag == "HR0868" & timestamp > "2022-05-16 09:30:00")) %>%
  filter(!(tag == "HR0869" & timestamp > "2021-06-02 12:22:00")) %>%
  filter(!(tag == "HR0871" & timestamp > "2021-06-15 08:37:00")) %>%
  filter(!(tag == "HR0872" & timestamp > "2021-06-15 06:51:00")) 

df_calcHR<-datHR %>% 
  drop_na(calchr)

datHR$hour <- strftime(datHR$timestamp, format="%Y-%m-%d %H")
datHR$hour<-as.POSIXct(datHR$hour, format="%Y-%m-%d %H")

datHR$weight<-as.numeric(datHR$weight)
datHR$length<-as.numeric(datHR$length)

datHR_hr<-datHR %>% 
  group_by(tag,hour) %>% 
  summarise(temp_m = mean(temp, na.rm = TRUE),
            calchr_m = mean(calchr, na.rm = TRUE)) %>% 
  drop_na(calchr_m)

datHR_hr$weight<-as.numeric(datHR_hr$weight)
datHR_hr$length<-as.numeric(datHR_hr$length)

hr_trends_2021<-datHR_hr %>% 
  filter(hour < "2021-07-01 00:00:00")

hr_trends_2022<-datHR_hr %>% 
  filter(hour > "2022-09-01 00:00:00")

####DO data####

#use this for DO data in analyses (see below)

datDO<-fread("D:\\2021 Recovered Tags - Manual Download\\datDO_cleaned_2023.csv")

#not using measurements from the back bay - fish appear to mostly reside in the larger bay during the winter and thus these data would not be appropriate

DO_data<-datDO %>%
  dplyr::select(Tag,depth = `ApproxDepth (ft)`, Location, DO = `Dissolved Oxygen`, temp = Temperature, hr) %>% 
  group_by(Tag, hr) %>% 
  summarise(DO_m = mean(DO,na.rm = TRUE),
            depth_m = mean(depth, na.rm = TRUE),
            temp_m = mean(temp,na.rm = TRUE)) %>% 
  na.omit()

DO_data$hr<-as.POSIXct(DO_data$hr, format = "%Y-%m-%d %H")
DO_data$nrow <- 1:length(DO_data$hr)

#need to switch from feet to meters
DO_data$depth_m<-DO_data$depth_m/3.28084

DO_data$date<-as.numeric(strftime(DO_data$hr, format = "%j"))

#producing a randomforest model to predict depth onto fish
DOtrain <- DO_data %>% sample_frac(0.7) # 70/30 split. 
DOtest <- DO_data[!(DO_data$nrow %in% DOtrain$nrow),]

glimpse(DO_data)

z3<-formula(DO_m~temp_m + depth_m+date) #this model performs best 

library(randomForest)
library(caret)

DOForest <- randomForest(formula=z3, data = DOtrain, replace=FALSE, na.action=na.omit,
                         importance=TRUE, do.trace=100, ntree=1000)

print(DOForest)

#assessing error
obs_pred_DO = data.frame(actualValue = DOtest$DO_m, 
                         predictedValue = predict(DOForest,DOtest),
                         timestamp = DOtest$hr)

obs_pred_DO$diff<-abs(obs_pred_DO$actualValue - obs_pred_DO$predictedValue)

avgerror<-obs_pred_DO %>% 
  summarise(error_m = mean(diff),
            error_sd = sd(diff))

DO_obs<-obs_pred_DO %>% 
  dplyr::select(DO = actualValue,timestamp)
DO_obs$type<-"Observed"

DO_pred<-obs_pred_DO %>% 
  dplyr::select(DO = predictedValue,timestamp)
DO_pred$type<-"Predicted"

DO_pred_obs<-rbind(DO_obs,DO_pred)

DO_pred_obs$timestamp<-as.POSIXct(DO_pred_obs$timestamp)
DO_pred_df<-merge(data.frame(date = seq(0,365,1)), data.frame(depth_m = seq(0,7, 0.1)))

#predicting onto fish using hourly data (the resolution of the DO loggers)
fish_temp<-df_hr_m %>% 
  dplyr::select(hour, id, temp_m) %>% 
  na.omit()

fish_temp$date<-as.numeric(strftime(fish_temp$hour, format = "%j"))
fish_temp<-fish_temp %>% 
  dplyr::select(date, temp_m) %>% 
  na.omit()

DO_pred_df<-DO_pred_df %>% 
  left_join(fish_temp, by = "date")

DO_pred_df$DO<-predict(DOForest,DO_pred_df)
DO_pred_df<-DO_pred_df %>% 
  mutate(DO = if_else(date > 300, DO, 
                      if_else(date < 150, DO, NA_real_))) %>% 
  drop_na(DO)

DO_pred_df<-DO_pred_df %>% 
  mutate(date_fall = if_else(date>299, date - 299, date + 65))

#plotting data

DOForest_plot_2021<-ggplot()+
  geom_point(data = DO_pred_obs, aes(x = timestamp, y = DO, colour = type), alpha = 0.6)+
  scale_y_continuous(breaks = seq(0,150,25), limits = c(0,150))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2020/2021")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.position = "none") +
  labs(y = "Dissolved Oxygen (%)", x = "Depth (m)")+
  scale_colour_manual(values = c("lightgrey","#3399FF"))

DOForest_plot_2022<-ggplot()+
  geom_point(data = DO_pred_obs, aes(x = timestamp, y = DO, colour = type), alpha = 0.6)+
  scale_y_continuous(breaks = seq(0,150,25), limits = c(0,150))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2021/2022")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.position = "none") +
  labs(y = "Dissolved Oxygen (%)", x = "Depth (m)")+
  scale_colour_manual(values = c("lightgrey","#3399FF"))

DOForest_plot_2023<-ggplot()+
  geom_point(data = DO_pred_obs, aes(x = timestamp, y = DO, colour = type), alpha = 0.6)+
  scale_y_continuous(breaks = seq(0,150,25), limits = c(0,150))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2022-09-20','2023-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2022/2023")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank())+
  #legend.position = "none") +
  labs(y = "Dissolved Oxygen (%)", x = "Depth (m)")+
  scale_colour_manual(label = c("Observed", "Predicted"), values = c("lightgrey","#3399FF"))+
  guides(color = guide_legend(override.aes = list(size = 3))) 


DOForest_plot_2021_depth<-ggplot()+
  geom_point(data = DO_data, aes(x = hr, y = DO_m,colour = depth_m),  alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,150,25), limits = c(0,150))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  #ggtitle("Winter 2020/2021")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.position = "none") +
  labs(y = "Dissolved Oxygen (%)", x = "Depth (m)")+
  scale_colour_viridis_c(option="magma", trans = 'reverse' )

DOForest_plot_2022_depth<-ggplot()+
  geom_point(data = DO_data, aes(x = hr, y = DO_m,colour = depth_m),  alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,150,25), limits = c(0,150))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  #ggtitle("Winter 2021/2022")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.position = "none") +
  labs(y = "Dissolved Oxygen (%)", x = "Depth (m)")+
  scale_colour_viridis_c(option="magma", trans = 'reverse' )

DOForest_plot_2023_depth<-ggplot()+
  geom_point(data = DO_data, aes(x = hr, y = DO_m,colour = depth_m),  alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,150,25), limits = c(0,150))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2022-09-20','2023-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  #ggtitle("Winter 2022/2023")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13)) +
  labs(y = "Dissolved Oxygen (%)", x = "Depth (m)", colour = "Depth (m)")+
  scale_colour_viridis_c(option="magma", trans = 'reverse')

(DOForest_plot_2021|DOForest_plot_2022|DOForest_plot_2023)/(DOForest_plot_2021_depth|DOForest_plot_2022_depth|DOForest_plot_2023_depth)

#### Heart Rate Trends ####

tempvHR.lmer<-lmer(log(calchr)~(temp)+(1|tag), data = df_calcHR) #this model performs best 

summary(tempvHR.lmer)
Anova(tempvHR.lmer)
r.squaredGLMM(tempvHR.lmer)

d <- broom.mixed::tidy(tempvHR.lmer)

tempHR <- 
  tibble(
    temp = (seq(0.1, 30, 0.1))) %>%
  mutate(
    lmer_m = tempvHR.lmer@beta[c(1,0)] + 
      tempvHR.lmer@beta[c(2,0)] * (temp), 
    lower = (tempvHR.lmer@beta[c(1,0)] - (d$std.error[c(1,0)])*1.96) + 
      (tempvHR.lmer@beta[c(2,0)] - (d$std.error[c(2,0)])*1.96) * (temp),
    upper = (tempvHR.lmer@beta[c(1,0)] + (d$std.error[c(1,0)]*1.96)) + 
      (tempvHR.lmer@beta[c(2,0)] + (d$std.error[c(2,0)]*1.96)) * (temp)
  )

tempHR$lmer_m<-exp(tempHR$lmer_m)
tempHR$upper<-exp(tempHR$upper)
tempHR$lower<-exp(tempHR$lower)

#plotting relationships

hr.temp.gg<-ggplot() +
  geom_point(data = df_calcHR,aes(x=temp,y=calchr), fill="lightgrey", colour = "black", size = 0.8, alpha = 0.3) +
  geom_line(data = tempHR, aes(x = temp, y = lmer_m), colour = "red", size = 1) + 
  geom_ribbon(data = tempHR, aes(ymin = lower, ymax = upper, x = temp, y = lmer_m), fill = "salmon", size = 0.8, alpha = 0.6) +
  # annotate(geom="text", x=1.55, y=100, label = "n = 21", size = 5)+
  annotate(geom="text", x=(4), y=100, label = "italic(R^2)[Marginal] == 0.58", parse = TRUE, size = 5)+
  annotate(geom="text", x=(4.75), y=90, label = "italic(R^2)[Conditional] == 0.69", parse = TRUE, size = 5)+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,20))+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,30,5), limits = c(0,30))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"))+
  labs(x = " Temperature (°C)", y = "Heart Rate (bpm)")


hr.time.gg<-ggplot() + 
  geom_point(data = datHR,aes(x=(as.Date(timestamp)),y=calchr),  fill="lightgrey", colour = "black", alpha = 0.3, size=0.8)+
  geom_smooth(data = datHR,aes(x=(as.Date(timestamp)),y=calchr, colour = "Heart Rate"), method = "gam", se = FALSE, size = 0.8, alpha = 1)+
  geom_smooth(data = datHR,aes(x=(as.Date(timestamp)),y=(temp*3.3333), colour = "Temperature"), method = "gam", se = FALSE, size = 0.8, alpha = 1)+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,20),
                     sec.axis = sec_axis(~./3.3333, breaks = seq(0,30,5),name="Temperature (°C) \n "))+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 week",
               date_labels = "%b")+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 20.5),
        axis.text = element_text(size = 20, colour = "black"),
        legend.position = c(0.15,0.9), 
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title.y.left =element_blank(),
        axis.text.y.left = element_blank(),
        axis.text.x = element_text(angle = 90))+
  labs(colour = "Fish", x = " \n Date", y = "Heart Rate (bpm)")+
  #annotate(geom="text", x= (as.Date("2020-11-23")), y=100, label = "n = 13", size = 4.1)+
  scale_colour_manual(name = "Legend", values = c("green3","blue")) 

hr.temp.gg|hr.time.gg

#Time Trends in Heart Rates

hr.temp.2021<-ggplot(hr_temp_trends_2021, aes(x=hour,y=temp_m)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.6, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,30,5), limits = c(0,30))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-10-15','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2020/2021")+
  theme(plot.title = element_text(size = 13),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = "Temperature (°C)") +
  scale_colour_manual(labels = c("Fish 1","Fish 2", "Fish 3", "Fish 4", "Fish 5", "Fish 6", "Fish 7", "Fish 8", "Fish 9",
                                 "Fish 10", "Fish 11", "Fish 12", "Fish 13", "Fish 14"), 
                      values = c("violetred", "#800000FF","#725663FF", "#D49464FF","#FFB547FF", "lightblue","#0072B2", "#5B8FA8FF", "cyan4","darkgreen", 
                                 "chartreuse3", "orangered2", "darkslategrey","#660099"))

hr.temp.2022<-ggplot(hr_temp_trends_2022, aes(x=hour,y=temp_m)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.6, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,30,5), limits = c(0,30))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2022-10-15','2023-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2022/2023")+
  theme(plot.title = element_text(size = 13),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = "Temperature (°C)") +
  scale_colour_manual(labels = c("Fish 15", "Fish 16", "Fish 17", "Fish 18", "Fish 19", "Fish 20", "Fish 21"), 
                      values = c("#CC79A7","burlywood", "#000000","#009E73","#D55E00", "#E69F00","#56B4E9"))

hr.time.2021<-ggplot(hr_trends_2021, aes(x=hour,y=calchr_m)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.6, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,90,20), limits = c(0,90))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-10-15','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  #ggtitle("Winter 2020/2021 - HRT")+
  theme(plot.title = element_text(size = 13),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"), 
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = "Heart Rate (bpm)") +
  scale_colour_manual(labels = c("Fish 1","Fish 2", "Fish 3", "Fish 4", "Fish 5", "Fish 6", "Fish 7", "Fish 8", "Fish 9",
                                 "Fish 10", "Fish 11", "Fish 12", "Fish 13", "Fish 14"), 
                      values = c("violetred", "#800000FF","#725663FF", "#D49464FF","#FFB547FF", "lightblue","#0072B2", "#5B8FA8FF", "cyan4","darkgreen", 
                                 "chartreuse3", "orangered2", "darkslategrey","#660099"))

hr.time.2022<-ggplot(hr_trends_2022, aes(x=hour,y=calchr_m)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.6, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,90,20), limits = c(0,90))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2022-10-15','2023-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  #ggtitle("Winter 2022/2023 - HRT")+
  theme(plot.title = element_text(size = 13),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = expression(italic("M")*O["2"]~(mgO[2]~"\U00B7"~kg^{"-1"}~"\U00B7"~hr^{"-1"}))) +
  scale_colour_manual(labels = c("Fish 15", "Fish 16", "Fish 17", "Fish 18", "Fish 19", "Fish 20", "Fish 21"), 
                      values = c("#CC79A7","burlywood", "#000000","#009E73","#D55E00", "#E69F00","#56B4E9"))

(hr.temp.2021|hr.temp.2022)/(hr.time.2021|hr.time.2022)

#Q10 measures
glimpse(datHR)
datHR$temp<-as.numeric(datHR$temp)

q20c<-datHR %>% 
  filter(temp > 19 & temp < 21) %>% 
  group_by(tag) %>% 
  summarise(temp_20 = mean(temp, na.rm = TRUE),
            calchr_20 = mean(calchr, na.rm = TRUE))

q15c<-datHR %>% 
  filter(temp > 14 & temp < 16) %>% 
  group_by(tag) %>% 
  summarise(temp_15 = mean(temp, na.rm = TRUE),
            calchr_15 = mean(calchr, na.rm = TRUE)) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .))

q10c<-datHR %>% 
  filter(temp > 9 & temp < 11) %>% 
  group_by(tag) %>% 
  summarise(temp_10 = mean(temp, na.rm = TRUE),
            calchr_10 = mean(calchr, na.rm = TRUE)) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 


q5c<-datHR %>% 
  filter(temp > 4 & temp < 6) %>% 
  group_by(tag) %>% 
  summarise(temp_5 = mean(temp, na.rm = TRUE),
            calchr_5 = mean(calchr, na.rm = TRUE)) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

Q10_df<-q5c %>% 
  left_join(q10c, by = "tag") %>% 
  left_join(q15c, by = "tag") %>% 
  left_join(q20c, by = "tag") %>% 
  mutate(q10v5 = ((calchr_5/calchr_10)^(10/(temp_5 - temp_10)))) %>% 
  mutate(q15v5 = ((calchr_5/calchr_15)^(10/(temp_5 - temp_15)))) %>% 
  mutate(q20v5 = ((calchr_5/calchr_20)^(10/(temp_5 - temp_20)))) %>% 
  mutate(q20v10 = ((calchr_10/calchr_20)^(10/(temp_10 - temp_20))))


####Assessing Timetrends in Data####

df_trends<-df_hr_m %>% 
  mutate(id = (if_else(id == "Tag_1", "1",
                       if_else(id == "Tag_2", "2",
                               if_else(id == "Tag_3", "3",
                                       if_else(id == "Tag_4", "4",
                                               if_else(id == "Tag_5", "5",
                                                       if_else(id == "Tag_10", "6",
                                                               if_else(id == "Tag_11", "7",
                                                                       if_else(id == "Tag_12", "8",
                                                                               if_else(id == "Tag_13", "9",
                                                                                       if_else(id == "Tag_14", "10",
                                                                                               if_else(id == "Tag_15", "11",
                                                                                                       if_else(id == "Tag_16", "12",
                                                                                                               if_else(id == "Tag_18", "13",
                                                                                                                       if_else(id == "Tag_19", "14",
                                                                                                                               if_else(id == "Tag_20", "15",
                                                                                                                                       if_else(id == "Tag_1_2022", "16",
                                                                                                                                               if_else(id == "Tag_2_2022", "17",
                                                                                                                                                       if_else(id == "Tag_3_2022", "18",
                                                                                                                                                               if_else(id == "Tag_4_2022", "19",
                                                                                                                                                                       if_else(id == "Tag_5_2022", "20",
                                                                                                                                                                               if_else(id == "Tag_12_2022", "21",
                                                                                                                                                                                       if_else(id == "Tag_13_2022", "22",
                                                                                                                                                                                               if_else(id == "Tag_14_2022", "23",
                                                                                                                                                                                                       if_else(id == "Tag_11_2022", "24",
                                                                                                                                                                                                               if_else(id == "Tag_201033_2022", "25",
                                                                                                                                                                                                                       NA_character_)))))))))))))))))))))))))))

df_trends$id<-as.integer(df_trends$id)

df_trends_2021<-df_trends %>% 
  filter(hour < "2021-07-20 00:00:00")

df_trends_2021<-df_trends_2021[!is.na(df_trends_2021$ODBA_ms2),]

df_trends_2022<-df_trends %>% 
  filter(hour > "2021-07-20 00:00:00" & hour < "2022-09-01 00:00:00")

df_trends_2022<-df_trends_2022[!is.na(df_trends_2022$ODBA_ms2),]

#predicting DO based on fish values

DO_fish<-df_trends %>% 
  rename(hr = hour)

DO_fish$date<-as.numeric(strftime(DO_fish$hr, format = "%j"))

DO_fish$DO_predict <- predict(DOForest,DO_fish)

DO_fish_2021<-DO_fish %>% 
  filter(hr < "2021-07-20 00:00:00")

DO_fish_2022<-DO_fish %>% 
  filter(hr > "2021-07-20 00:00:00")

DO_fish_2022<-DO_fish_2022%>%
  mutate(id = if_else(id == 23, NA_integer_, 
                      if_else(id == 13, NA_integer_, id))) %>% 
  drop_na(id)#error with 23

#converting DO to mg/L

# https://www.waterontheweb.org/under/waterquality/oxygen.html - formula here

DO_fish_2021$C<-exp(7.7117 - 1.31403*log(DO_fish_2021$temp_m + 45.93))

DO_fish_2021$Pw<-11.8571-(3840.7/DO_fish_2021$temp_m)-(216961/(DO_fish_2021$temp_m^2))

DO_fish_2021$O<-0.000975-((1.426*10^(-5))*DO_fish_2021$temp_m)+((6.436*10^(-8))*DO_fish_2021$temp_m)

DO_fish_2021$Cp<-(DO_fish_2021$C*DO_fish_2021$bar)*
  (((1-DO_fish_2021$Pw/DO_fish_2021$bar)*(1-DO_fish_2021$O*DO_fish_2021$bar))/
     ((1-DO_fish_2021$Pw)*(1-DO_fish_2021$O)))

DO_fish_2021$DO_mgL<-(DO_fish_2021$DO_predict*DO_fish_2021$Cp)/100

DO_fish_2022$bar<-(DO_fish_2022$depth_m*997.05*9.80665)/100000 + 1.013

max(DO_fish_2022$bar,na.rm = TRUE)

max(DO_fish_2022$depth_m, na.rm = TRUE)

DO_fish_2022<-DO_fish_2022 %>% 
  drop_na(DO_predict)
glimpse(DO_fish_2022)

DO_fish_2022$C<-exp(7.7117 - 1.31403*log(DO_fish_2022$temp_m + 45.93))

DO_fish_2022$Pw<-11.8571-(3840.7/DO_fish_2022$temp_m)-(216961/(DO_fish_2022$temp_m^2))

DO_fish_2022$O<-0.000975-((1.426*10^(-5))*DO_fish_2022$temp_m)+((6.436*10^(-8))*DO_fish_2022$temp_m)

DO_fish_2022$Cp<-(DO_fish_2022$C*DO_fish_2022$bar)*
  (((1-DO_fish_2022$Pw/DO_fish_2022$bar)*(1-DO_fish_2022$O*DO_fish_2022$bar))/
     ((1-DO_fish_2022$Pw)*(1-DO_fish_2022$O)))

DO_fish_2022$DO_mgL<-(DO_fish_2022$DO_predict*DO_fish_2022$Cp)/100


#plotting data

table(DO_fish_2021$id)

cols = c("1" = "violetred", "2" = "#800000FF", "3" = "#725663FF", "4" = "#D49464FF","5" = "#FFB547FF", "6" ="lightblue","7" ="#0072B2", "8" ="#5B8FA8FF","9" = "cyan4","10" ="darkgreen", 
         "11" = "chartreuse3", "12" = "orangered2", 
         "13" ="#CC79A7","14" ="burlywood","15" = "#000000","16" ="#009E73","17" ="#D55E00", "18" ="#660099","19" ="#E69F00",
         "20" ="#56B4E9", "21" ="#339900", "22" ="steelblue", "23" ="#FF6666", "24" ="#33ff99","25" ="#FF3333")

DO.time.2021<-ggplot(DO_fish_2021,
                     aes(x=hr,y=DO_mgL)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,15,3), limits = c(0,15))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "\n Date", y = expression("Dissolved Oxygen (mg L"^{-1}*")")) +
  scale_colour_manual(values = cols)

DO.time.2022<-ggplot(DO_fish_2022,
                     aes(x=hr,y=DO_mgL)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,15,3), limits = c(0,15))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "\nDate", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"}))) +
  scale_colour_manual(values = cols)


DO.time.2021|DO.time.2022

temp.time.2021<-ggplot(df_trends_2021, aes(x=hour,y=temp_m)) + 
  geom_point(colour = "lightgrey", size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,30,5), limits=c(0,30)) + 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2020/2021")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = "Temperature (°C)")+
  scale_colour_manual(values = cols)

table(df_trends_2022$id)

temp.time.2022<-ggplot(df_trends_2022, aes(x=hour,y=temp_m)) + 
  geom_point(colour = "lightgrey", size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,30,5), limits=c(0,30)) + 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2021/2022")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = "Temperature (°C)")+
  scale_colour_manual(values = cols)

max(df_trends_2021$ODBA_ms2)

table(df_trends_2021$id)

ODBA.time.2021<-ggplot(df_trends_2021, aes(x=hour,y=ODBA_ms2)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"}))) +
  scale_colour_manual(values = cols)

table(df_trends_2022$id)

ODBA.time.2022<-ggplot(df_trends_2022, aes(x=hour,y=ODBA_ms2)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"}))) +
  scale_colour_manual(values = cols)


Depth.time.2021<-ggplot(df_trends_2021,
                        aes(x=hour,y=depth_m)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_reverse(limits = c(6.2, 0), breaks = seq(0,6,1))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "\n Date", y = "Depth (m)") +
  scale_colour_manual(values = cols)

Depth.time.2022<-ggplot(df_trends_2022,
                        aes(x=hour,y=depth_m)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_reverse(limits = c(6.2, 0), breaks = seq(0,6,1))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  # ggtitle("Winter 2021/2022")+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "\nDate", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"}))) +
  scale_colour_manual(values = cols)

DO.time.2021<-ggplot(DO_fish_2021,
                     aes(x=hr,y=DO_predict)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,120,20), limits = c(0,120))+
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "\n Date", y = "Dissolved Oxygen (%)") +
  scale_colour_manual(values = cols)

DO.time.2022<-ggplot(DO_fish_2022,
                     aes(x=hr,y=DO_predict)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,120,20), limits = c(0,120))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "\nDate", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"}))) +
  scale_colour_manual(values = cols)

speed.time.2021<-ggplot(df_trends_2021, aes(x=hour,y=speed_bls)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,2.0,0.5), limits = c(0,2.0))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2020/2021")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = expression("Swimming Speed"~(BL~"\U00B7"~sec^{"-1"}))) +
  scale_colour_manual(values = cols)

speed.time.2022<-ggplot(df_trends_2022, aes(x=hour,y=speed_bls)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.8, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,2.0,0.5), limits = c(0,2.0))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2021/2022")+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) +
  labs(colour = "Fish", x = "Date", y = expression("speed"~(m~"\U00B7"~sec^{"-2"}))) +
  scale_colour_manual(values = cols)


speed.time.2021|speed.time.2022

((temp.time.2021|temp.time.2022)/(DO.time.2021|DO.time.2022)/plot_spacer())|((speed.time.2021|speed.time.2022)/(ODBA.time.2021|ODBA.time.2022)/(Depth.time.2021|Depth.time.2022))

####Diel cycles####

dielmonth_df<-dat %>% 
  left_join(dat_averages, by = "timestamp") %>% 
  mutate(temp_m = if_else(is.na(temp_m), temp_mean, temp_m)) %>% 
  mutate(speed_bls = if_else(speed_bls < 0, 0, speed_bls)) %>% 
  group_by(id,hour2,month_character) %>%
  summarise(timestamp = mean(timestamp),
            temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            ODBA_m = mean(ODBA_ms2),
            speed_bls_m = mean(speed_bls,na.rm = TRUE),
            ODBA_sum = sum(ODBA_ms2_sum),
            ODBA_max = max(ODBA_ms2_max),
            ODBA_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(id,fill=list(temp_m=NA,depth_m=NA,ODBA_m=NA, speed_bls_m = NA, ODBA_sum=NA, ODBA_max=NA, ODBA_min=NA))

diel_df<-dat %>% 
  left_join(dat_averages, by = "timestamp") %>% 
  mutate(temp_m = if_else(is.na(temp_m), temp_mean, temp_m)) %>% 
  mutate(speed_bls = if_else(speed_bls < 0, 0, speed_bls)) %>% 
  group_by(id,hour2) %>%
  summarise(timestamp = mean(timestamp),
            temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            ODBA_m = mean(ODBA_ms2),
            speed_bls_m = mean(speed_bls,na.rm = TRUE),
            ODBA_sum = sum(ODBA_ms2_sum),
            ODBA_max = max(ODBA_ms2_max),
            ODBA_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(id,fill=list(temp_m=NA,depth_m=NA,ODBA_m=NA, speed_bls_m = NA, ODBA_sum=NA, ODBA_max=NA, ODBA_min=NA))

#plotting the average diel pattern over the entire measurement period

diel_df$hour2<-as.integer(diel_df$hour2)

speed.diel.plot<-
  ggplot(diel_df,aes(x=hour2,y=speed_bls_m, group = hour2)) + 
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.3,1.0,0.1), limits = c(0.3,1.0))+
  annotate("text", x = 0, y = 0.88, label = "ab", size = 3.5)+
  annotate("text", x = 1, y = 0.865, label = "a", size = 3.5)+
  annotate("text", x = 6, y = 0.90, label = "ab", size = 3.5)+
  annotate("text", x = 7, y = 0.91, label = "abc", size = 3.5)+
  annotate("text", x = 8, y = 0.925, label = "ab\ncd", size = 3.5)+
  annotate("text", x = 9, y = 0.915, label = "cde", size = 3.5)+
  annotate("text", x = 10, y = 0.915, label = "de", size = 3.5)+
  annotate("text", x = 11, y = 0.915, label = "e", size = 3.5)+
  annotate("text", x = 12, y = 0.92, label = "e", size = 3.5)+
  annotate("text", x = 13, y = 0.925, label = "bcd", size = 3.5)+
  annotate("text", x = 14, y = 0.93, label = "abc", size = 3.5)+
  annotate("text", x = 15, y = 0.94, label = "abc", size = 3.5)+
  annotate("text", x = 16, y = 0.95, label = "abc", size = 3.5)+
  annotate("text", x = 22, y = 0.92, label = "ab", size = 3.5)+
  annotate("text", x = 23, y = 0.91, label = "ab", size = 3.5)+
  theme_classic() + 
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title.x=element_blank())+
  # axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("Speed"~(BL~"\U00B7"~sec^{"-1"})))+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 5)

ODBA.diel.plot<-
  ggplot(diel_df,aes(x=hour2,y=ODBA_m, group = hour2)) + 
  geom_boxplot(fill="white") +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.18,0.34,0.04), limits = c(0.18,0.36))+
  annotate("text", x = 0, y = 0.295, label = "ab", size = 3.5)+
  annotate("text", x = 1, y = 0.29, label = "a", size = 3.5)+
  annotate("text", x = 6, y = 0.33, label = "g", size = 3.5)+
  annotate("text", x = 7, y = 0.325, label = "fg", size = 3.5)+
  annotate("text", x = 8, y = 0.315, label = "cf", size = 3.5)+
  annotate("text", x = 9, y = 0.305, label = "cd", size = 3.5)+
  annotate("text", x = 10, y = 0.305, label = "cde", size = 3.5)+
  annotate("text", x = 11, y = 0.315, label = "bc\nde", size = 3.5)+
  annotate("text", x = 12, y = 0.315, label = "bc\nde", size = 3.5)+
  annotate("text", x = 13, y = 0.315, label = "cdf", size = 3.5)+
  annotate("text", x = 14, y = 0.315, label = "cdf", size = 3.5)+
  annotate("text", x = 15, y = 0.315, label = "cf", size = 3.5)+
  annotate("text", x = 16, y = 0.32, label = "cf", size = 3.5)+
  annotate("text", x = 22, y = 0.305, label = "ab\nde", size = 3.5)+
  annotate("text", x = 23, y = 0.295, label = "abe", size = 3.5)+
  theme_classic() + 
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title.x=element_blank())+
  # axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 5)

ODBA.diel.plot

depth.diel.plot<-ggplot(diel_df %>% filter(!is.na(depth_m)), #removed NAs here so that ggplot runs
                        aes(x=hour2,y=(0-depth_m+1), group = hour2)) + 
  geom_boxplot(fill="white") +
  geom_jitter(alpha=0.3,width=0.3,height=0,size=2) +
  labs(x="Time",y="ODBA") +
  scale_x_continuous(breaks = seq(0,23,2)) + 
  scale_y_continuous(limits = c(-3,0),
                     breaks = seq(-3,0,1),
                     labels = c("4", "3", "2", "1"))+
  annotate("text", x = 0, y = -0.05, label = "a", size = 3.5)+
  annotate("text", x = 1, y = -0.05, label = "a", size = 3.5)+
  annotate("text", x = 6, y = -0.3, label = "bcd", size = 3.5)+
  annotate("text", x = 7, y = -0.3, label = "d", size = 3.5)+
  annotate("text", x = 8, y = -0.3, label = "cd", size = 3.5)+
  annotate("text", x = 9, y = -0.3, label = "bcd", size = 3.5)+
  annotate("text", x = 10, y = -0.3, label = "bcd", size = 3.5)+
  annotate("text", x = 11, y = -0.3, label = "bcd", size = 3.5)+
  annotate("text", x = 12, y = -0.3, label = "bcd", size = 3.5)+
  annotate("text", x = 13, y = -0.3, label = "bcd", size = 3.5)+
  annotate("text", x = 14, y = -0.3, label = "bcd", size = 3.5)+
  annotate("text", x = 15, y = -0.3, label = "abc", size = 3.5)+
  annotate("text", x = 16, y = -0.3, label = "ab", size = 3.5)+
  annotate("text", x = 22, y = -0.05, label = "a", size = 3.5)+
  annotate("text", x = 23, y = -0.05, label = "a", size = 3.5)+
  theme_classic() + 
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black")) +
  labs(x = "Time (hr)", y = "Depth (m)")+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 5)

depth.diel.plot

speed.diel.plot/ODBA.diel.plot/depth.diel.plot

#plotting the average diel period by month 

#ODBA by month

table(dielmonth_df$month_character)

diel.Sept<-dielmonth_df %>% 
  filter(month_character == "K_Sep")

ODBA.dielSept.plot<-
  ggplot() + 
  geom_boxplot(data=diel.Sept,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Sept,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Sept,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=1.4, y=0.3, label = "n = 9 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.15,0.3,0.05), limits = c(0.15,0.3))+
  annotate("text", x = 0, y = 3.05, label = "a")+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("September")

#October

table(dielmonth_df$month_character)

diel.Oct<-dielmonth_df %>% 
  filter(month_character == "L_Oct")

ODBA.dielOct.plot<-
  ggplot() + 
  geom_boxplot(data=diel.Oct,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Oct,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Oct,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=1.5, y=0.3, label = "n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.15,0.3,0.05), limits = c(0.15,0.3))+
  annotate("text", x = 0, y = 3.05, label = "a")+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("October")

#November 
diel.Nov<-dielmonth_df %>% 
  filter(month_character == "A_Nov")

ODBA.dielNov.plot<-
  ggplot() + 
  geom_boxplot(data=diel.Nov,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Nov,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Nov,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=3.75, y=0.3, label = "n = 5 W1; n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.15,0.3,0.05), limits = c(0.15,0.3))+
  annotate("text", x = 0, y = 3.05, label = "a")+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("November")

#December

diel.Dec<-dielmonth_df %>% 
  filter(month_character == "B_Dec")

ODBA.dielDec.plot<-
  ggplot() + 
  geom_boxplot(data=diel.Dec,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Dec,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Dec,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=3.75, y=0.3, label = "n = 5 W1; n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.15,0.3,0.05), limits = c(0.15,0.3))+
  theme_classic()+
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("December")

#January

diel.Jan<-dielmonth_df %>% 
  filter(month_character == "C_Jan")

ODBA.dielJanplot<- 
  ggplot() + 
  geom_boxplot(data=diel.Jan,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Jan,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Jan,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.3, label = "n = 9 W1; n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.15,0.3,0.05), limits = c(0.15,0.3))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size =19),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("January")

#February
diel.Feb<-dielmonth_df %>% 
  filter(month_character == "D_Feb")

ODBA.dielFebplot<-
  ggplot() + 
  geom_boxplot(data=diel.Feb,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Feb,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Feb,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.3, label = "n = 12 W1; n = 13 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.15,0.3,0.05), limits = c(0.15,0.3))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("February")

#March

diel.Mar<-dielmonth_df %>% 
  filter(month_character == "E_Mar")

ODBA.dielMarplot<- 
  ggplot() + 
  geom_boxplot(data=diel.Mar,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Mar,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Mar,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.3, label = "n = 12 W1; n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.15,0.3,0.05), limits = c(0.15,0.3))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("March")

#April

diel.Apr<-dielmonth_df %>% 
  filter(month_character == "F_Apr")

ODBA.dielAprplot<-
  ggplot() + 
  geom_boxplot(data=diel.Apr,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Apr,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Apr,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.3, label = "n = 12 W1; n = 6 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.15,0.3,0.05), limits = c(0.15,0.3))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("April")


#May

diel.May<-dielmonth_df %>% 
  filter(month_character == "G_May")

ODBA.dielMayplot<- 
  ggplot() + 
  geom_boxplot(data=diel.May,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.May,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.May,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=3.75, y=0.51, label = "n = 12 W1; n = 1 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.2,0.5,0.1), limits = c(0.2,0.51))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("May")

#June

diel.Jun<-dielmonth_df %>% 
  filter(month_character == "H_Jun")

ODBA.dielJunplot<-
  ggplot() + 
  geom_boxplot(data=diel.Jun,aes(x=hour2,y=ODBA_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Jun,aes(x=as.numeric(hour2),y=ODBA_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Jun,aes(x=hour2,y=ODBA_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=3.75, y=0.51, label = "n = 6 W1; n = 1 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0.2,0.5,0.1), limits = c(0.2,0.51))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank()) +
  labs(x = "Time (hr)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("June")

#combining

fig_ODBAdielmonth <- patchworkGrob ((ODBA.dielSept.plot | ODBA.dielOct.plot)/
                                      (ODBA.dielNov.plot | ODBA.dielDec.plot)/
                                      (ODBA.dielJanplot | ODBA.dielFebplot)/ 
                                      (ODBA.dielMarplot | ODBA.dielAprplot)/
                                      (ODBA.dielMayplot | ODBA.dielJunplot))

annotate_figure(fig_ODBAdielmonth, left = text_grob(expression("ODBA"~(m~"\U00B7"~sec^{"-2"})), color = "black", size = 14, rot = 90),
                bottom = text_grob("Time (hr)", color = "black", size = 14))


#Swimming Speed by Month
#Sept
diel.Sept<-dielmonth_df %>% 
  filter(month_character == "K_Sep")

speed.dielSept.plot<-
  ggplot() + 
  geom_boxplot(data=diel.Sept,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Sept,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Sept,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=1.4, y=0.3, label = "n = 9 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("September")

#October

diel.Oct<-dielmonth_df %>% 
  filter(month_character == "L_Oct")

speed.dielOct.plot<-
  ggplot() + 
  geom_boxplot(data=diel.Oct,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Oct,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Oct,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=1.4, y=0.3, label = "n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("October")

#November

diel.Nov<-dielmonth_df %>% 
  filter(month_character == "A_Nov")

speed.dielNov.plot<-
  ggplot() + 
  geom_boxplot(data=diel.Nov,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Nov,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Nov,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=3.75, y=0.3, label = "n = 5 W1; n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("November")

#December

diel.Dec<-dielmonth_df %>% 
  filter(month_character == "B_Dec")

speed.dielDec.plot<-
  ggplot() + 
  geom_boxplot(data=diel.Dec,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Dec,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Dec,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=3.75, y=0.3, label = "n = 5 W1; n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic()+
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("December")

#January

diel.Jan<-dielmonth_df %>% 
  filter(month_character == "C_Jan")

speed.dielJanplot<- 
  ggplot() + 
  geom_boxplot(data=diel.Jan,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Jan,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Jan,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.3, label = "n = 9 W1; n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size =19),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("January")

#February
diel.Feb<-dielmonth_df %>% 
  filter(month_character == "D_Feb")

speed.dielFebplot<-
  ggplot() + 
  geom_boxplot(data=diel.Feb,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Feb,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Feb,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.3, label = "n = 12 W1; n = 13 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("February")

#March

diel.Mar<-dielmonth_df %>% 
  filter(month_character == "E_Mar")

speed.dielMarplot<- 
  ggplot() + 
  geom_boxplot(data=diel.Mar,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Mar,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Mar,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.3, label = "n = 12 W1; n = 12 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("March")

#April

diel.Apr<-dielmonth_df %>% 
  filter(month_character == "F_Apr")

speed.dielAprplot<-
  ggplot() + 
  geom_boxplot(data=diel.Apr,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Apr,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Apr,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.3, label = "n = 12 W1; n = 6 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("April")


#May

diel.May<-dielmonth_df %>% 
  filter(month_character == "G_May")

speed.dielMayplot<- 
  ggplot() + 
  geom_boxplot(data=diel.May,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.May,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.May,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=3.75, y=0.51, label = "n = 12 W1; n = 1 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("May")

#June

diel.Jun<-dielmonth_df %>% 
  filter(month_character == "H_Jun")

speed.dielJunplot<-
  ggplot() + 
  geom_boxplot(data=diel.Jun,aes(x=hour2,y=speed_bls_m, group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Jun,aes(x=as.numeric(hour2),y=speed_bls_m),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Jun,aes(x=hour2,y=speed_bls_m, group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=3.75, y=0.51, label = "n = 6 W1; n = 1 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,23,2)) +
  scale_y_continuous(breaks = seq(0,1.6,0.25), limits = c(0,1.6))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank()) +
  labs(x = "Time (hr)", y = expression("speed"~(m~"\U00B7"~sec^{"-2"})))+
  ggtitle("June")

#combining

fig_speeddielmonth <- patchworkGrob (((speed.dielNov.plot | speed.dielDec.plot)/
                                        (speed.dielJanplot | speed.dielFebplot)/ 
                                        (speed.dielMarplot | speed.dielAprplot)/
                                        (speed.dielMayplot | speed.dielJunplot)))

annotate_figure(fig_speeddielmonth, left = text_grob(expression("speed"~(m~"\U00B7"~sec^{"-2"})), color = "black", size = 14, rot = 90),
                bottom = text_grob("Time (hr)", color = "black", size = 14))


#Depth by month

#September 

diel.Sept_depth<-dielmonth_df_depth %>% 
  filter(month_character == "K_Sep")

depth.dielSept.plot<-
  ggplot() +
  geom_boxplot(data=diel.Sept_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Sept_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.Sept_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=1.4, y=0.25, label = "n = 6 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("September")

#October 

diel.Oct_depth<-dielmonth_df_depth %>% 
  filter(month_character == "L_Oct")

depth.dielOct.plot<-
  ggplot() +
  geom_boxplot(data=diel.Oct_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Oct_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.Oct_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=1.4, y=0.25, label = "n = 8 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("October")


#November 

diel.Nov_depth<-dielmonth_df_depth %>% 
  filter(month_character == "A_Nov")

depth.dielNov.plot<-
  ggplot() +
  geom_boxplot(data=diel.Nov_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Nov_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.Nov_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour = "black", size=0.7) +
  annotate(geom="text", x=4, y=0.25, label = "n = 5 W1; n = 8 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("November")


#December
diel.Dec_depth<-dielmonth_df_depth %>% 
  filter(month_character == "B_Dec")

depth.dielDec.plot<-
  ggplot() +
  geom_boxplot(data=diel.Dec_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Dec_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.Dec_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour = "black", size = 0.7) +
  annotate(geom="text", x=4, y=0.25, label = "n = 5 W1; n = 8 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("December")

#January

diel.Jan_depth<-dielmonth_df_depth %>% 
  filter(month_character == "C_Jan")

depth.dielJanplot<-
  ggplot() +
  geom_boxplot(data=diel.Jan_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Jan_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.Jan_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour ="black", size = 0.7) +
  annotate(geom="text", x=4, y=0.25, label = "n = 9 W1; n = 8 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("January")

#February

diel.Feb_depth<-dielmonth_df_depth %>% 
  filter(month_character == "D_Feb")

depth.dielFebplot<- 
  ggplot() +
  geom_boxplot(data=diel.Feb_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Feb_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.Feb_depth,aes(x=hour2,y=(0-depth_m), group = hour2),colour = "black", size = 0.7) +
  annotate(geom="text", x=4, y=0.25, label = "n = 11 W1; n = 8 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("February")

#March

diel.Mar_depth<-dielmonth_df_depth %>% 
  filter(month_character == "E_Mar")

depth.dielMarplot<- 
  ggplot() +
  geom_boxplot(data=diel.Mar_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Mar_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.Mar_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour = "black", size = 0.7) +
  annotate(geom="text", x=4, y=0.25, label = "n = 11 W1; n = 7 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("March")

#April

diel.Apr_depth<-dielmonth_df_depth %>% 
  filter(month_character == "F_Apr")

depth.dielAprplot<-
  ggplot() +
  geom_boxplot(data=diel.Apr_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Apr_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.Apr_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour = "black", size = 0.7) +
  annotate(geom="text", x=4, y=0.25, label = "n = 11 W1; n = 2 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"), 
        axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("April")

#May

diel.May_depth<-dielmonth_df_depth %>% 
  filter(month_character == "G_May")

depth.dielMayplot<-
  ggplot() +
  geom_boxplot(data=diel.May_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.May_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7,colour = "red") +
  geom_point(data=diel.May_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour = "black", size = 0.7) +
  annotate(geom="text", x=3.75, y=0.25, label = "n = 10 W1; n = 1 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("May")

#June

diel.Jun_depth<-dielmonth_df_depth %>% 
  filter(month_character == "H_Jun")

depth.dielJunplot<-
  ggplot() +
  geom_boxplot(data=diel.Jun_depth,aes(x=hour2,y=(0-depth_m), group = hour2),fill="white", outlier.shape = NA) +
  geom_smooth(data=diel.Jun_depth,aes(x=as.numeric(hour2),y=(0-depth_m)),method="loess",se = TRUE, fill = "salmon", alpha = 0.3, size = 0.7, colour = "red") +
  geom_point(data=diel.Jun_depth,aes(x=hour2,y=(0-depth_m), group = hour2), colour = "black", size = 0.7) +
  annotate(geom="text", x=3.75, y=0.25, label = "n = 4 W1; n = 1 W2", size = 4)+
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(limits = c(-5, 0.25),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  theme_classic() + 
  theme(plot.title = element_text(size = 14),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "\nTime (hr)", y = "Depth (m)\n")+
  ggtitle("June")

#combining
fig_depthdielmonth <- patchworkGrob (((depth.dielSept.plot | depth.dielOct.plot)/
                                        (depth.dielNov.plot | depth.dielDec.plot)/
                                        (depth.dielJanplot | depth.dielFebplot)/ 
                                        (depth.dielMarplot | depth.dielAprplot)/
                                        (depth.dielMayplot | depth.dielJunplot)))

annotate_figure(fig_depthdielmonth, left = text_grob("Depth (m)", color = "black", size = 14, rot = 90),
                bottom = text_grob("Time (hr)", color = "black", size = 14))


#### Histograms of Monthly Behaviour ####

#Nov 

dat.Nov<-dat %>% 
  filter(month_character == "A_Nov")

Nov.ODBA.hist<-ggplot(dat.Nov, aes (x = ODBA_ms2, y = ..count../1000)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "white")+
  geom_vline(xintercept = mean(dat.Nov$ODBA_ms2), col = "red", lwd = 0.7)+
  annotate(geom="text", x=0.91, y=285000/1000, label = "Mean = 0.21", size = 3.5, hjust = 1)+
  # annotate(geom="text", x= 0.91, y=260000/1000, label = "n = 5 W1\nn = 12 W2\nMean = 0.21", size = 3.5, hjust = 1)+
  #annotate(geom="text", x=0.91, y=225000, label =  "Mean = 0.212", size = 3.5, hjust = 1)+
  ggtitle("November")+
  scale_x_continuous(breaks = seq(0,1,0.2), limits=c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,300000,50000)/1000, limits = c(0, 300000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank())+
  labs(x = "\nAverage ODBA (g) min", y = "Count (Total Minutes)\n")

#Dec

dat.Dec<-dat %>% 
  filter(month_character == "B_Dec")

Dec.ODBA.hist<-ggplot(dat.Dec, aes (y = ..count../1000, x = ODBA_ms2)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "white")+
  ggtitle("December")+
  geom_vline(xintercept = mean(dat.Dec$ODBA_ms2), col = "red", lwd = 0.7)+
  annotate(geom="text", x=0.91, y=285000/1000, label = "Mean = 0.21", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=260000/1000, label = "n = 5 W1\nn = 12 W2\nMean = 0.21", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=225000/1000, label =  "Mean = 0.206", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,1,0.2), limits=c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,300000,50000)/1000, limits = c(0, 300000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.y=element_blank(),
        # axis.text.y = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank())+
  labs(x = "\nODBA (g)", y = "Count\n")

#Jan

dat.Jan<-dat %>% 
  filter(month_character == "C_Jan")

Jan.ODBA.hist<-ggplot(dat.Jan, aes (y = ..count../1000, x = ODBA_ms2)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "white")+
  ggtitle("January")+
  geom_vline(xintercept = mean(dat.Jan$ODBA_ms2), col = "red", lwd = 0.7)+
  annotate(geom="text", x=0.91, y=285000/1000, label = "Mean = 0.21", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=260000/1000, label = "n = 9 W1\nn = 12 W2\nMean = 0.21", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=225000/1000, label =  "Mean = 0.213", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,1,0.2), limits=c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,300000,50000)/1000, limits = c(0, 300000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank())+
  labs(x = "\nODBA (g)", y = "Count\n")

#Feb

dat.Feb<-dat %>% 
  filter(month_character == "D_Feb")

Feb.ODBA.hist<-ggplot(dat.Feb, aes (y = ..count../1000, x = ODBA_ms2)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "white")+
  ggtitle("February")+
  geom_vline(xintercept = mean(dat.Feb$ODBA_ms2), col = "red", lwd = 0.7)+
  annotate(geom="text", x=0.91, y=285000/1000, label = "Mean = 0.22", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=260000/1000, label = "n = 12 W1\nn = 13 W2\nMean = 0.22", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,1,0.2), limits=c(0,1), expand = c(0.005,0.005))+
  #annotate(geom="text", x=0.91, y=225000/1000, label =  "Mean = 0.217", size = 3.5, hjust = 1)+
  scale_y_continuous(breaks = seq(0,300000,50000)/1000, limits = c(0, 300000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.y=element_blank(),
        #axis.text.y = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank())+
  labs(x = "\nODBA (g)", y = "Count\n")

#Mar

dat.Mar<-dat %>% 
  filter(month_character == "E_Mar")

Mar.ODBA.hist<-ggplot(dat.Mar, aes (y = ..count../1000, x = ODBA_ms2)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "white")+
  ggtitle("March")+
  geom_vline(xintercept = mean(dat.Mar$ODBA_ms2), col = "red", lwd = 0.7)+
  annotate(geom="text", x=0.91, y=285000/1000, label = "Mean = 0.23", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=260000/1000, label = "n = 12 W1\nn = 12 W2\nMean = 0.23", size = 3.5, hjust = 1)+
  #annotate(geom="text", x=0.91, y=225000/1000, label =  "Mean = 0.227", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,1,0.2), limits=c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,300000,50000)/1000, limits = c(0, 300000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  #labs(x =  expression(atop("",paste("Average ODBA (g ⋅ min" ^-1, ")"))), y = "Count (Total Minutes)\n")
  labs(x = "\nODBA (g)", y = "Count (Total Minutes)\n")

#Apr

dat.Apr<-dat %>% 
  filter(month_character == "F_Apr")

Apr.ODBA.hist<-ggplot(dat.Apr, aes (y = ..count../1000, x = ODBA_ms2)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "white")+
  ggtitle("April")+
  geom_vline(xintercept = mean(dat.Apr$ODBA_ms2), col = "red", lwd = 0.7)+
  annotate(geom="text", x=0.91, y=190000/1000, label = "Mean = 0.33", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=260000/1000, label = "n = 12 W1\nn = 6 W2\nMean = 0.23", size = 3.5, hjust = 1)+
  #annotate(geom="text", x=0.91, y=225000/1000, label =  "Mean = 0.230", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,1,0.2), limits=c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,200000,50000)/1000, limits = c(0, 200000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.y=element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  #labs(x =  expression(atop("",paste("Average ODBA (g ⋅ min" ^-1, ")"))), y = "Count (Total Minutes)\n")
  labs(x = "\nODBA (g)", y = "Count (Total Minutes)\n")

#May

dat.May<-dat %>% 
  filter(month_character == "G_May")

May.ODBA.hist<-ggplot(dat.May, aes (y = ..count../1000, x = ODBA_ms2)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "white")+
  ggtitle("May")+
  geom_vline(xintercept = mean(dat.May$ODBA_ms2), col = "red", lwd = 0.7)+
  annotate(geom="text", x=0.91, y=57000/1000, label = "Mean = 0.32", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=260000/1000, label = "n = 12 W1\nn = 1 W2\nMean = 0.32", size = 3.5, hjust = 1)+
  #annotate(geom="text", x=0.91, y=225000/1000, label =  "Mean = 0.330", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,1,0.2), limits=c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,60000,10000)/1000, limits = c(0, 60000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title.y=element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  #labs(x =  expression(atop("",paste("Average ODBA (g ⋅ min" ^-1, ")"))), y = "Count (Total Minutes)\n")
  labs(x = "\nODBA (g)", y = "Count (Total Minutes)\n")

#Jun

dat.Jun<-dat %>% 
  filter(month_character == "H_Jun")

Jun.ODBA.hist<-ggplot(dat.Jun, aes (y = ..count../1000, x = ODBA_ms2)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "white")+
  ggtitle("June")+
  geom_vline(xintercept = mean(dat.Jun$ODBA_ms2), col = "red", lwd = 0.7)+
  annotate(geom="text", x=0.91, y=19000/1000, label = "Mean = 0.33", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=0.91, y=260000/1000, label = "n = 12 W1\nn = 1 W2\nMean = 0.33", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,1,0.2), limits=c(0,1), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,20000,5000)/1000, limits = c(0, 20000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.y=element_blank())+
  #labs(x =  expression(atop("",paste("Average ODBA (g ⋅ min" ^-1, ")"))), y = "Count (Total Minutes)\n")
  labs(x = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})), y = "Count (Total Minutes)\n")

#combining

require(grid)

OBDA.hist<-ggarrange(Nov.ODBA.hist+rremove("ylab")+rremove("xlab"), 
                     Dec.ODBA.hist+rremove("ylab")+rremove("xlab"),
                     Jan.ODBA.hist+rremove("ylab")+rremove("xlab"), 
                     Feb.ODBA.hist+rremove("ylab")+rremove("xlab"), 
                     Mar.ODBA.hist+rremove("ylab")+rremove("xlab"),
                     Apr.ODBA.hist+rremove("ylab")+rremove("xlab"),
                     May.ODBA.hist+rremove("ylab")+rremove("xlab"),
                     Jun.ODBA.hist+rremove("ylab")+rremove("xlab"),
                     ncol = 4, nrow = 2, widths = c(1.2,1,1,1,1,1,1,1))


annotate_figure(OBDA.hist, left = text_grob("Count (Total Minutes)", color = "black", size = 14, rot = 90),
                bottom = text_grob(expression("ODBA"~(m~"\U00B7"~sec^{"-2"})), color = "black", size = 14))

#swimming speed

dat.filled$speed_bls<-10^(1.3156+0.0277*dat.filled$temp_m+0.3423*log10(dat.filled$ODBA_ms2)-0.6505*log10(dat.filled$length))
dat.filled<-dat.filled %>% 
  mutate(speed_bls = if_else(speed_bls < 0, 0, speed_bls)) %>% 
  mutate(speed_bls = if_else(ODBA_ms2 < 0.1, 0, speed_bls))

#Nov 
dat.Nov<-dat.filled%>% 
  filter(month_character == "A_Nov")

Novspeed.hist<-ggplot(dat.Nov, aes (y = ..count../1000, x = speed_bls)) + 
  geom_histogram(binwidth = 0.1, color = "black", fill = "white", boundary = 0)+
  ggtitle("November")+
  geom_vline(xintercept = mean(dat.Nov$speed_bls), col = "red", lwd = 0.7)+
  annotate(geom="text", x=2.275, y=190000/1000, label = "Mean = 0.40", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=1.456, y=225000/1000, label =  "Mean = 0.206", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,2.5,0.5), limits=c(0,2.5), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,200000,50000)/1000, limits = c(0, 200000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "\nAverage ODBA (g) min", y = "Count (Total Minutes)\n")

#Dec

dat.Dec<-dat.filled%>% 
  filter(month_character == "B_Dec")

Decspeed.hist<-ggplot(dat.Dec, aes (y = ..count../1000, x = speed_bls)) + 
  geom_histogram(binwidth = 0.1, color = "black", fill = "white", boundary = 0)+
  ggtitle("December")+
  geom_vline(xintercept = mean(dat.Dec$speed_bls), col = "red", lwd = 0.7)+
  annotate(geom="text", x=2.275, y=380000/1000, label = "Mean = 0.34", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=1.456, y=225000/1000, label =  "Mean = 0.206", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,2.5,0.5), limits=c(0,2.5), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,400000,50000)/1000, limits = c(0, 400000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.y=element_blank(),
        #axis.text.y = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank())+
  labs(x = "\nODBA (g)", y = "Count\n")

#Jan

dat.Jan<-dat.filled%>% 
  filter(month_character == "C_Jan")

Janspeed.hist<-ggplot(dat.Jan, aes (y = ..count../1000, x = speed_bls)) + 
  geom_histogram(binwidth = 0.1, color = "black", fill = "white", boundary = 0)+
  ggtitle("January")+
  geom_vline(xintercept = mean(dat.Jan$speed_bls), col = "red", lwd = 0.7)+
  annotate(geom="text", x=2.275, y=475000/1000, label = "Mean = 0.35", size = 3.5, hjust = 1)+
  # annotate(geom="text", x=1.456, y=225000/1000, label =  "Mean = 0.213", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,2.5,0.5), limits=c(0,2.5), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,500000,50000)/1000, limits = c(0, 500000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank())+
  labs(x = "\nODBA (g)", y = "Count\n")

#Feb

dat.Feb<-dat.filled%>% 
  filter(month_character == "D_Feb")

Febspeed.hist<-ggplot(dat.Feb, aes (y = ..count../1000, x = speed_bls)) + 
  geom_histogram(binwidth = 0.1, color = "black", fill = "white", boundary = 0)+
  ggtitle("February")+
  geom_vline(xintercept = mean(dat.Feb$speed_bls), col = "red", lwd = 0.7)+
  annotate(geom="text", x=2.275, y=494000/1000, label = "Mean = 0.35", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,2.5,0.5), limits=c(0,2.5), expand = c(0.005,0.005))+
  #annotate(geom="text", x=1.456, y=225000/1000, label =  "Mean = 0.217", size = 3.5, hjust = 1)+
  scale_y_continuous(breaks = seq(0,500000,100000)/1000, limits = c(0, 520000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.y=element_blank(),
        #axis.text.y = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank())+
  labs(x = "\nODBA (g)", y = "Count\n")

#Mar

dat.Mar<-dat.filled%>% 
  filter(month_character == "E_Mar")

Marspeed.hist<-ggplot(dat.Mar, aes (y = ..count../1000, x = speed_bls)) + 
  geom_histogram(binwidth = 0.1, color = "black", fill = "white", boundary = 0)+
  ggtitle("March")+
  geom_vline(xintercept = mean(dat.Mar$speed_bls), col = "red", lwd = 0.7)+
  annotate(geom="text", x=2.275, y=494000/1000, label = "Mean = 0.37", size = 3.5, hjust = 1)+
  #annotate(geom="text", x=1.456, y=225000/1000, label =  "Mean = 0.227", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,2.5,0.5), limits=c(0,2.5), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,520000,100000)/1000, limits = c(0, 520000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  #labs(x =  expression(atop("",paste("Average ODBA (g ⋅ min" ^-1, ")"))), y = "Count (Total Minutes)\n")
  labs(x = "\nODBA (g)", y = "Total Minutes X 1000")

#Apr

dat.Apr<-dat.filled%>% 
  filter(month_character == "F_Apr")

Aprspeed.hist<-ggplot(dat.Apr, aes (y = ..count../1000, x = speed_bls)) + 
  geom_histogram(binwidth = 0.1, color = "black", fill = "white", boundary = 0)+
  ggtitle("April")+
  geom_vline(xintercept = mean(dat.Apr$speed_bls,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=2.275, y=142500/1000, label = "Mean = 0.56", size = 3.5, hjust = 1)+
  #annotate(geom="text", x=1.456, y=225000/1000, label =  "Mean = 0.230", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,2.5,0.5), limits=c(0,2.5), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,150000,25000)/1000, limits = c(0, 150000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.y=element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  #labs(x =  expression(atop("",paste("Average ODBA (g ⋅ min" ^-1, ")"))), y = "Count (Total Minutes)\n")
  labs(x = "\nODBA (g)", y = "Count (Total Minutes)\n")

#May

dat.May<-dat.filled%>% 
  filter(month_character == "G_May")

Mayspeed.hist<-ggplot(dat.May, aes (y = ..count../1000, x = speed_bls)) + 
  geom_histogram(binwidth = 0.1, color = "black", fill = "white", boundary = 0)+
  ggtitle("May")+
  geom_vline(xintercept = mean(dat.May$speed_bls,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=2.275, y=47500/1000, label = "Mean = 0.90", size = 3.5, hjust = 1)+
  #annotate(geom="text", x=1.456, y=225000/1000, label =  "Mean = 0.330", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,2.5,0.5), limits=c(0,2.5), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,50000,10000)/1000, limits = c(0, 50000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title.y=element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  #labs(x =  expression(atop("",paste("Average ODBA (g ⋅ min" ^-1, ")"))), y = "Count (Total Minutes)\n")
  labs(x = "\nODBA (g)", y = "Count (Total Minutes)\n")

#Jun

dat.Jun<-dat.filled%>% 
  filter(month_character == "H_Jun")

Junspeed.hist<-ggplot(dat.Jun, aes (y = ..count../1000, x = speed_bls)) + 
  geom_histogram(binwidth = 0.1, color = "black", fill = "white", boundary = 0)+
  ggtitle("June")+
  geom_vline(xintercept = mean(dat.Jun$speed_bls,na.rm =TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=2.275, y=13650/1000, label = "Mean = 1.24", size = 3.5, hjust = 1)+
  #annotate(geom="text", x=1.456, y=225000/1000, label =  "Mean = 0.358", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,2.5,0.5), limits=c(0,2.5), expand = c(0.005,0.005))+
  scale_y_continuous(breaks = seq(0,15000,2000)/1000, limits = c(0, 15000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.y=element_blank())+
  #labs(x =  expression(atop("",paste("Average ODBA (g ⋅ min" ^-1, ")"))), y = "Count (Total Minutes)\n")
  labs(x = expression("Speed"~(BL~"\U00B7"~sec^{"-1"})), y = "Count (Total Minutes)\n")

#combining

require(grid)

speed.hist<-ggarrange(Novspeed.hist+rremove("ylab")+rremove("xlab"), 
                      Decspeed.hist+rremove("ylab")+rremove("xlab"),
                      Janspeed.hist+rremove("ylab")+rremove("xlab"), 
                      Febspeed.hist+rremove("ylab")+rremove("xlab"), 
                      Marspeed.hist+rremove("ylab")+rremove("xlab"),
                      Aprspeed.hist+rremove("ylab")+rremove("xlab"),
                      Mayspeed.hist+rremove("ylab")+rremove("xlab"),
                      Junspeed.hist+rremove("ylab")+rremove("xlab"),
                      ncol = 4, nrow = 2, widths = c(1.2,1,1,1,1,1,1,1))


annotate_figure(speed.hist, left = text_grob("Count (Total Minutes)", color = "black", size = 14, rot = 90),
                bottom = text_grob(expression("Speed"~(BL~"\U00B7"~sec^{"-1"})), color = "black", size = 14))

#Depth 

#Nov 

Nov.d.hist<-ggplot(dat.Nov, aes (y = ..count../1000, x = depth_m)) + 
  geom_histogram(binwidth = 0.25, color = "black", fill = "white")+
  ggtitle("November")+
  geom_vline(xintercept = mean(dat.Nov$depth_m,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=6.4, y=47500/1000, label = "Mean = 3.04", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,7,1), limits=c(-0.25,7), expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(0,50000,10000)/1000, limits = c(0, 50000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "\nDepth (m)", y = "Count (Total Minutes)\n")

#Dec

Dec.d.hist<-ggplot(dat.Dec, aes (y = ..count../1000, x = depth_m)) + 
  geom_histogram(binwidth = 0.25, color = "black", fill = "white")+
  ggtitle("December")+
  geom_vline(xintercept = mean(dat.Dec$depth_m,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=6.4, y=47500/1000, label = "Mean = 3.73", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,7,1), limits=c(-0.25,7), expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(0,50000,10000)/1000, limits = c(0, 50000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "\nDepth (m)", y = "Count (Total Minutes)\n")

#Jan

Jan.d.hist<-ggplot(dat.Jan, aes (y = ..count../1000, x = depth_m)) + 
  geom_histogram(binwidth = 0.25, color = "black", fill = "white")+
  ggtitle("January")+
  geom_vline(xintercept = mean(dat.Jan$depth_m,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=6.4, y=85500/1000, label = "Mean = 3.82", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,7,1), limits=c(-0.25,7), expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(0,90000,20000)/1000, limits = c(0, 90000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "\nDepth (m)", y = "Count (Total Minutes)\n")

#Feb

Feb.d.hist<-ggplot(dat.Feb, aes (y = ..count../1000, x = depth_m)) + 
  geom_histogram(binwidth = 0.25, color = "black", fill = "white")+
  ggtitle("February")+
  geom_vline(xintercept = mean(dat.Feb$depth_m,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=6.4, y=85500/1000, label = "Mean = 3.17", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,7,1), limits=c(-0.25,7), expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(0,90000,20000)/1000, limits = c(0, 90000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "\nDepth (m)", y = "Count (Total Minutes)\n")

#Mar

Mar.d.hist<-ggplot(dat.Mar, aes (y = ..count../1000, x = depth_m)) + 
  geom_histogram(binwidth = 0.25, color = "black", fill = "white")+
  ggtitle("March")+
  geom_vline(xintercept = mean(dat.Mar$depth_m,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=6.4, y=85500/1000, label = "Mean = 2.19", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,7,1), limits=c(-0.25,7), expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(0,90000,20000)/1000, limits = c(0, 90000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "\nDepth (m)", y = "Count (Total Minutes)\n")

#Apr

Apr.d.hist<-ggplot(dat.Apr, aes (y = ..count../1000, x = depth_m)) + 
  geom_histogram(binwidth = 0.25, color = "black", fill = "white")+
  ggtitle("April")+
  geom_vline(xintercept = mean(dat.Apr$depth_m,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=6.4, y=57000/1000, label = "Mean = 0.99", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,7,1), limits=c(-0.25,7), expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(0,60000,10000)/1000, limits = c(0, 60000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = "\nDepth (m)", y = "Count (Total Minutes)\n")

#May

May.d.hist<-ggplot(dat.May, aes (y = ..count../1000, x = depth_m)) + 
  geom_histogram(binwidth = 0.25, color = "black", fill = "white")+
  ggtitle("May")+
  geom_vline(xintercept = mean(dat.May$depth_m,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=6.4, y=28500/1000, label = "Mean = 0.98", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,7,1), limits=c(-0.25,7), expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(0,30000,5000)/1000, limits = c(0, 30000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = "\nDepth (m)", y = "Count (Total Minutes)\n")

#Jun

Jun.d.hist<-ggplot(dat.Jun, aes (y = ..count../1000, x = depth_m)) + 
  geom_histogram(binwidth = 0.25, color = "black", fill = "white")+
  ggtitle("June")+
  geom_vline(xintercept = mean(dat.Jun$depth_m,na.rm = TRUE), col = "red", lwd = 0.7)+
  annotate(geom="text", x=6.4, y=9500/1000, label = "Mean = 0.84", size = 3.5, hjust = 1)+
  scale_x_continuous(breaks = seq(0,7,1), limits=c(-0.25,7), expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(0,10000,2000)/1000, limits = c(0, 10000)/1000, expand = c(0,0)) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title.y = element_blank())+
  #axis.text.y = element_blank())+
  labs(x = "Depth (m)", y = "Count (Total Minutes)\n")

#combining

require(grid)

depth.hist<-ggarrange(Nov.d.hist+ rremove("ylab"), 
                      Dec.d.hist, 
                      Jan.d.hist, 
                      Feb.d.hist, 
                      Mar.d.hist+rremove("ylab")+rremove("xlab"),
                      Apr.d.hist+rremove("xlab"), 
                      May.d.hist+rremove("xlab"), 
                      Jun.d.hist+rremove("xlab"), 
                      ncol = 4, nrow = 2, widths = c(1.2,1,1,1,1,1,1,1))

annotate_figure(depth.hist, left = text_grob("Count (Total Minutes)\n", color = "black", size = 13, rot = 90),
                bottom = text_grob("Depth (m)", color = "black", size = 13))


#combining all histograms

all.hist<-ggarrange(Novspeed.hist+rremove("ylab")+rremove("xlab"), 
                    Nov.ODBA.hist+rremove("ylab")+rremove("xlab"),
                    Nov.d.hist+rremove("ylab")+rremove("xlab"), 
                    Decspeed.hist+rremove("ylab")+rremove("xlab"),
                    Dec.ODBA.hist+rremove("ylab")+rremove("xlab"),
                    Dec.d.hist+rremove("ylab")+rremove("xlab"),
                    Janspeed.hist+rremove("ylab")+rremove("xlab"), 
                    Jan.ODBA.hist+rremove("ylab")+rremove("xlab"),
                    Jan.d.hist+rremove("ylab")+rremove("xlab"), 
                    Febspeed.hist+rremove("ylab")+rremove("xlab"),
                    Feb.ODBA.hist+rremove("ylab")+rremove("xlab"), 
                    Feb.d.hist+rremove("ylab")+rremove("xlab"), 
                    Marspeed.hist+rremove("ylab")+rremove("xlab"),
                    Mar.ODBA.hist+rremove("ylab")+rremove("xlab"),
                    Mar.d.hist+rremove("ylab")+rremove("xlab"),
                    Aprspeed.hist+rremove("ylab")+rremove("xlab"),
                    Apr.ODBA.hist+rremove("ylab")+rremove("xlab"),
                    Apr.d.hist+rremove("ylab")+rremove("xlab"),
                    Mayspeed.hist+rremove("ylab")+rremove("xlab"),
                    May.ODBA.hist+rremove("ylab")+rremove("xlab"),
                    May.d.hist+rremove("ylab")+rremove("xlab"),
                    Junspeed.hist,
                    Jun.ODBA.hist,
                    Jun.d.hist,
                    ncol = 3, nrow = 8, heights = c(1,1,1,1,1,1,1,1,
                                                    1,1,1,1,1,1,1,1,
                                                    1,1,1,1,1,1,1,1), align = "h")

annotate_figure(all.hist, left = text_grob("Count (Total Minutes)", color = "black", size = 14, rot = 90))


(Novspeed.hist|Nov.ODBA.hist|Nov.d.hist)/ 
  (Decspeed.hist|Dec.ODBA.hist|Dec.d.hist)/
  (Janspeed.hist|Jan.ODBA.hist|Jan.d.hist)/
  (Febspeed.hist|Feb.ODBA.hist|Feb.d.hist)/
  (Marspeed.hist|Mar.ODBA.hist|Mar.d.hist)/
  (Aprspeed.hist|Apr.ODBA.hist|Apr.d.hist)/
  (Mayspeed.hist|May.ODBA.hist|May.d.hist)/
  (Junspeed.hist|Jun.ODBA.hist|Jun.d.hist)

####Predicting MR of Fish####
#heart rate calc

#predicting MR based on Reeve et al. 2024
hr_trends_2021$MR<-10^(9.4369-3.4368*log10(hr_trends_2021$length)+0.9120*log10(hr_trends_2021$calchr_m))
hr_trends_2022$MR<-10^(9.4369-3.4368*log10(hr_trends_2022$length)+0.9120*log10(hr_trends_2022$calchr_m))

#plotting
MR_hr.time.2021<-ggplot(hr_trends_2021, aes(x=hour,y=MR)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.6, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,400,50), limits = c(0,400))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2020/2021 - HRT")+
  theme(plot.title = element_text(size = 13),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"), 
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = expression(italic("M")*O["2"]~(mgO[2]~"\U00B7"~kg^{"-1"}~"\U00B7"~hr^{"-1"}))) +
  scale_colour_manual(labels = c("Fish 1","Fish 2", "Fish 3", "Fish 4", "Fish 5", "Fish 6", "Fish 7", "Fish 8", "Fish 9",
                                 "Fish 10", "Fish 11", "Fish 12", "Fish 13", "Fish 14"), 
                      values = c("violetred", "#800000FF","#725663FF", "#D49464FF","#FFB547FF", "lightblue","#0072B2", "#5B8FA8FF", "cyan4","darkgreen", 
                                 "chartreuse3", "orangered2", "darkslategrey","#660099"))

MR_hr.time.2022<-ggplot(hr_trends_2022, aes(x=hour,y=MR)) + 
  geom_point(colour = "lightgrey",size=0.8)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.6, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,400,50), limits = c(0,400))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2022-09-20','2023-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2022/2023 - HRT")+
  theme(plot.title = element_text(size = 13),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = expression(italic("M")*O["2"]~(mgO[2]~"\U00B7"~kg^{"-1"}~"\U00B7"~hr^{"-1"}))) +
  scale_colour_manual(labels = c("Fish 15", "Fish 16", "Fish 17", "Fish 18", "Fish 19", "Fish 20", "Fish 21"), 
                      values = c("#CC79A7","burlywood", "#000000","#009E73","#D55E00", "#E69F00","#56B4E9"))

MR_hr.time.2021|MR_hr.time.2022

#techno calculations - using dat filled to estimate temp for those that had broken tags (note that this is also averaged to the hour)

glimpse(dat_filled)

dat_filled_trends<-dat_filled %>% 
  mutate(id = (if_else(id == "Tag_1", "1",
                       if_else(id == "Tag_2", "2",
                               if_else(id == "Tag_3", "3",
                                       if_else(id == "Tag_4", "4",
                                               if_else(id == "Tag_5", "5",
                                                       if_else(id == "Tag_10", "6",
                                                               if_else(id == "Tag_11", "7",
                                                                       if_else(id == "Tag_12", "8",
                                                                               if_else(id == "Tag_13", "9",
                                                                                       if_else(id == "Tag_14", "10",
                                                                                               if_else(id == "Tag_15", "11",
                                                                                                       if_else(id == "Tag_16", "12",
                                                                                                               if_else(id == "Tag_18", "13",
                                                                                                                       if_else(id == "Tag_19", "14",
                                                                                                                               if_else(id == "Tag_20", "15",
                                                                                                                                       if_else(id == "Tag_1_2022", "16",
                                                                                                                                               if_else(id == "Tag_2_2022", "17",
                                                                                                                                                       if_else(id == "Tag_3_2022", "18",
                                                                                                                                                               if_else(id == "Tag_4_2022", "19",
                                                                                                                                                                       if_else(id == "Tag_5_2022", "20",
                                                                                                                                                                               if_else(id == "Tag_12_2022", "21",
                                                                                                                                                                                       if_else(id == "Tag_13_2022", "22",
                                                                                                                                                                                               if_else(id == "Tag_14_2022", "23", 
                                                                                                                                                                                                       if_else(id == "Tag_11_2022", "24",
                                                                                                                                                                                                               if_else(id == "Tag_201033_2022", "25",
                                                                                                                                                                                                                       NA_character_)))))))))))))))))))))))))))

dat_filled_trends$id<-as.integer(dat_filled_trends$id)
dat_filled_trends$weight<-as.numeric(dat_filled_trends$weight)

dat_trends_2021<-dat_filled_trends %>% 
  filter(timestamp < "2021-07-20 00:00:00")
dat_trends_2021<-dat_trends_2021[!is.na(dat_trends_2021$ODBA_ms2),]

dat_trends_2022<-dat_filled_trends %>% 
  filter(timestamp > "2021-07-20 00:00:00" & timestamp < "2022-09-01 00:00:00")
dat_trends_2022<-dat_trends_2022[!is.na(dat_trends_2022$ODBA_ms2),]

dat_trends_2021$weight<-as.numeric(dat_trends_2021$weight)

dat_trends_2021$MR<-10^(4.4489+0.0426*(dat_trends_2021$temp_m)-1.1446*log10(dat_trends_2021$length)+0.3321*log10(dat_trends_2021$ODBA_ms2))
dat_trends_2021<-dat_trends_2021 %>% 
  mutate(MR = if_else(ODBA_ms2 < 0.1, 10^(4.4489+0.0426*(temp_m)-1.1446*log10(length)+0.3321*log10(0.1)), MR))

dat_trends_2022$weight<-as.numeric(dat_trends_2022$weight)

dat_trends_2022$MR<-10^(4.4489+0.0426*(dat_trends_2022$temp_m)-1.1446*log10(dat_trends_2022$length)+0.3321*log10(dat_trends_2022$ODBA_ms2))
dat_trends_2022<-dat_trends_2022 %>% 
  mutate(MR = if_else(ODBA_ms2 < 0.1, 10^(4.4489+0.0426*(temp_m)-1.1446*log10(length)+0.3321*log10(0.1)), MR))

#plotting

cols = c("1" = "violetred", "2" = "#800000FF", "3" = "#725663FF", "4" = "#D49464FF","5" = "#FFB547FF", "6" ="lightblue","7" ="#0072B2", "8" ="#5B8FA8FF","9" = "cyan4","10" ="darkgreen", 
         "11" = "chartreuse3", "12" = "orangered2", 
         "13" ="#CC79A7","14" ="burlywood","15" = "#000000","16" ="#009E73","17" ="#D55E00", "18" ="#660099","19" ="#E69F00",
         "20" ="#56B4E9", "21" ="#339900", "22" ="steelblue", "23" ="#FF6666", "24" ="#33ff99","25" ="#FF3333")

MR.time.2021<-ggplot(dat_trends_2021, aes(x=timestamp,y=MR)) + 
  geom_point(colour = "lightgrey",size=0.6)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.6, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,400,50), limits = c(0,400))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2020-09-20','2021-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2020/2021 - Axy-5")+
  theme(plot.title = element_text(size = 13),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"), 
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank(), 
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = expression(italic("M")*O["2"]~(mgO[2]~"\U00B7"~kg^{"-1"}~"\U00B7"~hr^{"-1"}))) +
  scale_colour_manual(values = cols)

MR.time.2022<-ggplot(dat_trends_2022, aes(x=timestamp,y=MR)) + 
  geom_point(colour = "lightgrey",size=0.6)+
  geom_smooth(aes(colour = factor(id)), method = "gam", se = FALSE, size = 0.6, alpha = 0.8)+
  scale_y_continuous(breaks = seq(0,400,50), limits = c(0,400))+ 
  scale_x_datetime(date_breaks = "1 month", date_minor_breaks = "1 week",
                   date_labels = "%b",
                   limits=as.POSIXct(c('2021-09-20','2022-07-01'), format = "%Y-%m-%d"))+
  theme_classic()+
  ggtitle("Winter 2021/2022 - Axy-5")+
  theme(plot.title = element_text(size = 13),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(colour = "Fish", x = "Date", y = expression(italic("M")*O["2"]~(mgO[2]~"\U00B7"~kg^{"-1"}~"\U00B7"~hr^{"-1"}))) +
  scale_colour_manual(values = cols)


(MR.time.2021|MR.time.2022)/
  (MR_hr.time.2021|MR_hr.time.2022)

#calculating MR Q10

glimpse(dat_trends_2022)

dat_trends_2022$id<-as.character(dat_trends_2022$id)
dat_trends_2021$id<-as.character(dat_trends_2021$id)

glimpse(hr_trends_2022)

MR_combined<-rbind((hr_trends_2021 %>% dplyr::select(id = tag, timestamp = hour, temp_m, MR)),
                   (hr_trends_2022 %>% dplyr::select(id = tag, timestamp = hour, temp_m, MR)),
                   (dat_trends_2021 %>% dplyr::select(id, timestamp, temp_m, MR)),
                   (dat_trends_2022 %>% dplyr::select(id, timestamp, temp_m, MR)))

q20c<-MR_combined %>% 
  filter(temp_m > 19 & temp_m < 21) %>% 
  group_by(id) %>% 
  summarise(temp_20 = mean(temp_m, na.rm = TRUE),
            MR_20 = mean(MR, na.rm = TRUE))

q15c<-MR_combined %>% 
  filter(temp_m > 14 & temp_m < 16) %>% 
  group_by(id) %>% 
  summarise(temp_15 = mean(temp_m, na.rm = TRUE),
            MR_15 = mean(MR, na.rm = TRUE)) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .))

q10c<-MR_combined %>% 
  filter(temp_m > 9 & temp_m < 11) %>% 
  group_by(id) %>% 
  summarise(temp_10 = mean(temp_m, na.rm = TRUE),
            MR_10 = mean(MR, na.rm = TRUE)) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 


q5c<-MR_combined %>% 
  filter(temp_m > 4 & temp_m < 6) %>% 
  group_by(id) %>% 
  summarise(temp_5 = mean(temp_m, na.rm = TRUE),
            MR_5 = mean(MR, na.rm = TRUE)) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

Q10_df<-q5c %>% 
  left_join(q10c, by = "id") %>% 
  left_join(q15c, by = "id") %>% 
  left_join(q20c, by = "id") %>% 
  mutate(q10v5 = ((MR_5/MR_10)^(10/(temp_5 - temp_10)))) %>% 
  mutate(q15v5 = ((MR_5/MR_15)^(10/(temp_5 - temp_15)))) %>% 
  mutate(q20v5 = ((MR_5/MR_20)^(10/(temp_5 - temp_20)))) %>% 
  mutate(q20v10 = ((MR_10/MR_20)^(10/(temp_10 - temp_20))))

####Time Inactive/Bursting####

#Inactive, bursting, and critical swimming behaviour based on Reeve et al. 2024
glimpse(dat)

dat<-dat %>% 
  mutate(pred_Ucrit = (10^(2.0351+0.0276*temp_m-0.9138*log10(length))))

swim_behaviour_month_id<-dat %>% 
  group_by(id, month_character) %>% 
  summarise(speed_m = mean(speed_bls, na.rm = TRUE),
            speed_sd = sd(speed_bls, na.rm = TRUE),
            temp_avg = mean(temp_m, na.rm = TRUE),
            count_speed = sum(!(is.na(speed_bls))),
            count_burst = sum(if_else(ODBA_ms2 >= 1.5, 1, 0), na.rm = TRUE),#changed from pred_Ucrit to ODBA_m2 > 1.5
            count_crit = sum(if_else(speed_bls>=pred_Ucrit, 1, 0), na.rm = TRUE),
            count_inactive = sum(if_else(ODBA_ms2 < 0.1, 1, 0), na.rm = TRUE)) %>% #try different speeds here - right now assuming 0.25 BL/s as essentially inactive (equates to ~10cm/s or less)
  mutate(prop_burst = (count_burst/count_speed)*100) %>% 
  mutate(prop_inactive = (count_inactive/count_speed)*100) %>% 
  mutate(prop_crit = (count_crit/count_speed)*100) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  drop_na(speed_sd) %>% 
  drop_na(temp_avg)

swim_behaviour_month<-swim_behaviour_month_id %>% 
  group_by(month_character) %>% 
  summarise(speed_mean = mean(speed_m, na.rm = TRUE),
            speed_sd = sd(speed_m, na.rm = TRUE),
            temp_mean = mean(temp_avg),
            prop_burst_m = mean(prop_burst),
            prop_burst_sd = sd(prop_burst),
            prop_inactive_m = mean(prop_inactive),
            prop_inactive_sd = sd(prop_inactive),
            prop_crit_m = mean(prop_inactive),
            prop_crit_sd = sd(prop_inactive),
            sum_total = sum(count_speed)) 

#fixing month values to present properly on figure
swim_behaviour_month_id<-swim_behaviour_month_id %>% 
  mutate(month_character = if_else(month_character == "K_Sep", "1_Sep",
                                   if_else(month_character == "L_Oct", "2_Oct", month_character))) %>% 
  filter(!(month_character %in% c("I_Jul", "J_Aug"))) %>% 
  drop_na(month_character)

prop_inactive_month<-
  ggplot(data = test_speed_month_id,aes(x=month_character,y=prop_inactive)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0, 5)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 13, angle = 90),
        #axis.text.x = element_blank(),
        legend.position = "none")+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)+
  labs(x = "Month", y = "Proportion of Time\n Inactive (%)", size = 6)

prop_bursting_month<-
  ggplot(data = test_speed_month_id,aes(x=month_character,y=prop_burst)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 13, angle = 90),
        legend.position = "none")+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)+
  labs(x = "Month", y = "Proportion of Time\n Burst Swimming (%)", size = 6)

prop_crit_month<- #include statistics but do not include figure.... 
  ggplot(data = test_speed_month_id,aes(x=month_character,y=prop_crit)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 13, angle = 90),
        legend.position = "none")+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)+
  labs(x = "Month", y = "Proportion of Time\n Critical Swimming (%)", size = 6)

prop_inactive_month/prop_bursting_month

####Calculating Averages by Month####

glimpse(dat_filled_MR)

dat_filled_MR$date<-as.numeric(strftime(dat_filled_MR$timestamp, format = "%j"))

dat_filled_MR$DO <- predict(DOForest,dat_filled_MR)

tech_averages_month_id<-dat_filled_MR %>% 
  group_by(id, month_character) %>% 
  summarise(speed_avg = mean(speed_bls, na.rm = TRUE),
            speed_sd = sd(speed_bls, na.rm = TRUE),
            temp_avg = mean(temp_m, na.rm = TRUE),
            temp_sd = sd(temp_m, na.rm = TRUE),
            DO_avg = mean(DO, na.rm = TRUE),
            DO_sd = sd(DO,na.rm = TRUE),
            depth_avg = mean(depth_m, na.rm = TRUE),
            depth_sd = sd(depth_m, na.rm = TRUE),
            ODBA_avg = mean(ODBA_ms2, na.rm = TRUE),
            ODBA_sd = sd(ODBA_ms2, na.rm = TRUE),
            MR_avg = mean(MR, na.rm = TRUE),
            MR_sd = sd(MR, na.rm = TRUE),
            count_speed = sum(if_else(speed_bls > -2, 1, 0), na.rm = TRUE),
            count_burst = sum(if_else(ODBA_ms2 >= 1.5, 1, 0), na.rm = TRUE), #updated ODBA > 1.5
            count_inactive = sum(if_else(speed_bls <= 0.25, 1, 0), na.rm = TRUE)) %>% #updated speed_bls < 0.25
  mutate(prop_burst = (count_burst/count_speed)*100) %>% 
  mutate(prop_inactive = (count_inactive/count_speed)*100)%>% 
  mutate(month_character = if_else(month_character == "K_Sep", "1_Sep",
                                   if_else(month_character == "L_Oct", "2_Oct", month_character))) %>% 
  filter(!(month_character %in% c("I_Jul", "J_Aug"))) %>% 
  drop_na(month_character)

tech_averages_month_id$bar<-(tech_averages_month_id$depth_avg*997.05*9.80665)/100000 + 1.013

tech_averages_month_id$C<-exp(7.7117 - 1.31403*log(tech_averages_month_id$temp_avg + 45.93))

tech_averages_month_id$Pw<-11.8571-(3840.7/tech_averages_month_id$temp_avg)-(216961/(tech_averages_month_id$temp_avg^2))

tech_averages_month_id$O<-0.000975-((1.426*10^(-5))*tech_averages_month_id$temp_avg)+((6.436*10^(-8))*tech_averages_month_id$temp_avg)

tech_averages_month_id$Cp<-(tech_averages_month_id$C*tech_averages_month_id$bar)*
  (((1-tech_averages_month_id$Pw/tech_averages_month_id$bar)*(1-tech_averages_month_id$O*tech_averages_month_id$bar))/
     ((1-tech_averages_month_id$Pw)*(1-tech_averages_month_id$O)))

tech_averages_month_id$DO_mgL<-(tech_averages_month_id$DO_avg*tech_averages_month_id$Cp)/100

tech_averages_month_id<-tech_averages_month_id %>% 
  left_join(weight_df, by = "id")

tech_averages_month<-tech_averages_month_id %>% 
  group_by(month_character) %>% 
  summarise(speed_m = mean(speed_avg, na.rm = TRUE),
            speed_sd = sd(speed_avg, na.rm = TRUE),
            temp_m = mean(temp_avg, na.rm = TRUE),
            temp_sd = sd(temp_avg, na.rm = TRUE),
            DO_m = mean(DO_avg, na.rm = TRUE),
            DO_sd = sd(DO_avg,na.rm = TRUE),
            DO_mgL_m = mean(DO_mgL, na.rm = TRUE),
            DO_mgL_sd = sd(DO_mgL, na.rm = TRUE),
            depth_m = mean(depth_avg, na.rm = TRUE),
            depth_sd = sd(depth_avg, na.rm = TRUE),
            ODBA_m = mean(ODBA_avg, na.rm = TRUE),
            ODBA_sd = sd(ODBA_avg, na.rm = TRUE),
            MR_m = mean(MR_avg, na.rm = TRUE),
            MR_sd = sd(MR_avg, na.rm = TRUE),
            prop_burst_m = mean(prop_burst),
            prop_burst_sd = sd(prop_burst),
            prop_inactive_m = mean(prop_inactive),
            prop_inactive_sd = sd(prop_inactive),
            sum_total = sum(count_speed)) 

table(tech_averages_month_id$month_character)

glimpse(dat)

#heart rate data

glimpse(datHR)

datHR$month <- strftime(datHR$timestamp, format="%m") 

datHR<-datHR%>% 
  mutate(month_character=case_when(month=="01"~'C_Jan',
                                   month=="02"~'D_Feb',
                                   month=="03"~'E_Mar',
                                   month=="04"~'F_Apr',
                                   month=="05"~'G_May',
                                   month=="06"~'H_Jun',
                                   month=="07"~'I_Jul',
                                   month=="08"~'J_Aug',
                                   month=="09"~'1_Sep',
                                   month=="10"~'2_Oct',
                                   month=="11"~'A_Nov',
                                   month=="12"~'B_Dec',
                                   TRUE~'Other')) 

table(datHR$month_character)

datHR<-datHR %>% 
  left_join(hr_weight, by = "tag")

datHR$weight<-as.numeric(datHR$weight)
datHR$length<-as.numeric(datHR$length)
datHR$calchr<-as.numeric(datHR$calchr)

glimpse(datHR)

min(datHR$calchr)

datHR$MR<-10^(9.4369-3.4368*log10(datHR$length)+0.9120*log10(datHR$calchr))

max(datHR$MR, na.rm = TRUE)

dat_month_data<-datHR %>% 
  dplyr::select(tag, month_character, hour = timestamp) %>% 
  distinct()

datHR_hour<-datHR %>% 
  group_by(tag,hour) %>% 
  summarise(temp_m = mean(temp, na.rm =TRUE),
            calchr_m = mean(calchr, na.rm = TRUE),
            MR_m = mean(MR, na.rm =TRUE),
            weight = mean(weight, na.rm = TRUE)) %>% 
  left_join(dat_month_data, by = c("tag", "hour"))

glimpse(datHR_hour)

star_averages_month_id<-datHR %>% 
  group_by(tag, month_character) %>% 
  summarise(temp_avg = mean(temp, na.rm = TRUE),
            temp_sd = sd(temp, na.rm = TRUE),
            calchr_m = mean(calchr, na.rm = TRUE),
            calchr_sd = sd(calchr, na.rm = TRUE),
            MR_avg = mean(MR, na.rm = TRUE),
            MR_sd = sd(MR, na.rm = TRUE)) %>% 
  filter(!(month_character %in% c("I_Jul", "J_Aug")))

table(star_averages_month_id$tag, star_averages_month_id$month_character)

star_averages_month<-star_averages_month_id %>% 
  group_by(month_character) %>% 
  summarise(temp_mean = mean(temp_avg, na.rm = TRUE),
            temp_sd = sd(temp_avg, na.rm = TRUE),
            calchr_mean = mean(calchr_m, na.rm = TRUE),
            calchr_sd = sd(calchr_m, na.rm = TRUE),
            MR_mean = mean(MR_avg, na.rm = TRUE),
            MR_sd = sd(MR_avg, na.rm = TRUE)) 


speed_month<-
  ggplot(data = tech_averages_month_id,aes(x=month_character,y=speed_avg)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(breaks = seq(0.25, 1.5, 0.25), limits = c(0.25, 1.5)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        #axis.text.x = element_text(size = 13, angle = 90),
        axis.text.x = element_blank(),
        legend.position = "none")+
  labs(x = "Month", y = expression("Speed"~(BL~"\U00B7"~sec^{"-1"})), size = 6)+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)

ODBA_month<-
  ggplot(data = tech_averages_month_id,aes(x=month_character,y=ODBA_avg)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  #scale_y_continuous(breaks = seq(0, 0.3, 0.05), limits = c(0, 0.3)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  labs(x = "Month", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})), size = 6)+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)

depth_month<-
  ggplot(data = tech_averages_month_id,aes(x=month_character,y=(0-depth_avg))) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(limits = c(-5, 0),
                     breaks = seq(-5,0,1),
                     labels = c("5", "4", "3", "2", "1", "0"))+
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  labs(x = "Month", y = "Depth (m)", size = 6)+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)

tech_MR_month<-
  ggplot(data = tech_averages_month_id,aes(x=month_character,y=MR_avg)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(breaks = seq(0, 250, 50), limits = c(0, 250)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  labs(x = "Month", y = expression("Acc"~italic("M")*O["2"]~(mgO[2]~"\U00B7"~kg^{"-1"}~"\U00B7"~hr^{"-1"})), size = 6)+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)

star_MR_month<-
  ggplot(data = star_averages_month_id,aes(x=month_character,y=MR_avg)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(breaks = seq(25, 100, 25), limits = c(10, 100)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  labs(x = "Month", y = expression("HR"~italic("M")*O["2"]~(mgO[2]~"\U00B7"~kg^{"-1"}~"\U00B7"~hr^{"-1"})), size = 6)+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)

HR_month<-
  ggplot(data = star_averages_month_id,aes(x=month_character,y=calchr_m)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  labs(x = "Month", y = "Heart Rate (bpm)", size = 6)+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)

DO_month<-
  ggplot(data = tech_averages_month_id,aes(x=month_character,y=DO_mgL)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 12)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  labs(x = "Month", y = expression("Dissolved Oxygen"~"("*mg~L^{"-1"}*")"), size = 6)+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)

min(tech_averages_month_id$temp_avg)
max(tech_averages_month_id$temp_avg)

Temp_month<-
  ggplot(data = tech_averages_month_id,aes(x=month_character,y=temp_avg)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width=0.3,height=0,size=2) +
  scale_y_continuous(breaks = seq(4, 24, 4), limits = c(3, 24)) + 
  scale_x_discrete(labels = c("Sep","Oct","Nov", "Dec","Jan","Feb","Mar","Apr","May","June"))+ 
  scale_colour_manual(values = c("black")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  labs(x = "Month", y = "Temperature (°C)", size = 6)+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 4, alpha = 0.8)

(Temp_month/tech_MR_month/ODBA_month/speed_month/prop_bursting_month) | (DO_month/star_MR_month/HR_month/depth_month/prop_inactive_month)

####Calculating Averages by Month/Year####

glimpse(dat_filled_MR)

tech_averages_month_id_W1<-dat_filled_MR %>%
  mutate(pred_Ucrit = (10^(2.0351+0.0276*temp_m-0.9138*log10(length)))) %>% 
  mutate(year = if_else(id %in% c("Tag_1", "Tag_10", "Tag_11", "Tag_12", "Tag_16", "Tag_13", 
                                  "Tag_14", "Tag_15", "Tag_16", "Tag_2", "Tag_3", "Tag_4", "Tag_5"), "year_1", "year_2")) %>% 
  filter(year == "year_1") %>% 
  group_by(id, month_character) %>% 
  summarise(speed_avg = mean(speed_bls, na.rm = TRUE),
            speed_sd = sd(speed_bls, na.rm = TRUE),
            temp_avg = mean(temp_m, na.rm = TRUE),
            temp_sd = sd(temp_m, na.rm = TRUE),
            DO_avg = mean(DO, na.rm = TRUE),
            DO_sd = sd(DO,na.rm = TRUE),
            depth_avg = mean(depth_m, na.rm = TRUE),
            depth_sd = sd(depth_m, na.rm = TRUE),
            ODBA_avg = mean(ODBA_ms2, na.rm = TRUE),
            ODBA_sd = sd(ODBA_ms2, na.rm = TRUE),
            MR_avg = mean(MR, na.rm = TRUE),
            MR_sd = sd(MR, na.rm = TRUE),
            count_speed = sum(!(is.na(speed_bls))),
            count_burst = sum(if_else(ODBA_ms2 >= 1.5, 1, 0), na.rm = TRUE),#changed from pred_Ucrit to ODBA_m2 > 1.5
            count_crit = sum(if_else(speed_bls>=pred_Ucrit, 1, 0), na.rm = TRUE),
            count_inactive = sum(if_else(ODBA_ms2 < 0.1, 1, 0), na.rm = TRUE)) %>% #try different speeds here - right now assuming 0.25 BL/s as essentially inactive (equates to ~10cm/s or less)
  mutate(prop_burst = (count_burst/count_speed)*100) %>% 
  mutate(prop_inactive = (count_inactive/count_speed)*100) %>% 
  mutate(prop_crit = (count_crit/count_speed)*100) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  drop_na(speed_sd) %>% 
  mutate(month_character = if_else(month_character == "K_Sep", "1_Sep",
                                   if_else(month_character == "L_Oct", "2_Oct", month_character))) %>% 
  filter(!(month_character %in% c("I_Jul", "J_Aug"))) %>% 
  drop_na(month_character)

tech_averages_month_id_W1$bar<-(tech_averages_month_id_W1$depth_avg*997.05*9.80665)/100000 + 1.013

max(tech_averages_month_id_W1$bar,na.rm = TRUE)

max(tech_averages_month_id_W1$depth_avg, na.rm = TRUE)

tech_averages_month_id_W1$bar<-(tech_averages_month_id_W1$depth_avg*997.05*9.80665)/100000 + 1.013

tech_averages_month_id_W1$C<-exp(7.7117 - 1.31403*log(tech_averages_month_id_W1$temp_avg + 45.93))

tech_averages_month_id_W1$Pw<-11.8571-(3840.7/tech_averages_month_id_W1$temp_avg)-(216961/(tech_averages_month_id_W1$temp_avg^2))

tech_averages_month_id_W1$O<-0.000975-((1.426*10^(-5))*tech_averages_month_id_W1$temp_avg)+((6.436*10^(-8))*tech_averages_month_id_W1$temp_avg)

tech_averages_month_id_W1$Cp<-(tech_averages_month_id_W1$C*tech_averages_month_id_W1$bar)*
  (((1-tech_averages_month_id_W1$Pw/tech_averages_month_id_W1$bar)*(1-tech_averages_month_id_W1$O*tech_averages_month_id_W1$bar))/
     ((1-tech_averages_month_id_W1$Pw)*(1-tech_averages_month_id_W1$O)))

tech_averages_month_id_W1$DO_mgL<-(tech_averages_month_id_W1$DO_avg*tech_averages_month_id_W1$Cp)/100

tech_averages_month_W1<-tech_averages_month_id_W1 %>% 
  group_by(month_character) %>% 
  summarise(speed_m = mean(speed_avg, na.rm = TRUE),
            speed_sd = sd(speed_avg, na.rm = TRUE),
            temp_m = mean(temp_avg, na.rm = TRUE),
            temp_sd = sd(temp_avg, na.rm = TRUE),
            DO_m = mean(DO_avg, na.rm = TRUE),
            DO_sd = sd(DO_avg,na.rm = TRUE),
            DO_mgL_m = mean(DO_mgL, na.rm = TRUE),
            DO_mgL_sd = sd(DO_mgL, na.rm = TRUE),
            depth_m = mean(depth_avg, na.rm = TRUE),
            depth_sd = sd(depth_avg, na.rm = TRUE),
            ODBA_m = mean(ODBA_avg, na.rm = TRUE),
            ODBA_sd = sd(ODBA_avg, na.rm = TRUE),
            MR_m = mean(MR_avg, na.rm = TRUE),
            MR_sd = sd(MR_avg, na.rm = TRUE),
            prop_burst_m = mean(prop_burst, na.rm = TRUE),
            prop_burst_sd = sd(prop_burst, na.rm = TRUE),
            prop_crit_m = mean(prop_crit, na.rm = TRUE),
            prop_crit_sd = sd(prop_crit, na.rm = TRUE),
            prop_inactive_m = mean(prop_inactive, na.rm = TRUE),
            prop_inactive_sd = sd(prop_inactive, na.rm = TRUE),
            sum_total = sum(count_speed, na.rm = TRUE)) 

table(tech_averages_month_id_W1$month_character)

#year 2

tech_averages_month_id_W2<-dat_filled_MR %>%
  mutate(pred_Ucrit = (10^(2.0351+0.0276*temp_m-0.9138*log10(length)))) %>% 
  mutate(year = if_else(id %in% c("Tag_1", "Tag_10", "Tag_11", "Tag_12", "Tag_16", "Tag_13", 
                                  "Tag_14", "Tag_15", "Tag_16", "Tag_2", "Tag_3", "Tag_4", "Tag_5"), "year_1", "year_2")) %>% 
  filter(year == "year_2") %>% 
  group_by(id, month_character) %>% 
  summarise(speed_avg = mean(speed_bls, na.rm = TRUE),
            speed_sd = sd(speed_bls, na.rm = TRUE),
            temp_avg = mean(temp_m, na.rm = TRUE),
            temp_sd = sd(temp_m, na.rm = TRUE),
            DO_avg = mean(DO, na.rm = TRUE),
            DO_sd = sd(DO,na.rm = TRUE),
            depth_avg = mean(depth_m, na.rm = TRUE),
            depth_sd = sd(depth_m, na.rm = TRUE),
            ODBA_avg = mean(ODBA_ms2, na.rm = TRUE),
            ODBA_sd = sd(ODBA_ms2, na.rm = TRUE),
            MR_avg = mean(MR, na.rm = TRUE),
            MR_sd = sd(MR, na.rm = TRUE),
            count_speed = sum(!(is.na(speed_bls))),
            count_burst = sum(if_else(ODBA_ms2 >= 1.5, 1, 0), na.rm = TRUE),#changed from pred_Ucrit to ODBA_m2 > 1.5
            count_crit = sum(if_else(speed_bls>=pred_Ucrit, 1, 0), na.rm = TRUE),
            count_inactive = sum(if_else(ODBA_ms2 < 0.1, 1, 0), na.rm = TRUE)) %>% #try different speeds here - right now assuming 0.25 BL/s as essentially inactive (equates to ~10cm/s or less)
  mutate(prop_burst = (count_burst/count_speed)*100) %>% 
  mutate(prop_inactive = (count_inactive/count_speed)*100) %>% 
  mutate(prop_crit = (count_crit/count_speed)*100) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  drop_na(speed_sd) %>% 
  mutate(month_character = if_else(month_character == "K_Sep", "1_Sep",
                                   if_else(month_character == "L_Oct", "2_Oct", month_character))) %>% 
  filter(!(month_character %in% c("I_Jul", "J_Aug"))) %>% 
  drop_na(month_character)

tech_averages_month_id_W2$bar<-(tech_averages_month_id_W2$depth_avg*997.05*9.80665)/100000 + 1.013

max(tech_averages_month_id_W2$bar,na.rm = TRUE)

max(tech_averages_month_id_W2$depth_avg, na.rm = TRUE)

tech_averages_month_id_W2$bar<-(tech_averages_month_id_W2$depth_avg*997.05*9.80665)/100000 + 1.013

tech_averages_month_id_W2$C<-exp(7.7117 - 1.31403*log(tech_averages_month_id_W2$temp_avg + 45.93))

tech_averages_month_id_W2$Pw<-11.8571-(3840.7/tech_averages_month_id_W2$temp_avg)-(216961/(tech_averages_month_id_W2$temp_avg^2))

tech_averages_month_id_W2$O<-0.000975-((1.426*10^(-5))*tech_averages_month_id_W2$temp_avg)+((6.436*10^(-8))*tech_averages_month_id_W2$temp_avg)

tech_averages_month_id_W2$Cp<-(tech_averages_month_id_W2$C*tech_averages_month_id_W2$bar)*
  (((1-tech_averages_month_id_W2$Pw/tech_averages_month_id_W2$bar)*(1-tech_averages_month_id_W2$O*tech_averages_month_id_W2$bar))/
     ((1-tech_averages_month_id_W2$Pw)*(1-tech_averages_month_id_W2$O)))

tech_averages_month_id_W2$DO_mgL<-(tech_averages_month_id_W2$DO_avg*tech_averages_month_id_W2$Cp)/100

tech_averages_month_W2<-tech_averages_month_id_W2 %>% 
  group_by(month_character) %>% 
  summarise(speed_m = mean(speed_avg, na.rm = TRUE),
            speed_sd = sd(speed_avg, na.rm = TRUE),
            temp_m = mean(temp_avg, na.rm = TRUE),
            temp_sd = sd(temp_avg, na.rm = TRUE),
            DO_m = mean(DO_avg, na.rm = TRUE),
            DO_sd = sd(DO_avg,na.rm = TRUE),
            DO_mgL_m = mean(DO_mgL, na.rm = TRUE),
            DO_mgL_sd = sd(DO_mgL, na.rm = TRUE),
            depth_m = mean(depth_avg, na.rm = TRUE),
            depth_sd = sd(depth_avg, na.rm = TRUE),
            ODBA_m = mean(ODBA_avg, na.rm = TRUE),
            ODBA_sd = sd(ODBA_avg, na.rm = TRUE),
            MR_m = mean(MR_avg, na.rm = TRUE),
            MR_sd = sd(MR_avg, na.rm = TRUE),
            prop_burst_m = mean(prop_burst, na.rm = TRUE),
            prop_burst_sd = sd(prop_burst, na.rm = TRUE),
            prop_crit_m = mean(prop_crit, na.rm = TRUE),
            prop_crit_sd = sd(prop_crit, na.rm = TRUE),
            prop_inactive_m = mean(prop_inactive, na.rm = TRUE),
            prop_inactive_sd = sd(prop_inactive, na.rm = TRUE),
            sum_total = sum(count_speed, na.rm = TRUE)) 

table(tech_averages_month_id_W2$month_character)

#inactivity/bursting 

speed_month_id_W1<-dat %>% 
  mutate(pred_Ucrit = (10^(2.0351+0.0276*temp_m-0.9138*log10(length)))) %>% 
  mutate(year = if_else(id %in% c("Tag_1", "Tag_10", "Tag_11", "Tag_12", "Tag_16", "Tag_13", 
                                  "Tag_14", "Tag_15", "Tag_16", "Tag_2", "Tag_3", "Tag_4", "Tag_5"), "year_1", "year_2")) %>% 
  filter(year == "year_1") %>%
  group_by(id, month_character) %>% 
  summarise(speed_m = mean(speed_bls, na.rm = TRUE),
            speed_sd = sd(speed_bls, na.rm = TRUE),
            temp_avg = mean(temp_m, na.rm = TRUE),
            count_speed = sum(!(is.na(speed_bls))),
            count_burst = sum(if_else(ODBA_ms2 >= 1.5, 1, 0), na.rm = TRUE),#changed from pred_Ucrit to ODBA_m2 > 1.5
            count_crit = sum(if_else(speed_bls>=pred_Ucrit, 1, 0), na.rm = TRUE),
            count_inactive = sum(if_else(ODBA_ms2 < 0.1, 1, 0), na.rm = TRUE)) %>% #try different speeds here - right now assuming 0.25 BL/s as essentially inactive (equates to ~10cm/s or less)
  mutate(prop_burst = (count_burst/count_speed)*100) %>% 
  mutate(prop_inactive = (count_inactive/count_speed)*100) %>% 
  mutate(prop_crit = (count_crit/count_speed)*100) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  drop_na(speed_sd) %>% 
  drop_na(temp_avg)

speed_month_W1<-test_speed_month_id_W1 %>% 
  group_by(month_character) %>% 
  summarise(speed_mean = mean(speed_m, na.rm = TRUE),
            speed_sd = sd(speed_m, na.rm = TRUE),
            temp_mean = mean(temp_avg),
            prop_burst_m = mean(prop_burst),
            prop_burst_sd = sd(prop_burst),
            prop_crit_m = mean(prop_crit),
            prop_crit_sd = sd(prop_crit),
            prop_inactive_m = mean(prop_inactive),
            prop_inactive_sd = sd(prop_inactive),
            sum_total = sum(count_speed)) 

speed_month_id_W2<-dat %>% 
  mutate(pred_Ucrit = (10^(2.0351+0.0276*temp_m-0.9138*log10(length)))) %>% 
  mutate(year = if_else(id %in% c("Tag_1", "Tag_10", "Tag_11", "Tag_12", "Tag_16", "Tag_13", 
                                  "Tag_14", "Tag_15", "Tag_16", "Tag_2", "Tag_3", "Tag_4", "Tag_5"), "year_1", "year_2")) %>% 
  filter(year == "year_2") %>%
  group_by(id, month_character) %>% 
  summarise(speed_m = mean(speed_bls, na.rm = TRUE),
            speed_sd = sd(speed_bls, na.rm = TRUE),
            temp_avg = mean(temp_m, na.rm = TRUE),
            count_speed = sum(!(is.na(speed_bls))),
            count_burst = sum(if_else(ODBA_ms2 >= 1.5, 1, 0), na.rm = TRUE),#changed from pred_Ucrit to ODBA_m2 > 1.5
            count_crit = sum(if_else(speed_bls>=pred_Ucrit, 1, 0), na.rm = TRUE),
            count_inactive = sum(if_else(ODBA_ms2 < 0.1, 1, 0), na.rm = TRUE)) %>% #try different speeds here - right now assuming 0.25 BL/s as essentially inactive (equates to ~10cm/s or less)
  mutate(prop_burst = (count_burst/count_speed)*100) %>% 
  mutate(prop_inactive = (count_inactive/count_speed)*100) %>% 
  mutate(prop_crit = (count_crit/count_speed)*100) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  drop_na(speed_sd)  %>% 
  drop_na(temp_avg)

speed_month_W2<-test_speed_month_id_W2 %>% 
  group_by(month_character) %>% 
  summarise(speed_mean = mean(speed_m, na.rm = TRUE),
            speed_sd = sd(speed_m, na.rm = TRUE),
            temp_mean = mean(temp_avg),
            prop_burst_m = mean(prop_burst),
            prop_burst_sd = sd(prop_burst),
            prop_crit_m = mean(prop_crit),
            prop_crit_sd = sd(prop_crit),
            prop_inactive_m = mean(prop_inactive),
            prop_inactive_sd = sd(prop_inactive),
            sum_total = sum(count_speed)) 

table(test_speed_month_id$id)
table(test_speed_month_id_W2$month_character)
table(test_speed_month_id_W1$month_character)

#star oddi

glimpse(datHR)

star_averages_month_id_W1<-datHR %>% 
  mutate(year = if_else(tag %in% c("HR0855", "HR0858_2", "HR0859_2", "HR0860_2", "HR0861_2", "HR1006", "HR1010"), "year_2", "year_1")) %>% 
  filter(year == "year_1") %>%
  group_by(tag, month_character) %>% 
  summarise(temp_avg = mean(temp, na.rm = TRUE),
            temp_sd = sd(temp, na.rm = TRUE),
            calchr_m = mean(calchr, na.rm = TRUE),
            calchr_sd = sd(calchr, na.rm = TRUE),
            MR_avg = mean(MR, na.rm = TRUE),
            MR_sd = sd(MR, na.rm = TRUE),
            count_HR = sum(if_else(calchr > 0, 1, 0), na.rm = TRUE),
            count_all = sum(if_else(hr > 0, 1, 0), na.rm = TRUE)) %>%
  mutate(prop_ECG = count_HR/(count_all*0.16666666666666)*100) %>% 
  filter(!(month_character %in% c("I_Jul", "J_Aug")))

star_averages_month_W1<-star_averages_month_id_W1 %>% 
  group_by(month_character) %>% 
  summarise(temp_mean = mean(temp_avg, na.rm = TRUE),
            temp_sd = sd(temp_avg, na.rm = TRUE),
            calchr_mean = mean(calchr_m, na.rm = TRUE),
            calchr_sd = sd(calchr_m, na.rm = TRUE),
            MR_mean = mean(MR_avg, na.rm = TRUE),
            MR_sd = sd(MR_avg, na.rm = TRUE),
            prop_ECG_m = mean(prop_ECG),
            prop_ECG_sd = sd(prop_ECG))


star_averages_month_id_W2<-datHR %>% 
  mutate(year = if_else(tag %in% c("HR0855", "HR0858_2", "HR0859_2", "HR0860_2", "HR0861_2", "HR1006", "HR1010"), "year_2", "year_1")) %>% 
  filter(year == "year_2") %>%
  group_by(tag, month_character) %>% 
  summarise(temp_avg = mean(temp, na.rm = TRUE),
            temp_sd = sd(temp, na.rm = TRUE),
            calchr_m = mean(calchr, na.rm = TRUE),
            calchr_sd = sd(calchr, na.rm = TRUE),
            MR_avg = mean(MR, na.rm = TRUE),
            MR_sd = sd(MR, na.rm = TRUE),
            count_HR = sum(if_else(calchr > 0, 1, 0), na.rm = TRUE),
            count_all = sum(if_else(QI > 0, 1, 0), na.rm = TRUE)) %>%
  mutate(prop_ECG = count_HR/(count_all)*100) %>% 
  filter(!(month_character %in% c("I_Jul", "J_Aug")))

star_averages_month_W2<-star_averages_month_id_W2 %>% 
  group_by(month_character) %>% 
  summarise(temp_mean = mean(temp_avg, na.rm = TRUE),
            temp_sd = sd(temp_avg, na.rm = TRUE),
            calchr_mean = mean(calchr_m, na.rm = TRUE),
            calchr_sd = sd(calchr_m, na.rm = TRUE),
            MR_mean = mean(MR_avg, na.rm = TRUE),
            MR_sd = sd(MR_avg, na.rm = TRUE),
            prop_ECG_m = mean(prop_ECG),
            prop_ECG_sd = sd(prop_ECG))
#### Assessing Post-Release Behaviour ####

postrelease.act<-dat %>% #filtering to the first 3 weeks (studies suggest that surgical effects may last up to 3wks)
  filter (numdays < 22)

postavgtemp<-postrelease.act %>% 
  group_by(numdays) %>% 
  summarise(timestamp = mean(timestamp),
            temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            ODBA_ms2 = mean(ODBA_ms2),
            ODBA_sum = sum(ODBA_ms2_sum),
            ODBA_max = max(ODBA_ms2_max),
            ODBA_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(numdays,fill=list(temp_m=NA,depth_m=NA,ODBA_ms2=NA, ODBA_sum=NA, ODBA_max=NA, ODBA_min=NA))

avgpostcomparison<-dat %>% #using the fourth week as a reference to compare against
  filter(numdays > 21 & numdays < 30) %>% 
  group_by(id) %>% 
  summarise(timestamp = mean(timestamp),
            temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            ODBA_ms2 = mean(ODBA_ms2),
            ODBA_sum = sum(ODBA_ms2_sum),
            ODBA_max = max(ODBA_ms2_max),
            ODBA_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(id,fill=list(temp_m=NA,depth_m=NA,ODBA_ms2=NA, ODBA_sum=NA, ODBA_max=NA, ODBA_min=NA)) %>% 
  mutate(highlight = "TRUE",.after="id") %>%
  mutate(numdays = 22,.before = "id")

postrelease.act<-postrelease.act %>% 
  group_by(id,numdays) %>%
  summarise(timestamp = mean(timestamp),
            temp_m = mean(temp_m, na.rm = TRUE),
            depth_m = mean(depth_m, na.rm = TRUE),
            ODBA_ms2 = mean(ODBA_ms2),
            ODBA_sum = sum(ODBA_ms2_sum),
            ODBA_max = max(ODBA_ms2_max),
            ODBA_min = min(ODBA_ms2_min))%>% 
  ungroup() %>%
  complete(numdays,id,fill=list(temp_m=NA,depth_m=NA,ODBA_ms2=NA, ODBA_sum=NA, ODBA_max=NA, ODBA_min=NA)) %>%
  mutate(highlight = "FALSE",.after="id") %>%
  union_all(avgpostcomparison)

#statistics

#stats post release behaviour

glimpse(postrelease.act)

postrelease.act$numdays<-as.character(postrelease.act$numdays)

hist(postrelease.act$ODBA_ms2)
hist(postrelease.act$depth_m)
ggqqplot(postrelease.act$ODBA_ms2)


ggqqplot(postrelease.act$depth_m)

plotdist(postrelease.act$ODBA_m, histo = TRUE, demp = TRUE)
descdist(postrelease.act$ODBA_m, boot = 1000)

postODBAvnumday.lmer<-lmer(ODBA_ms2~numdays*temp_m + (1|id), data = postrelease.act)
postODBAvnumday.lmer1<-lmer(ODBA_ms2~numdays+temp_m + (1|id), data = postrelease.act)
postODBAvnumday.lmer2<-lmer(ODBA_ms2~numdays + (1|id), data = postrelease.act)

AIC(postODBAvnumday.lmer,postODBAvnumday.lmer1,postODBAvnumday.lmer2)

Anova(postODBAvnumday.lmer2, type = 2)

postdepthvnumday.lmer<-lmer(depth_m~numdays + (1|id), data = postrelease.act)
postdepthvnumday.lmer1<-lmer(depth_m~numdays*temp_m + (1|id), data = postrelease.act)
postdepthvnumday.lmer2<-lmer(depth_m~numdays+temp_m + (1|id), data = postrelease.act)
postdepthvtemp.lmer<-lmer(depth_m~temp_m + (1|id), data = postrelease.act)
AIC(postdepthvnumday.lmer,postdepthvnumday.lmer1,postdepthvnumday.lmer2, postdepthvtemp.lmer)
#postdepthvnumday.lmer is best based on AIC - no interaction therefore Anova is used (i.e. type II) rather than anova (type III)

Anova(postdepthvnumday.lmer2, type = 2)
postdepth.emm = emmeans(postdepthvnumday.lmer,~numdays)
options(max.print = 100000)
pairs(postdepth.emm, adjust = "fdr")
post.relposthoc<-cld((emmeans(postdepthvnumday.lmer, ~numdays)), Letters = letters, adjust = "fdr")

#post release behaviour 

postrelease.act$numdays<-as.character(postrelease.act$numdays)

postODBAvnumday.lmer<-lmer(ODBA_ms2~numdays*temp_m + (1|id), data = postrelease.act)
postODBAvnumday.lmer1<-lmer(ODBA_ms2~numdays+temp_m + (1|id), data = postrelease.act)
postODBAvnumday.lmer2<-lmer(ODBA_ms2~numdays + (1|id), data = postrelease.act)

AIC(postODBAvnumday.lmer,postODBAvnumday.lmer1,postODBAvnumday.lmer2)

Anova(postODBAvnumday.lmer1, type = 2)

#diel cycle

glimpse(diel_df)

max(diel_df$hour2)

glimpse(dat)

diel_df$weight<-as.numeric(diel_df$weight)
diel_df$hour2<-as.character(diel_df$hour2)

table(diel_df$hour2)

overal_diel_ODBA<-lmer(ODBA_ms2~hour2*temp_m*weight + (1|id), data = diel_df)
overal_diel_ODBA1<-lmer(ODBA_ms2~hour2+temp_m*weight + (1|id), data = diel_df)
overal_diel_ODBA2<-lmer(ODBA_ms2~hour2*temp_m+weight + (1|id), data = diel_df)
overal_diel_ODBA3<-lmer(ODBA_ms2~hour2+temp_m+weight + (1|id), data = diel_df)
overal_diel_ODBA4<-lmer(ODBA_ms2~hour2*temp_m + (1|id), data = diel_df)
overal_diel_ODBA5<-lmer(ODBA_ms2~hour2+temp_m + (1|id), data = diel_df)
overal_diel_ODBA6<-lmer(ODBA_ms2~hour2*weight + (1|id), data = diel_df)
overal_diel_ODBA7<-lmer(ODBA_ms2~hour2+weight + (1|id), data = diel_df)
overal_diel_ODBA8<-lmer(ODBA_m~hour2 + (1|id), data = diel_df)

AIC(overal_diel_ODBA,overal_diel_ODBA1,overal_diel_ODBA2,overal_diel_ODBA3,overal_diel_ODBA4,
    overal_diel_ODBA5,overal_diel_ODBA6,overal_diel_ODBA7,overal_diel_ODBA8)

Anova(overal_diel_ODBA8)

postODBA.emm <- emmeans(overal_diel_ODBA8,~hour2)
pairs(postODBA.emm, adjust = "fdr")
diel_ODBA.posthoc<-cld(postODBA.emm, Letters = letters, sort = FALSE)

diel_df$hour2<-as.character(diel_df$hour2)

diel_dat<-diel_df %>% 
  drop_na(speed_bls_m)

diel_dat$hour2<-as.character(diel_dat$hour2)

overal_diel_speed8<-lmer(speed_bls_m~hour2 + (1|id), data = diel_dat)
anova(overal_diel_speed8)

postspeed.emm <- emmeans(overal_diel_speed8,~hour2)
pairs(postspeed.emm, adjust = "fdr")

diel_speed.posthoc<-cld(postspeed.emm, Letters = letters, sort = FALSE)

overal_diel_depth8<-lmer(depth_m~hour2 + (1|id), data = diel_df)
Anova(overal_diel_depth8)

postdepth.emm <- emmeans(overal_diel_depth8,~hour2)
pairs(postdepth.emm, adjust = "fdr")

diel_depth.posthoc<-cld(postdepth.emm, Letters = letters, sort = FALSE)

#plotting

postrelease.act.plot<-
  ggplot() +
  geom_boxplot(data = postrelease.act,aes(x=numdays,y=ODBA_ms2, group = numdays, colour = highlight),
               fill="white", outlier.shape = NA) +
  geom_jitter(data = postrelease.act,aes(x=numdays,y=ODBA_ms2, group = numdays, colour = highlight),
              alpha = 0.3, width=0.3,height=0,size=2) +
  labs(x="Days Post-Surgery",y=expression("ODBA"~(m~"\U00B7"~sec^{"-2"})), size = 6) +
  scale_y_continuous(breaks = seq(0.15, 0.3, 0.05), limits = c(0.15, 0.3)) + 
  scale_x_continuous(breaks = seq(1,22,1), limits = c(0.5,22.5),
                     labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22-29"))+ 
  scale_colour_manual(values = c("black","darkred")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"), 
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 5)

#average depth

postrelease.depth.plot<-ggplot() +
  geom_boxplot(data = postrelease.act,aes(x=numdays,y=(0-depth_m), group = numdays, colour = highlight),
               fill="white", outlier.shape = NA) +
  geom_jitter(data = postrelease.act,aes(x=numdays,y=(0-depth_m), group = numdays, colour = highlight),
              alpha = 0.3, width=0.3,height=0,size=2) +
  labs(x="Days Post-Surgery",y="Depth (m) ", size = 6) +
  scale_y_continuous(limits = c(-6, 0),
                     breaks = seq(-6,0,1),
                     labels = c("6", "5", "4", "3", "2", "1", "0"))+
  scale_x_continuous(breaks = seq(1,22,1),
                     limits = c(0.5,22.5),
                     labels = c("","2","","4","","6","","8","","10","","12","","14","","16","","18","","20","","22-29"))+
  scale_colour_manual(values = c("black","darkred")) +
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"),
        legend.position = "none")+
  stat_summary(fun = mean, colour = "darkred", geom= "point", shape = 18, size = 5)


####Assessing Relationships Between Temp. + DO and Behaviour####

library(lme4)
library(lmerTest)
library(MuMIn) 

#simulating data 
mass <- expand.grid(weight= seq(500, 1500, 25))
temps <- expand.grid(temp_m= seq(3, 27, 0.1))
sim_data<- merge(mass, temps)

sim_weight<-data.frame(weight = (1:100))

#temp vs depth
tempdepth_df<-dat_filled %>% 
  dplyr::select(id,timestamp,depth_m,temp_m,weight) %>% 
  na.omit()

tempdepth_df$weight<-as.numeric(tempdepth_df$weight)

require(MuMIn)
global_model_w_tempdepth<-lmer((depth_m)~log(temp_m)*log(weight)+(1|id), data = tempdepth_df) #log-log relationship is better
options(na.action = "na.fail")
combinations_w_tempdepth<-dredge(global_model_w_tempdepth)
print(combinations_w_tempdepth)

red_model_temp_depth<-lmer(log(depth_m)~log(temp_m):log(weight)+(1|id), data = tempdepth_df)  
red_model_temp_depth1<-lmer(log(depth_m)~log(temp_m)+(1|id), data = tempdepth_df) #going with this model - best explains the data
red_model_temp_depth2<-lmer(log(depth_m)~log(temp_m)+log(weight)+(1|id), data = tempdepth_df) 

AIC(red_model_temp_depth,red_model_temp_depth1,red_model_temp_depth2)

performance::check_model(red_model_temp_depth1)
require("ggResidpanel")
ggResidpanel::resid_panel(red_model_temp_depth1, smoother = TRUE, qqbands = TRUE, type = "pearson")

anova(red_model_temp_depth1)
r.squaredGLMM(red_model_temp_depth1)
confint(red_model_temp_depth1)

sim_data$pred_depth_temp<-exp(predict(red_model_temp_depth1, sim_data, re.form= ~0))
tempdepth_df$pred<-exp(predict(red_model_temp_depth1, tempdepth_df, re.form = ~0))

#temp v ODBA

tempODBA_df<-dat_filled %>% 
  dplyr::select(id,timestamp,ODBA_ms2,speed_bls,temp_m, weight, length) %>% 
  na.omit()

tempODBA_df$weight<-as.numeric(tempODBA_df$weight)

glimpse(tempODBA_df)

require(MuMIn)
global_model_w_tempODBA<-lmer((ODBA_ms2)~temp_m*log(weight)+(1|id), data = tempODBA_df) 
options(na.action = "na.fail")
combinations_w_tempODBA<-dredge(global_model_w_tempODBA)
print(combinations_w_tempODBA)

test_model<-lmer((ODBA_ms2)~temp_m+(1|id), data = tempODBA_df) 
test_model1<-lmer((ODBA_ms2)~temp_m+log(weight)+(1|id), data = tempODBA_df) 
test_model2<-lmer((ODBA_ms2)~temp_m+log(length)+(1|id), data = tempODBA_df)  

AIC(test_model,test_model1,test_model2)

anova(test_model)
r.squaredGLMM(test_model)
confint(test_model)

performance::check_model(red_model_temp_ODBA)
require("ggResidpanel")
ggResidpanel::resid_panel(red_model_temp_ODBA, smoother = TRUE, qqbands = TRUE, type = "pearson")

anova(red_model_temp_ODBA)
r.squaredGLMM(red_model_temp_ODBA)
confint(red_model_temp_ODBA)

glimpse(sim_data)

sim_data$pred_ODBA_temp<-(predict(test_model, sim_data, re.form= ~0))
tempODBA_df$pred<-(predict(test_model, tempODBA_df, re.form= ~0))

#DO v ODBA

DOODBA_df<-dat_filled %>% 
  dplyr::select(id,hr = timestamp,ODBA_ms2,temp_m,depth_m,weight, length) %>% 
  na.omit()

DOODBA_df$date<-as.numeric(strftime(DOODBA_df$hr, format = "%j"))
DOODBA_df$weight<-as.numeric(DOODBA_df$weight)

DOODBA_df$DO_m <- predict(DOForest,DOODBA_df)

# https://www.waterontheweb.org/under/waterquality/oxygen.html - formula here

DOODBA_df$bar<-(DOODBA_df$depth_m*997.05*9.80665)/100000 + 1.013

DOODBA_df$C<-exp(7.7117 - 1.31403*log(DOODBA_df$temp_m + 45.93))

DOODBA_df$Pw<-11.8571-(3840.7/DOODBA_df$temp_m)-(216961/(DOODBA_df$temp_m^2))

DOODBA_df$O<-0.000975-((1.426*10^(-5))*DOODBA_df$temp_m)+((6.436*10^(-8))*DOODBA_df$temp_m)

DOODBA_df$Cp<-(DOODBA_df$C*DOODBA_df$bar)*
  (((1-DOODBA_df$Pw/DOODBA_df$bar)*(1-DOODBA_df$O*DOODBA_df$bar))/
     ((1-DOODBA_df$Pw)*(1-DOODBA_df$O)))

DOODBA_df$DO_mgL<-(DOODBA_df$DO_m*DOODBA_df$Cp)/100

DOODBA_df$DO_m_scaled<-DOODBA_df$DO_m/100

DOvODBA.lmer3<-lmer((ODBA_ms2)~(DO_mgL)+(1|id), data = DOODBA_df) #going with linear - no weight is best
DOvODBA.lmer4<-lmer((ODBA_ms2)~(DO_mgL)*log(weight)+(1|id), data = DOODBA_df) 
DOvODBA.lmer5<-lmer((ODBA_ms2)~(DO_mgL)+log(weight)+(1|id), data = DOODBA_df) 
DOvODBA.lmer6<-lmer((ODBA_ms2)~(DO_mgL):log(weight)+(1|id), data = DOODBA_df) 
DOvODBA.lmer7<-lmer((ODBA_ms2)~(DO_mgL)+(1|id), data = DOODBA_df) 

AIC(DOvODBA.lmer3,DOvODBA.lmer4,DOvODBA.lmer5,DOvODBA.lmer6,DOvODBA.lmer7)

anova(DOvODBA.lmer3)
summary(DOvODBA.lmer3)
r.squaredGLMM(DOvODBA.lmer3)

DOODBA_df$pred<-predict(DOvODBA.lmer3, DOODBA_df, re.form=~0)

performance::check_model(DOvODBA.lmer3)
require("ggResidpanel")
ggResidpanel::resid_panel(DOvODBA.lmer3, smoother = TRUE, qqbands = TRUE, type = "pearson")

#DO v depth

DOdepth_df<-dat_filled %>% 
  dplyr::select(id,hr = timestamp,ODBA_ms2,temp_m,depth_m, weight) %>% 
  na.omit() 

DOdepth_df$date<-as.numeric(strftime(DOdepth_df$hr, format = "%j"))

DOdepth_df$DO_m <- predict(DOForest,DOdepth_df)

DOdepth_df$bar<-(DOdepth_df$depth_m*997.05*9.80665)/100000 + 1.013

DOdepth_df$C<-exp(7.7117 - 1.31403*log(DOdepth_df$temp_m + 45.93))

DOdepth_df$Pw<-11.8571-(3840.7/DOdepth_df$temp_m)-(216961/(DOdepth_df$temp_m^2))

DOdepth_df$O<-0.000975-((1.426*10^(-5))*DOdepth_df$temp_m)+((6.436*10^(-8))*DOdepth_df$temp_m)

DOdepth_df$Cp<-(DOdepth_df$C*DOdepth_df$bar)*
  (((1-DOdepth_df$Pw/DOdepth_df$bar)*(1-DOdepth_df$O*DOdepth_df$bar))/
     ((1-DOdepth_df$Pw)*(1-DOdepth_df$O)))

DOdepth_df$DO_mgL<-(DOdepth_df$DO_m*DOdepth_df$Cp)/100

DOvdepth.lmer2<-lmer(log(depth_m)~(DO_mgL)+(1|id), data = DOdepth_df)
DOvdepth.lmer3<-lmer(log(depth_m)~(DO_mgL)+log(weight)+(1|id), data = DOdepth_df)
DOvdepth.lmer4<-lmer(log(depth_m)~(DO_mgL):log(weight)+(1|id), data = DOdepth_df)
DOvdepth.lmer5<-lmer(log(depth_m)~(DO_mgL)+(1|id), data = DOdepth_df)

AIC(DOvdepth.lmer2,DOvdepth.lmer3,DOvdepth.lmer4,DOvdepth.lmer5) #temperature alone is best

anova(DOvdepth.lmer2)
summary(DOvdepth.lmer2)
r.squaredGLMM(DOvdepth.lmer2)

DOdepth_df$pred<-exp(predict(DOvdepth.lmer2,DOdepth_df,re.form = ~0))

performance::check_model(DOvdepth.lmer2)
require("ggResidpanel")
ggResidpanel::resid_panel(DOvdepth.lmer2, smoother = TRUE, qqbands = TRUE, type = "pearson")

#swimming speed vs temperature

glimpse(tempODBA_df)

sum(is.na(tempODBA_df$speed_bls))

#log speed wont work due to zeros - adding small value - linear relationship appears better anyways 

tempODBA_df$speed_bls<-tempODBA_df$speed_bls+0.0001

tempvspeed.lmer<-lmer(log(speed_bls)~(temp_m)*log(weight)+(1|id), data = tempODBA_df)
tempvspeed.lmer1<-lmer(log(speed_bls)~(temp_m)+log(weight)+(1|id), data = tempODBA_df)
tempvspeed.lmer2<-lmer(log(speed_bls)~(temp_m):log(weight)+(1|id), data = tempODBA_df)
tempvspeed.lmer3<-lmer(log(speed_bls)~(temp_m)+(1|id), data = tempODBA_df)

AIC(tempvspeed.lmer,tempvspeed.lmer1,tempvspeed.lmer2,tempvspeed.lmer3)

anova(tempvspeed.lmer1)
r.squaredGLMM(tempvspeed.lmer1) 
summary(tempvspeed.lmer1)#poly is marginally better (hard to compare by AIC) but an exponential relationship is more often used for swimming v temp

performance::check_model(tempvspeed.lmer1)
require("ggResidpanel")
ggResidpanel::resid_panel(tempvspeed.lmer1, smoother = TRUE, qqbands = TRUE, type = "pearson")

sim_data$pred_speed_temp<-(predict(tempvspeed.lmer1, sim_data, re.form = ~0))
tempODBA_df$pred_speed<-exp(predict(tempvspeed.lmer1, tempODBA_df, re.form = ~0))

#DO v Speed
DOODBA_df$speed_bls<-10^(-0.3846)+0.0273*DOODBA_df$temp_m+0.0452*(log10(DOODBA_df$ODBA_ms2)*log10(DOODBA_df$weight)*log10(DOODBA_df$length))

DOODBA_df<-DOODBA_df %>% 
  mutate(speed_bls = if_else(speed_bls < 0, 0, speed_bls))



DOvspeed.lmer<-lmer((speed_bls)~log(DO_mgL)*log(weight)+(1|id), data = DOODBA_df)
DOvspeed.lmer1<-lmer((speed_bls)~log(DO_mgL):log(weight)+(1|id), data = DOODBA_df)
DOvspeed.lmer2<-lmer((speed_bls)~log(DO_mgL)+log(weight)+(1|id), data = DOODBA_df)
DOvspeed.lmer3<-lmer((speed_bls)~log(DO_mgL)+(1|id), data = DOODBA_df)

AIC(DOvspeed.lmer,DOvspeed.lmer1,DOvspeed.lmer2,DOvspeed.lmer3)

anova(DOvspeed.lmer3)
summary(DOvspeed.lmer3)
r.squaredGLMM(DOvspeed.lmer3)

performance::check_model(DOvspeed.lmer3)
require("ggResidpanel")
ggResidpanel::resid_panel(DOvspeed.lmer3, smoother = TRUE, qqbands = TRUE, type = "pearson")

DOODBA_df$pred_speed<-(predict(DOvspeed.lmer3, DOODBA_df, re.form = ~0))

#plotting 

library("jcolors")

temp.x.depth_plot<-ggplot()+
  geom_point(data = tempdepth_df,aes(x=temp_m,y=depth_m), fill="lightgrey", colour = "lightgrey", size = 0.8, alpha = 0.5) +
  #geom_point(data = tempdepth_df,aes(x=temp_m,y=depth_m, colour = weight), size = 0.8, alpha = 0.5) +
  #geom_smooth(data = sim_data,aes(x=temp_m,y=pred_depth_temp),colour = "red", size = 0.8)+
  geom_smooth(data = tempdepth_df,aes(x=temp_m,y=pred),colour = "red", size = 0.8)+
  #geom_point(data = sim_data,aes(x=temp_m,y=pred_depth_temp, colour = weight), size = 1.5)+
  # geom_line(data = tempdepth, aes(x = temp_m, y = lmer_m), colour = "red", size = 1) + 
  # geom_ribbon(data = tempdepth, aes(ymin = lower, ymax = upper, x = temp_m, y = lmer_m), fill = "salmon", size = 0.8, alpha = 0.6) +
  scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 1)+
  theme_classic()+
  #scale_colour_viridis_c(option = "plasma")+
  scale_y_reverse(limits = c(6, -0.5), breaks = seq(0,6,1))+
  scale_x_continuous(breaks = seq(0,30,4), limits = c(2,30))+
  annotate(geom="text", x=(6.5), y=-0.25, label = "italic(R^2)[Marginal] == 0.31", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(18.5), y=-0.25, label = "italic(R^2)[Conditional] == 0.48", parse = TRUE, size = 4.5)+
  #annotate(geom="text", x=(1.85), y=3.2, label = "italic(R)^2[Marginal] == 0.72", parse = TRUE, size = 4)+
  #annotate(geom="text", x=(2), y=(-0.1), label = "n == 16", parse = TRUE, size = 8)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = "Temperature (°C)", y = "Depth (m)")

max(tempODBA_df$ODBA_ms2)
min(tempODBA_df$ODBA_ms2)

glimpse(tempODBA_df)

temp.x.ODBA_plot<-ggplot()+
  geom_point(data = tempODBA_df,aes(x=temp_m,y=ODBA_ms2), fill="lightgrey", colour = "lightgrey", size = 0.8, alpha = 0.5) +
  geom_smooth(data = tempODBA_df,aes(x=temp_m,y=pred),colour = "red", size = 0.8)+
  # geom_line(data = tempODBA, aes(x = temp_m, y = lmer_m), colour = "red", size = 1) + 
  # geom_ribbon(data = tempODBA, aes(ymin = lower, ymax = upper, x = temp_m, y = lmer_m), fill = "salmon", size = 0.8, alpha = 0.6) +
  theme_classic()+
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0,1.3,0.25))+
  scale_x_continuous(breaks = seq(0,30,4), limits = c(2,30))+
  annotate(geom="text", x=(6.5), y=1.27, label = "italic(R^2)[Marginal] == 0.16", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(18.5), y=1.2725, label = "italic(R^2)[Conditional] == 0.27", parse = TRUE, size = 4.5)+
  #annotate(geom="text", x=(1.85), y=3.2, label = "italic(R)^2[Marginal] == 0.72", parse = TRUE, size = 4)+
  #annotate(geom="text", x=(2), y=(-0.1), label = "n == 16", parse = TRUE, size = 8)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"))+
  labs(x = "Temperature (°C)", y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))

DO.x.ODBA_plot<-ggplot()+
  geom_point(data = DOODBA_df,aes(x=DO_mgL,y=ODBA_ms2), fill="lightgrey", colour = "lightgrey", size = 0.8, alpha = 0.5) +
  geom_smooth(data = DOODBA_df,aes(x=DO_mgL,y=pred),colour = "red", size = 0.8)+
  # geom_line(data = DOODBA, aes(x = DO_m, y = lmer_m), colour = "red", size = 1) + 
  # geom_ribbon(data = DOODBA, aes(ymin = lower, ymax = upper, x = DO_m, y = lmer_m), fill = "salmon", size = 0.8, alpha = 0.6) +
  theme_classic()+
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0,1.3,0.25))+
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12))+
  annotate(geom="text", x=(2), y=1.27, label = "italic(R^2)[Marginal] == 0.01", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(7), y=1.27, label = "italic(R^2)[Conditional] == 0.19", parse = TRUE, size = 4.5)+
  #annotate(geom="text", x=(1.85), y=3.2, label = "italic(R)^2[Marginal] == 0.72", parse = TRUE, size = 4)+
  #annotate(geom="text", x=(2), y=(-0.1), label = "n == 16", parse = TRUE, size = 8)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"))+
  labs(x = expression("Dissolved Oxygen"~(mg~"\U00B7"~L^{"-1"})), y = expression("ODBA"~(m~"\U00B7"~sec^{"-2"})))

DO.x.depth_plot<-ggplot()+
  geom_point(data = DOdepth_df,aes(x=DO_mgL,y=depth_m), fill="lightgrey", colour = "lightgrey", size = 0.8, alpha = 0.5) +
  geom_smooth(data = DOdepth_df,aes(x=DO_mgL,y=pred),colour = "red", size = 0.8)+
  # geom_line(data = DOdepth, aes(x = DO_m, y = lmer_m), colour = "red", size = 1) + 
  # geom_ribbon(data = DOdepth, aes(ymin = lower, ymax = upper, x = DO_m, y = lmer_m), fill = "salmon", size = 0.8, alpha = 0.6) +
  theme_classic()+
  scale_y_reverse(limits = c(6, -0.5), breaks = seq(0,6,1))+
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12))+
  annotate(geom="text", x=(2), y=-0.25, label = "italic(R^2)[Marginal] == 0.08", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(7), y=-0.25, label = "italic(R^2)[Conditional] == 0.33", parse = TRUE, size = 4.5)+
  #annotate(geom="text", x=(1.85), y=3.2, label = "italic(R)^2[Marginal] == 0.72", parse = TRUE, size = 4)+
  # annotate(geom="text", x=(2), y=(-0.1), label = "n == 16", parse = TRUE, size = 8)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = expression("Dissolved Oxygen"~(mg~"\U00B7"~L^{"-1"})), y = "Depth (m)")

temp.x.speed_plot<-ggplot()+
  geom_point(data = tempODBA_df,aes(x=temp_m,y=speed_bls), fill="lightgrey", colour = "lightgrey", size = 0.8, alpha = 0.5) +
  geom_smooth(data = tempODBA_df,aes(x=temp_m,y=pred_speed),colour = "red", size = 0.8)+
  # geom_line(data = tempspeed, aes(x = temp_m, y = lmer_m), colour = "red", size = 1) + 
  # geom_ribbon(data = tempspeed, aes(ymin = lower, ymax = upper, x = temp_m, y = lmer_m), fill = "salmon", size = 0.8, alpha = 0.6) +
  theme_classic()+
  scale_y_continuous(breaks = seq(0,2,0.5))+
  coord_cartesian(ylim = c(0, 2), xlim = c(2,30))+
  scale_x_continuous(breaks = seq(0,30,4))+
  annotate(geom="text", x=(6.5), y=1.88, label = "italic(R^2)[Marginal] == 0.92", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(18.5), y=1.88, label = "italic(R^2)[Conditional] == 0.93", parse = TRUE, size = 4.5)+
  #annotate(geom="text", x=(1.85), y=3.2, label = "italic(R)^2[Marginal] == 0.72", parse = TRUE, size = 4)+
  #annotate(geom="text", x=(2), y=(-0.1), label = "n == 16", parse = TRUE, size = 8)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = "Temperature (°C)", y = expression("Swimming Speed"~(BL~"\U00B7"~sec^{"-1"})))

DO.x.speed_plot<-ggplot()+
  geom_point(data = DOODBA_df,aes(x=DO_mgL,y=speed_bls), fill="lightgrey", colour = "lightgrey", size = 0.8, alpha = 0.5) +
  geom_smooth(data = DOODBA_df,aes(x=DO_mgL,y=pred_speed),colour = "red", size = 0.8)+
  # geom_line(data = DOspeed, aes(x = DO_m, y = lmer_m), colour = "red", size = 1) +
  # geom_ribbon(data = DOspeed, aes(ymin = lower, ymax = upper, x = DO_m, y = lmer_m), fill = "salmon", size = 0.8, alpha = 0.6) +
  # geom_smooth(data = DOODBA_df,aes(x=DO_m,y=speed_bls), method = "loess", se = FALSE)+
  theme_classic()+
  scale_y_continuous(breaks = seq(0,1.25,0.25))+
  coord_cartesian(ylim = c(0, 1.25), xlim = c(0,12))+
  scale_x_continuous(breaks = seq(0,12,2))+
  annotate(geom="text", x=(2), y=1.2, label = "italic(R^2)[Marginal] == 0.23", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(7), y=1.2, label = "italic(R^2)[Conditional] == 0.32", parse = TRUE, size = 4.5)+
  #annotate(geom="text", x=(1.85), y=3.2, label = "italic(R)^2[Marginal] == 0.72", parse = TRUE, size = 4)+
  #annotate(geom="text", x=(2), y=(-0.1), label = "n == 16", parse = TRUE, size = 8)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = "Dissolved Oxygen (%)", y = expression("Swimming Speed"~(BL~"\U00B7"~sec^{"-1"})))


(temp.x.depth_plot | temp.x.ODBA_plot)/(DO.x.depth_plot | DO.x.ODBA_plot)

hr.temp.gg<-ggplot() +
  geom_point(data = df_calcHR,aes(x=temp,y=calchr), fill="lightgrey", colour = "lightgrey", size = 0.8, alpha = 0.5) +
  geom_line(data = tempHR, aes(x = temp, y = lmer_m), colour = "red", size = 1) + 
  #geom_ribbon(data = tempHR, aes(ymin = lower, ymax = upper, x = temp, y = lmer_m), fill = "salmon", size = 0.8, alpha = 0.6) +
  # annotate(geom="text", x=1.55, y=100, label = "n = 21", size = 5)+
  annotate(geom="text", x=(6.5), y=100, label = "italic(R^2)[Marginal] == 0.58", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(18.5), y=100, label = "italic(R^2)[Conditional] == 0.68", parse = TRUE, size = 4.5)+
  coord_cartesian(ylim = c(0, 100), xlim = c(2,30))+
  scale_y_continuous(breaks = seq(0,100,20))+
  scale_x_continuous(breaks = seq(0,30,4))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = " Temperature (°C)", y = "Heart Rate (bpm)")

(hr.temp.gg/temp.x.depth_plot/temp.x.speed_plot/temp.x.ODBA_plot)|(plot_spacer()/DO.x.depth_plot/DO.x.speed_plot/DO.x.ODBA_plot)



####Data for Bioenergetics Simulations####

#use data averaged to the hour, add up each hour within a day 
#convert to mgO2 / g / hr then to calories using the Oxycal coefficient 3.24 cal mg 02-1

dat_averages<-dat %>% 
  group_by(timestamp) %>% 
  summarise(depth_mean = mean(depth_m,na.rm = TRUE),
            temp_mean = mean(temp_m, na.rm = TRUE))

dat_filled_MR<-dat %>% 
  dplyr::select(id,timestamp,temp_m,hour,depth_m,ODBA_ms2, speed_bls) %>% 
  left_join(dat_averages, by = "timestamp") %>% 
  mutate(depth_m = if_else(is.na(depth_m), depth_mean, depth_m)) %>% 
  mutate(temp_m = if_else(is.na(temp_m), temp_mean, temp_m)) %>% 
  group_by(id,hour) %>% 
  summarise(mean_depth = mean(depth_m),
            mean_temp = mean(temp_m),
            mean_ODBA = mean(ODBA_ms2),
            mean_speed_bls = mean(speed_bls)) %>% 
  dplyr::select(id,timestamp = hour,temp_m = mean_temp, depth_m = mean_depth, ODBA_ms2 = mean_ODBA, speed_bls= mean_speed_bls)

dat_month<-dat %>% 
  dplyr::select(id, timestamp, month_character)

dat_filled_MR<-dat_filled_MR %>% 
  left_join(weight_df, by = "id") %>% 
  left_join(dat_month, by = c("id", "timestamp"))

dat_filled_MR$speed_bls<-10^(1.3156+0.0277*dat_filled_MR$temp_m+0.3423*log10(dat_filled_MR$ODBA_ms2)-0.6505*log10(dat_filled_MR$length))

dat_filled_MR<-dat_filled_MR %>% 
  mutate(speed_bls = if_else(speed_bls < 0, 0, speed_bls)) %>% 
  mutate(speed_bls = if_else(ODBA_ms2 < 0.1, 0, speed_bls))

#predicting MO2 then converting to gO2/g/day
glimpse(dat_filled_MR)

dat_filled_MR$weight<-as.numeric(dat_filled_MR$weight)
dat_filled_MR$weight<-as.numeric(dat_filled_MR$length)

dat_filled_MR$weight<-as.numeric(dat_filled_MR$weight)

dat_filled_MR$MR<-10^(4.4489+0.0426*(dat_filled_MR$temp_m)-1.1446*log10(dat_filled_MR$length)+0.3321*log10(dat_filled_MR$ODBA_ms2))

dat_filled_MR<-dat_filled_MR %>% 
  mutate(MR = if_else(ODBA_ms2 < 0.1, 10^(4.4489+0.0426*(temp_m)-1.1446*log10(length)+0.3321*log10(ODBA_ms2)), MR))

dat_filled_MR$MR_bio<-((dat_filled_MR$MR*24)/1000)/1000

#speed needs to be in cm/s not BL/s for bioenergetics models
dat_filled_MR$speed_abs<-dat_filled_MR$speed_bls*(dat_filled_MR$length/10)

dat_filled_MR$date<-as.Date(dat_filled_MR$timestamp)

dat_filled_MR_avg<-dat_filled_MR %>% 
  group_by(id, date) %>% 
  summarise(MR_bio_day = mean(MR_bio,na.rm = TRUE),
            swim_day = mean(speed_abs, na.rm = TRUE),
            temp_m = mean(temp_m, na.rm = TRUE),
            weight = mean(weight)) %>% 
  mutate(MR_oxycal = MR_bio_day*13560*weight)#oxycal coefficient  to convert joules to calories # see Rice et al 1983

#heart rate

glimpse(datHR_hr)

datHR_hr$MR<-10^(9.4369-3.4368*log10(datHR_hr$length)+0.9120*log10(datHR_hr$calchr_m))
datHR_hr$MR_bio<-((datHR_hr$MR*24)/1000)/1000 #g O2/ g/day

datHR_hr$date<-as.Date(datHR_hr$hour)

datHR_bio<-datHR_hr %>% 
  group_by(id = tag, date) %>% 
  summarise(MR_bio_day = mean(MR_bio,na.rm = TRUE),
            temp_m = mean(temp_m, na.rm = TRUE),
            weight = mean(weight)) %>% 
  mutate(MR_oxycal = MR_bio_day*13560*weight)


dat_filled_MR_combined<-dat_filled_MR_avg %>% 
  rbind(datHR_bio) %>% 
  dplyr::select(-c(weight)) %>% 
  left_join(dat_morphs, by = "id") %>% 
  mutate(fall_day1 = if_else(start_date < as.Date("2021-03-03"), as.numeric(difftime(start_date, as.Date("2020-09-01"), units = "days")),
                             if_else(start_date > as.Date("2021-03-03") & start_date < as.Date("2022-09-01"),as.numeric(difftime(start_date, as.Date("2021-09-01"), units = "days")),
                                     as.numeric(difftime(start_date, as.Date("2022-09-01"), units = "days"))))) %>% 
  mutate(std_end_date = if_else(start_date < as.Date("2021-03-03"), as.numeric(difftime(recap_date, as.Date("2020-09-01"), units = "days")),
                                if_else(start_date > as.Date("2021-03-03") & start_date < as.Date("2022-09-01"),as.numeric(difftime(recap_date, as.Date("2021-09-01"), units = "days")),
                                        as.numeric(difftime(recap_date, as.Date("2022-09-01"), units = "days"))))) %>%
  mutate(std_diff = std_end_date - fall_day1)


#selecting fish with fall release and spring recapture to standardize weight loss

winter_loss<-dat_morphs %>% 
  filter(tagged_duration < 365 & release_day > 50 & recap_day < 200) %>% 
  summarise(weight_loss_m = mean(weight_change,na.rm = TRUE),
            weight_loss_sd = sd(weight_change,na.rm = TRUE),
            prop_weight_loss_m = mean(prop_weight_change,na.rm = TRUE),
            prop_weight_loss_sd = sd(prop_weight_change,na.rm = TRUE))


#Calculating Bioenergetic Coefficient####

#first trying to determine ACT from field data by relating absolute swimming speed to temperature and weight

glimpse(dat_filled_MR)

min(dat_filled_MR$speed_abs,na.rm = TRUE)

ACT_lmer1<-lmer(log(speed_abs+0.0001)~log(weight)+temp_m+(1|id), data = dat_filled_MR)#speed cannot be zero so I need to shift

anova(ACT_lmer1)

summary(ACT_lmer1) #rewritten as 1.4343*(W^0.3162)*(e^0.05918*temp), ACT = 1.4343, BACT = 0.06904, RK4 = 0.3162
confint(ACT_lmer1)
r.squaredGLMM(ACT_lmer1) #0.61, 0.61

mass <- expand.grid(weight= seq(500, 1500, 25))
temps <- expand.grid(temp_m= seq(3, 27, 0.1))
sim_data<- merge(mass, temps)

sim_data$test_Vel<-0.3792*sim_data$weight^(0.5366)*exp(0.05918*sim_data$temp_m)

pred<-dat_filled_MR
pred$pred<-exp(predict(ACT_lmer1, re.form = ~0, test)) #looks pretty good

ggplot()+
  geom_point(data = dat_filled_MR, mapping = aes(x = temp_m, y = speed_abs))+
  geom_smooth(data = test, mapping = aes(x = temp_m, y = pred))+
  theme_classic()

act.x.speed_plot<-ggplot()+
  geom_point(data = dat_filled_MR, mapping = aes(x = temp_m, y = speed_abs), fill="lightgrey", colour = "lightgrey", size = 0.8, alpha = 0.5) +
  geom_smooth(data = pred, mapping = aes(x = temp_m, y = pred),colour = "red", size = 0.8)+
  theme_classic()+
  scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,4))+
  annotate(geom="text", x=(8.8), y=97, label = "italic(R^2)[Marginal] == 0.61", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(9.5), y=90, label = "italic(R^2)[Conditional] == 0.61", parse = TRUE, size = 4.5)+
  #annotate(geom="text", x=(1.85), y=3.2, label = "italic(R)^2[Marginal] == 0.72", parse = TRUE, size = 4)+
  #annotate(geom="text", x=(2), y=(-0.1), label = "n == 16", parse = TRUE, size = 8)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"))+
  labs(x = "Temperature (°C)", y = expression("Swimming Speed"~~(cm~"\U00B7"~sec^{"-1"})))

act.x.speed_sim<-ggplot(sim_data, aes(x = temp_m, y = test_Vel, colour = weight/1000))+
  geom_point()+
  theme_classic()+
  scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,4))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "Temperature (°C)", y = expression("Swimming Speed"~~(cm~"\U00B7"~sec^{"-1"})), colour = "Mass (kg)")

act.x.speed_plot|act.x.speed_sim

####Bioenergetics Modeling####
####Bioenergetics Analyses####
#Loading in Data 
oxycal <- 13560 #standard in FB4
EDP <- 4000 #prey energy density (consider altering)

param <- fread("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/Parameters_official_updatedLMB.csv")
lmb <- param %>% filter(Species=="Largemouth bass (adult)")
lmb_updated<-param %>% filter(Species=="Largemouth bass (updated)")
lmb$ED <- as.numeric(lmb$ED)
lmb_updated$ED <- as.numeric(lmb_updated$ED)
lmb$RTL <- as.numeric(lmb$RTL)
lmb_updated$RTL <- as.numeric(lmb_updated$RTL)

glimpse(lmb_updated)

#data 
fb <- merge(data.frame(temp=seq(1,30,1)), data.frame(weight=seq(500,2000,10)))

fb_ACT <- merge(data.frame(temp=seq(1,30,1)), data.frame(weight=seq(500,2000,10)))

fb_updated <- merge(data.frame(temp=seq(1,30,1)), data.frame(weight=seq(500,2000,10)))
head(fb)


lmb_updated$RTL<-20
lmb_updated$RK1<-5.7058
lmb_updated$RB<-lmb$RB

lmb_act<-lmb

lmb_act$RTL<-lmb_updated$RTL
lmb_act$RK1<-lmb_updated$RK1
lmb_act$RK4<-lmb_updated$RK4
lmb_act$ACT<-lmb_updated$ACT
lmb_act$BACT<-lmb_updated$BACT

####Updated Model####

#metabolism 
#EQ1

#below are respiration equations added from FB4

fb_updated$ft <- exp(lmb_updated$RQ*fb_updated$temp)

glimpse(lmb_updated)

glimpse(fb_updated)

#Equation 1 ### Temperature function equation 1 (Hanson et al. 1997; Stewart et al. 1983)

fb_updated<-fb_updated %>% 
  mutate(VEL = if_else(temp <= lmb_updated$RTL, (lmb_updated$ACT * weight ^ lmb_updated$RK4 * exp(lmb_updated$BACT * temp)), #this describes when RTL is noted
                       (lmb_updated$RK1 * weight ^ lmb_updated$RK4 * exp(lmb_updated$RK5 * temp))))

fb_updated$ACTIVITY<-exp(lmb_updated$RTO*fb_updated$VEL)

fb_updated$Rmax <- lmb_updated$RA * fb_updated$weight ^ lmb_updated$RB

fb_updated$Met.J <- fb_updated$Rmax * fb_updated$ft *fb_updated$weight*fb_updated$ACTIVITY * oxycal

#metabolism values
ggplot(data = fb_updated, aes(x = temp, y = ft))+geom_point()+
  ggplot(data = fb_updated, aes(temp, Rmax, col=weight))+geom_point()+
  ggplot(data = fb_updated, aes(temp, Met.J, col=weight))+geom_point()+
  ggplot(data = fb_updated, aes(temp, ACTIVITY, col=weight))+geom_point()

#consumption
#EQ2 
lmb_updated$CY <- log(lmb_updated$CQ) * (lmb_updated$CTM - lmb_updated$CTO + 2)
lmb_updated$CZ <- log(lmb_updated$CQ) * (lmb_updated$CTM - lmb_updated$CTO)
lmb_updated$CX <- (lmb_updated$CZ^2 * (1+(1+40/lmb_updated$CY)^0.5)^2)/400

fb_updated$V <- (lmb_updated$CTM - fb_updated$temp) / (lmb_updated$CTM - lmb_updated$CTO)
fb_updated$ft_c <- ifelse(fb_updated$temp < lmb_updated$CTM, fb_updated$V^lmb_updated$CX * exp(lmb_updated$CX * (1 - fb_updated$V)), 0)
fb_updated$Cmax <- lmb_updated$CA * (fb_updated$weight) ^lmb_updated$CB
fb_updated$Cons.p <- 0.75 #max consumption at 1
fb_updated$C <- fb_updated$Cmax * fb_updated$Cons.p * fb_updated$ft_c
fb_updated$Cons.g <- fb_updated$C*fb_updated$weight
fb_updated$Cons.J <- fb_updated$Cons.g*EDP #prey energy density

#consumption values:
ggplot(fb_updated, aes(temp, ft_c))+geom_point()+
  ggplot(fb_updated, aes(temp, Cmax, col=weight))+geom_point()+
  ggplot(fb_updated, aes(temp, C, col=weight))+geom_point()+
  ggplot(fb_updated, aes(temp, Cons.J, col=weight))+geom_point()

#wastes
#EG equation 1
fb_updated$Eg = lmb_updated$FA*fb_updated$Cons.J 
fb_updated$Ex = lmb_updated$UA*(fb_updated$Cons.J-fb_updated$Eg)

#SDA
fb_updated$SDA <- lmb_updated$SDA *(fb_updated$Cons.J-fb_updated$Eg) 

#growth
fb_updated$growth.J <- fb_updated$Cons.J-fb_updated$Met.J-fb_updated$SDA-fb_updated$Eg-fb_updated$Ex
fb_updated$growth.g <- fb_updated$growth.J/lmb_updated$ED

test<-fb_updated %>% 
  filter(temp < 10)

#outputs
head(fb_updated)
ggplot(fb_updated, aes(temp, Cons.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_updated, aes(temp, Met.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_updated, aes(temp, VEL, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_updated, aes(temp, SDA, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_updated, aes(temp, Eg, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_updated, aes(temp, Ex, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_updated, aes(temp, growth.g, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()

#outputs plots 
cons.g_updated<-ggplot(fb_updated, aes(temp, Cons.g, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
 scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Consumption"~~(g~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

cons_updated<-ggplot(fb_updated, aes(temp, Cons.J, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Consumption"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")


met_updated<-ggplot(fb_updated, aes(temp, Met.J, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
 scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Respiration"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

vel_updated<-ggplot(fb_updated, aes(temp, VEL, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = "Temperature (°C)", y = expression("Swimming Speed"~~(cm~"\U00B7"~sec^{"-1"})), colour = "Mass (kg)")

sda_updated<-ggplot(fb_updated, aes(temp, SDA, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
   scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("SDA"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

eg_updated<-ggplot(fb_updated, aes(temp, Eg, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Egestion"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

ex_updated<-ggplot(fb_updated, aes(temp, Ex, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Excretion"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

growth.g_updated<-ggplot(fb_updated, aes(temp, growth.g, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Growth"~~(g~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

(cons.g_updated|cons_updated|met_updated|vel_updated)/(sda_updated|eg_updated|ex_updated|growth.g_updated)


#grow function 
#grow a fish over period of time 

#set up a data frame with temp and fish size
temp <- read.csv("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/Temperature.csv")
temp.loess <- loess(temperature~day, data=temp)
sim <- data.frame(day=seq(1:366))
sim$temp <- predict(temp.loess, sim)
ggplot(sim, aes(day, temp))+geom_point()

#set up weights at start, grow function will update them over time
sim$weight <- NA
sim$weight[1] <- 500
sim$E <- NA
sim$E[1] <- sim$weight[1]*as.numeric(lmb_updated$ED)
head(sim)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow(data=sim, meta=lmb_updated, ndays=365,  p=0.75) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(day, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(day, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(day, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(day, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()


#use fit.p function to estimate consumption by known growth rate 
fit.p_updated <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow(data=sim, meta=lmb_updated, ndays=365,  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow(data=sim, meta=lmb_updated, ndays=365,  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_updated(IW = 500, 
                   FW = 600, 
                   p = 0.2,  #starting value
                   W.tol = 0.0001, #tolerance for difference from value
                   max.iter = 25) #max iterations

p #can now plug this value into the grow function 

####Old Model####
fb$ft <- exp(lmb$RQ*fb$temp)

fb<-fb %>%
  mutate(VEL = if_else(temp <= lmb$RTL, (lmb$ACT * weight ^ lmb$RK4 * exp(lmb$BACT * temp)), #this describes when RTL is noted
                       (lmb$RK1 * weight ^ lmb$RK4 * exp(lmb$RK5 * temp))))

fb$ACTIVITY<-exp(lmb$RTO*fb$VEL)

fb$Rmax <- lmb$RA * fb$weight ^ lmb$RB

fb$Met.J <- fb$Rmax * fb$ft *fb$weight*fb$ACTIVITY * oxycal


#metabolism values
ggplot(data = fb, aes(x = temp, y = ft))+geom_point()+
  ggplot(data = fb, aes(temp, Rmax, col=weight))+geom_point()+
  ggplot(data = fb, aes(temp, Met.J, col=weight))+geom_point()

#consumption
#EQ2 
lmb$CY <- log(lmb$CQ) * (lmb$CTM - lmb$CTO + 2)
lmb$CZ <- log(lmb$CQ) * (lmb$CTM - lmb$CTO)
lmb$CX <- (lmb$CZ^2 * (1+(1+40/lmb$CY)^0.5)^2)/400

fb$V <- (lmb$CTM - fb$temp) / (lmb$CTM - lmb$CTO)
fb$ft_c <- ifelse(fb$temp < lmb$CTM, fb$V^lmb$CX * exp(lmb$CX * (1 - fb$V)), 0)
fb$Cmax <- lmb$CA * (fb$weight) ^lmb$CB
fb$Cons.p <- 0.75 #max consumption at 1
fb$C <- fb$Cmax * fb$Cons.p * fb$ft_c
fb$Cons.g <- fb$C*fb$weight
fb$Cons.J <- fb$Cons.g*EDP #prey energy density

#consumption values:
ggplot(fb, aes(temp, ft_c))+geom_point()+
  ggplot(fb, aes(temp, Cmax, col=weight))+geom_point()+
  ggplot(fb, aes(temp, C, col=weight))+geom_point()+
  ggplot(fb, aes(temp, Cons.J, col=weight))+geom_point()

#wastes
#EG equation 1
fb$Eg = lmb$FA*fb$Cons.J 
fb$Ex = lmb$UA*(fb$Cons.J-fb$Eg)

#SDA
fb$SDA <- lmb$SDA *(fb$Cons.J-fb$Eg) 

#growth
fb$growth.J <- fb$Cons.J-fb$Met.J-fb$SDA-fb$Eg-fb$Ex
fb$growth.g <- fb$growth.J/lmb$ED

#outputs
head(fb)
ggplot(fb, aes(temp, Cons.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb, aes(temp, Met.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb, aes(temp, VEL, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb, aes(temp, SDA, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb, aes(temp, Eg, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb, aes(temp, Ex, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb, aes(temp, growth.g, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()

#outputs plots 
cons.g_original<-ggplot(fb, aes(temp, Cons.g, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  # scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Consumption"~~(g~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

cons_original<-ggplot(fb, aes(temp, Cons.J, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  # scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Consumption"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")


met_original<-ggplot(fb, aes(temp, Met.J, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  # scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Respiration"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

vel_original<-ggplot(fb, aes(temp, VEL, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  # scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  labs(x = "Temperature (°C)", y = expression("Swimming Speed"~~(cm~"\U00B7"~sec^{"-1"})), colour = "Mass (kg)")

sda_original<-ggplot(fb, aes(temp, SDA, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  # scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("SDA"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

eg_original<-ggplot(fb, aes(temp, Eg, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  # scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Egestion"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

ex_original<-ggplot(fb, aes(temp, Ex, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  # scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Excretion"~~(J~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

growth.g_original<-ggplot(fb, aes(temp, growth.g, col=weight/1000))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  # scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(), 
        legend.position = "none")+
  labs(x = "Temperature (°C)", y = expression("Growth"~~(g~"\U00B7"~day^{"-1"})), colour = "Mass (kg)")

(cons.g_original|cons_original|met_original|vel_original)/(sda_original|eg_original|ex_original|growth.g_original)

ggsave(("D:/2021 Recovered Tags - Manual Download/Figures-2024/LMB_original_growth_outputs2024.png"),dpi = 600, width = 13, height = 6, units = "in")


#grow function 
#grow a fish over period of time 

#set up a data frame with temp and fish size
temp <- read.csv("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/Temperature.csv")
temp.loess <- loess(temperature~day, data=temp)
sim <- data.frame(day=seq(1:366))
sim$temp <- predict(temp.loess, sim)
ggplot(sim, aes(day, temp))+geom_point()

#set up weights at start, grow function will update them over time
sim$weight <- NA
sim$weight[1] <- 500
sim$E <- NA
sim$E[1] <- sim$weight[1]*as.numeric(lmb$ED)
head(sim)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
output_lmb <- grow(data=sim, meta=lmb, ndays=365,  p=0.5865) #proportion of cmax
energetics_lmb <- as.data.frame(output_lmb[[1]])
FW <- output_lmb[[2]]
head(energetics_lmb)
FW

ggplot(energetics_lmb, aes(day, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb, aes(day, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb, aes(day, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb, aes(day, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()


#use fit.p function to estimate consumption by known growth rate 
#use fit.p function to estimate consumption by known growth rate 
fit.p_lmb <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow(data=sim, meta=lmb, ndays=365,  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow(data=sim, meta=lmb, ndays=365,  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_lmb(IW = 800, 
               FW = 1000, 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations

p #can now plug this value into the grow function 

#set up a data frame with temp and fish size
temp <- read.csv("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/Temperature.csv")
temp.loess <- loess(temperature~day, data=temp)
sim <- data.frame(day=seq(1:366))
sim$temp <- predict(temp.loess, sim)
ggplot(sim, aes(day, temp))+geom_point()

#set up weights at start, grow function will update them over time
sim$weight <- NA
sim$weight[1] <- 500
sim$E <- NA
sim$E[1] <- sim$weight[1]*as.numeric(lmb$ED)
head(sim)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
output_lmb <- grow(data=sim, meta=lmb, ndays=365,  p=0.5865) #proportion of cmax
energetics_lmb <- as.data.frame(output_lmb[[1]])
FW <- output_lmb[[2]]
head(energetics_lmb)
FW

ggplot(energetics_lmb, aes(day, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb, aes(day, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb, aes(day, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb, aes(day, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()


#use fit.p function to estimate consumption by known growth rate 
#use fit.p function to estimate consumption by known growth rate 
fit.p_lmb <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow(data=sim, meta=lmb, ndays=365,  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow(data=sim, meta=lmb, ndays=365,  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_lmb(IW = 800, 
               FW = 1000, 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations

####Old Model w Act####
fb_ACT$ft <- exp(lmb_act$RQ*fb_ACT$temp)

fb_ACT<-fb_ACT %>%
  mutate(VEL = if_else(temp <= lmb_act$RTL, (lmb_act$ACT * weight ^ lmb_act$RK4 * exp(lmb_act$BACT * temp)), #this describes when RTL is noted
                       (lmb_act$RK1 * weight ^ lmb_act$RK4 * exp(lmb_act$RK5 * temp))))

fb_ACT$ACTIVITY<-exp(lmb_act$RTO*fb_ACT$VEL)

fb_ACT$Rmax <- lmb_act$RA * fb_ACT$weight ^ lmb_act$RB

fb_ACT$Met.J <- fb_ACT$Rmax * fb_ACT$ft *fb_ACT$weight*fb_ACT$ACTIVITY * oxycal


#metabolism values
ggplot(data = fb_ACT, aes(x = temp, y = ft))+geom_point()+
  ggplot(data = fb_ACT, aes(temp, Rmax, col=weight))+geom_point()+
  ggplot(data = fb_ACT, aes(temp, Met.J, col=weight))+geom_point()

#consumption
#EQ2 
lmb_act$CY <- log(lmb_act$CQ) * (lmb_act$CTM - lmb_act$CTO + 2)
lmb_act$CZ <- log(lmb_act$CQ) * (lmb_act$CTM - lmb_act$CTO)
lmb_act$CX <- (lmb_act$CZ^2 * (1+(1+40/lmb_act$CY)^0.5)^2)/400

fb_ACT$V <- (lmb_act$CTM - fb_ACT$temp) / (lmb_act$CTM - lmb_act$CTO)
fb_ACT$ft_c <- ifelse(fb_ACT$temp < lmb_act$CTM, fb_ACT$V^lmb_act$CX * exp(lmb_act$CX * (1 - fb_ACT$V)), 0)
fb_ACT$Cmax <- lmb_act$CA * (fb_ACT$weight) ^lmb_act$CB
fb_ACT$Cons.p <- 0.3 #max consumption at 1
fb_ACT$C <- fb_ACT$Cmax * fb_ACT$Cons.p * fb_ACT$ft_c
fb_ACT$Cons.g <- fb_ACT$C*fb_ACT$weight
fb_ACT$Cons.J <- fb_ACT$Cons.g*EDP #prey energy density

#consumption values:
ggplot(fb_ACT, aes(temp, ft_c))+geom_point()+
  ggplot(fb_ACT, aes(temp, Cmax, col=weight))+geom_point()+
  ggplot(fb_ACT, aes(temp, C, col=weight))+geom_point()+
  ggplot(fb_ACT, aes(temp, Cons.J, col=weight))+geom_point()

#wastes
#EG equation 1
fb_ACT$Eg = lmb_act$FA*fb_ACT$Cons.J 
fb_ACT$Ex = lmb_act$UA*(fb_ACT$Cons.J-fb_ACT$Eg)

#SDA
fb_ACT$SDA <- lmb_act$SDA *(fb_ACT$Cons.J-fb_ACT$Eg) 

#growth
fb_ACT$growth.J <- fb_ACT$Cons.J-fb_ACT$Met.J-fb_ACT$SDA-fb_ACT$Eg-fb_ACT$Ex
fb_ACT$growth.g <- fb_ACT$growth.J/lmb_act$ED

#outputs
head(fb_ACT)
ggplot(fb_ACT, aes(temp, Cons.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_ACT, aes(temp, Met.J, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_ACT, aes(temp, VEL, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_ACT, aes(temp, SDA, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_ACT, aes(temp, Eg, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_ACT, aes(temp, Ex, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(fb_ACT, aes(temp, growth.g, col=weight))+geom_point()+scale_color_viridis_c()+theme_bw()

####Applying Models to Individuals####
Warner_data<-fread("D:\\2021 Recovered Tags - Manual Download\\Bioenergetics Analysis\\bioenergetic_temperature_MR_speed_updated2024.csv")

glimpse(Warner_data)

#adding numdays
Warner_data<-Warner_data %>% 
  mutate(numdays = as.numeric(difftime(date,start_date,units = "days")))

#fish caught within 30days of tag recording ending

Warner_close<-Warner_data %>% 
  filter(close_recap == "yes")

#Tag_1####

glimpse(lmb_updated)

Warner_Tag_1<-Warner_data %>% 
  filter(id == "Tag_1") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_1 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date, swim_day)

glimpse(data)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.3164
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, swim_day, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()


#comparing with updated bioenergetics model

data<-Warner_Tag_1 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J_recorded = MR_oxycal, ndays = numdays, date)#need to change ndays to day

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.709) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, VEL, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.8200718) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g),
            cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_1"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_1"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_1"

wint_cons_combined<-rbind(wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)


#Tag_2 ####

Warner_Tag_2<-Warner_data %>% 
  filter(id == "Tag_2") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_2 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.3309
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.6098251) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.790192) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

Tag_2_ind_cons<-ggplot(energetics_ind, aes(date, Cons.J, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Predicted Metabolism")+
  scale_y_continuous(breaks = seq(0,25000,5000), limits = c(3000,25000), labels = label_comma())+
  annotate("text", x = as.Date("2021-03-01"), y = 25000, label = "p [Consumption] == 0.33", parse = TRUE, size = 3.5)+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Consumption"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_2_ind_met<-ggplot(energetics_ind, aes(date, Met.J, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Predicted Metabolism")+
  scale_y_continuous(breaks = seq(0,17500,2500), limits = c(4000,17500), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Respiration"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_2_ind_weight<-ggplot(energetics_ind, aes(date, weight, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Predicted Metabolism")+
  scale_y_continuous(breaks = seq(920,1040,20), limits = c(930,1040), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = "Mass (g)", x = "Date", colour = "Temperature (°C)")

Tag_2_upd_cons<-ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Updated Coefficients")+
  scale_y_continuous(breaks = seq(0,40000,5000), limits = c(6000,40000), labels = label_comma())+
  annotate("text", x = as.Date("2021-03-01"), y = 40000, label = "p [Consumption] == 0.57", parse = TRUE, size = 3.5)+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Consumption"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_2_upd_met<-ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Updated Coefficients")+
  scale_y_continuous(breaks = seq(0,25000,5000), limits = c(5000,25000), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Respiration"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_2_upd_weight<-ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Updated Coefficients")+
  scale_y_continuous(breaks = seq(940,1040,10), limits = c(940,1040), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = "Mass (g)", x = "Date", colour = "Temperature (°C)")

Tag_2_orig_cons<-ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Original Model")+
  annotate("text", x = as.Date("2021-03-01"), y = 50000, label = "p [Consumption] == 0.79", parse = TRUE, size = 3.5)+
  scale_y_continuous(breaks = seq(0,50000,10000), limits = c(5000,50000), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Consumption"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_2_orig_met<-ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Original Model")+
  scale_y_continuous(breaks = seq(10000,16000,2000), limits = c(10000,16000), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Respiration"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_2_orig_weight<-ggplot(energetics_lmb_original, aes(date, weight, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Original Model")+
  scale_y_continuous(breaks = seq(860,1040,40), limits = c(870,1040), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  labs(y = "Mass (g)", x = "Date", colour = "Temp. (°C)")


(Tag_2_orig_cons|Tag_2_orig_met|Tag_2_orig_weight)/(Tag_2_upd_cons|Tag_2_upd_met|Tag_2_upd_weight)/(Tag_2_ind_cons|Tag_2_ind_met|Tag_2_ind_weight)

ggsave(("D:/2021 Recovered Tags - Manual Download/Figures-2024/fish1_bioenergetic_example_updated2024.png"),dpi = 600, width = 12, height = 8, units = "in")

#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_2"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_2"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_2"

wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)


#Tag_3 ####

Warner_Tag_3<-Warner_data %>% 
  filter(id == "Tag_3") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_3 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.4602 - this model suggests most of the weight lost was in the spring - this is an issue... 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.4228531) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.6501209) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_3"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_3"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_3"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)


#Tag_4####

Warner_Tag_4<-Warner_data %>% 
  filter(id == "Tag_4") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations 

data<-Warner_Tag_4 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.6476 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.6042061) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.7530842) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_4"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_4"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_4"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)

#Tag_5 ####

Warner_Tag_5<-Warner_data %>% 
  filter(id == "Tag_5") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations 

data<-Warner_Tag_5 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.5103 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.4836349) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.6900005) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()


#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_5"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_5"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_5"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)


#Tag_10 ####

Warner_Tag_10<-Warner_data %>% 
  filter(id == "Tag_10") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations 

data<-Warner_Tag_10 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0.2615456) #proportion of cmax = 0.6571 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.484318) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.4272884) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_10"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_10"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_10"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)

#Tag_11 ####

Warner_Tag_11<-Warner_data %>% 
  filter(id == "Tag_11") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations 

data<-Warner_Tag_11 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.8095 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.5119871) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.4830769) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()


#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_11"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_11"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_11"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)

#Tag_12 ####

Warner_Tag_12<-Warner_data %>% 
  filter(id == "Tag_12") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations 

data<-Warner_Tag_12 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.877 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

energetics_ind<-energetics_ind %>% filter(ndays<132)

Tag_12_ind_cons<-ggplot(energetics_ind, aes(date, Cons.J, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Predicted Metabolism")+
  scale_y_continuous(breaks = seq(0,100000,20000), limits = c(5000,100000), labels = label_comma())+
  annotate("text", x = as.Date("2021-03-01"), y = 100000, label = "p [Consumption] == 0.88", parse = TRUE, size = 3.5)+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Consumption"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_12_ind_met<-ggplot(energetics_ind, aes(date, Met.J, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Predicted Metabolism")+
  scale_y_continuous(breaks = seq(0,80000,20000), limits = c(5000,80000), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Respiration"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_12_ind_weight<-ggplot(energetics_ind, aes(date, weight, col=temp))+
  geom_point()+
  scale_color_viridis_c()+
  theme_classic()+
  ggtitle("Predicted Metabolism")+
  scale_y_continuous(breaks = seq(810,850,10), limits = c(810,852), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = "Mass (g)", x = "Date", colour = "Temperature (°C)")

Tag_12_ind_cons|Tag_12_ind_met|Tag_12_ind_weight

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.6307693) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

energetics_lmb_updated<-energetics_lmb_updated %>% filter(ndays<132)

Tag_12_upd_cons<-ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+
  geom_point()+
  scale_color_viridis_c(breaks=seq(5,25,5), limits = c(4,25))+
  theme_classic()+
  ggtitle("Updated Coefficients")+
  scale_y_continuous(breaks = seq(0,80000,10000), limits = c(4000,80000), labels = label_comma())+
  annotate("text", x = as.Date("2021-03-01"), y = 80000, label = "p [Consumption] == 0.50", parse = TRUE, size = 3.5)+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Consumption"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_12_upd_met<-ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+
  geom_point()+
  scale_color_viridis_c(breaks=seq(5,25,5), limits = c(4,25))+
  theme_classic()+
  ggtitle("Updated Coefficients")+
  scale_y_continuous(breaks = seq(0,40000,10000), limits = c(5000,40000), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Respiration"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_12_upd_weight<-ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+
  geom_point()+
  scale_color_viridis_c(breaks=seq(5,25,5), limits = c(4,25))+
  theme_classic()+
  ggtitle("Updated Coefficients")+
  scale_y_continuous(breaks = seq(800,850,10), limits = c(800,852), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = "Mass (g)", x = "Date", colour = "Temperature (°C)")

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.6113695) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

energetics_lmb_original<-energetics_lmb_original %>% filter(ndays<132)

Tag_12_orig_cons<-ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+
  geom_point()+
  scale_color_viridis_c(breaks=seq(5,25,5), limits = c(4,25))+
  theme_classic()+
  ggtitle("Original Model")+
  annotate("text", x = as.Date("2021-03-01"), y = 70000, label = "p [Consumption] == 0.61", parse = TRUE, size = 3.5)+
  scale_y_continuous(breaks = seq(0,70000,10000), limits = c(5000,70000), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Consumption"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_12_orig_met<-ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+
  geom_point()+
  scale_color_viridis_c(breaks=seq(5,25,5), limits = c(4,25))+
  theme_classic()+
  ggtitle("Original Model")+
  scale_y_continuous(breaks = seq(10000,20000,2500), limits = c(9500,20000), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  labs(y = expression("Respiration"~~(J~"\U00B7"~day^{"-1"})), x = "Date", colour = "Temperature (°C)")

Tag_12_orig_weight<-ggplot(energetics_lmb_original, aes(date, weight, col=temp))+
  geom_point()+
  scale_color_viridis_c(breaks=seq(5,25,5), limits = c(4,25))+
  theme_classic()+
  ggtitle("Original Model")+
  scale_y_continuous(breaks = seq(730,850,20), limits = c(740,852), labels = label_comma())+
  theme(plot.title = element_text(size = 12),
        axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  labs(y = "Mass (g)", x = "Date", colour = "Temp. (°C)")

(Tag_12_orig_cons|Tag_12_orig_met|Tag_12_orig_weight)/(Tag_12_upd_cons|Tag_12_upd_met|Tag_12_upd_weight)/(Tag_12_ind_cons|Tag_12_ind_met|Tag_12_ind_weight)

#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_12"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_12"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_12"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)

#Tag_13 ####

Warner_Tag_13<-Warner_data %>% 
  filter(id == "Tag_13") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations

data<-Warner_Tag_13 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.7328 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.52156) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.4809) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_13"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_13"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_13"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)

#Tag_14 ####

Warner_Tag_14<-Warner_data %>% 
  filter(id == "Tag_14") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations

data<-Warner_Tag_14 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

glimpse(data)
glimpse(sim)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.8339 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.4842) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.4657871) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"Tag_14"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"Tag_14"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"Tag_14"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)

#Heart Rate Tags ####

#HR1010 is the only one with consistent MR data 

#HR1010 ####

Warner_HR1010<-Warner_data %>% 
  filter(id == "HR1010") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations

data<-Warner_HR1010 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data$ndays<-data$ndays - min(data$ndays)

glimpse(data)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

#grow function
source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")

fit.p_ind <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_ind(IW = data$weight_start[1], 
               FW = data$weight_end[1], 
               p = 0.2,  #starting value
               W.tol = 0.0001, #tolerance for difference from value
               max.iter = 25) #max iterations
p

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p) #proportion of cmax = 0.5344 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

ggplot(energetics_ind, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Cons.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_ind, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with updated bioenergetics model

glimpse(data)

data$weight[1] <- as.numeric(data$weight_start[1])
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)
head(data)

#grow function
#use fit.p function to estimate consumption by known growth rate 
fit.p_sim_data <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data(IW = data$weight_start[1], 
                    FW = data$weight_end[1], 
                    p = 0.2,  #starting value
                    W.tol = 0.0001, #tolerance for difference from value
                    max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0.5362835) #proportion of cmax
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

ggplot(energetics_lmb_updated, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_updated, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()

#comparing with original bioenergetics model

glimpse(lmb)

fit.p_sim_data_orig <- function(p, IW, FW, W.tol, max.iter) {
  W      <- IW    # Initial weight
  n.iter <- 0     # Counter for number of iterations
  p.max  <- 1  # current max
  p.min  <- 0  # current min
  
  output <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
  W.p <- output[[2]]
  
  while((n.iter <= max.iter) & (abs(W.p-FW) > W.tol)) {
    n.iter <- n.iter + 1
    if(W.p > FW) {p.max <- p} else {p.min <- p}
    p <- (p.min + p.max)/2 #p.min + (p.max - p.min)/2
    g <-  grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p)
    W.p <- g[[2]]
  }
  return(p)
}  

p <- fit.p_sim_data_orig(IW = data$weight_start[1], 
                         FW = data$weight_end[1], 
                         p = 0.2,  #starting value
                         W.tol = 0.0001, #tolerance for difference from value
                         max.iter = 25) #max iterations
p

source("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/functions_updated.R")
outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0.6562) #proportion of cmax
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

ggplot(energetics_lmb_original, aes(date, Cons.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Met.J, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, Growth.g, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()+
  ggplot(energetics_lmb_original, aes(date, weight, col=temp))+geom_point()+scale_color_viridis_c()+theme_bw()


#average winter consumption 

glimpse(energetics_ind)

wint_cons_ind<-energetics_ind %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_ind$model<-"Ind"
wint_cons_ind$Tag<-"HR1010"

wint_cons_lmb_updated<-energetics_lmb_updated %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_updated$model<-"lmb_updated"
wint_cons_lmb_updated$Tag<-"HR1010"

wint_cons_lmb_original<-energetics_lmb_original %>% 
  filter(temp < 10) %>% 
  summarise(cons_m = mean(Cons.g),
            cons_sum = sum(Cons.g), cons_dur = (max(ndays) - min(ndays)))

wint_cons_lmb_original$model<-"lmb_original"
wint_cons_lmb_original$Tag<-"HR1010"


wint_cons_combined<-rbind(wint_cons_combined,wint_cons_ind,wint_cons_lmb_updated,wint_cons_lmb_original)

####Predicting Winter Consumption####

glimpse(wint_cons_combined)

glimpse(Warner_close)

Warner_weights<-Warner_close %>% 
  dplyr::select(Tag = id, weight_start) %>% 
  group_by(Tag) %>% 
  summarise(weight = mean(weight_start))

wint_cons_combined<-wint_cons_combined %>% 
  left_join(Warner_weights, by = "Tag")

ind_cons<-wint_cons_combined %>% 
  filter(model == "Ind") %>% 
  summarise(cons_mean = mean(cons_m),
            cons_sd = sd(cons_m),
            dur_mean = mean(cons_dur),
            weight_m = mean(weight)) %>% 
  mutate(cons_err = cons_sd/sqrt(11))

upd_cons<-wint_cons_combined %>% 
  filter(model == "lmb_updated") %>% 
  summarise(cons_mean = mean(cons_m),
            cons_sd = sd(cons_m),
            dur_mean = mean(cons_dur),
            weight_m = mean(weight)) %>% 
  mutate(cons_err = cons_sd/sqrt(11))

orig_cons<-wint_cons_combined %>% 
  filter(model == "lmb_original") %>% 
  summarise(cons_mean = mean(cons_m),
            cons_sd = sd(cons_m),
            dur_mean = mean(cons_dur),
            weight_m = mean(weight)) %>% 
  mutate(cons_err = cons_sd/sqrt(11))

#Estimating Losses Assuming Zero Consumption Below 10C -----------

#Tag_1 ####

#start weight = 1310, obs end = 1210 

#using calculated MR, pred end = 975.74, pred end_red = 1063
#using updated model, pred end = 986.73, 
#using original model, pred end = 882.21, 

glimpse(lmb_updated)

Warner_Tag_1<-Warner_data %>% 
  filter(id == "Tag_1") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_1 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#%loss 

(IW_obs-FW)/IW_obs

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

(IW_obs-FW)/IW_obs

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

(IW_obs-FW)/IW_obs


#Tag_2 -----

#start weight = 1040, obs end = 935 

#using calculated MR, pred end = 761.1, pred end_red = 835
#using updated model, pred end = 776.59, 
#using original model, pred end = 681.18, 

#note 3 days at start before 10C or less, starts at 12C

glimpse(lmb_updated)

Warner_Tag_2<-Warner_data %>% 
  filter(id == "Tag_2") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_2 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

(IW_obs-FW)/IW_obs

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW


#Tag_3 -----

#start weight = 1015, obs end = 829 

#using calculated MR, pred end = 758.98, pred end_red = 826.99
#using updated model, pred end = 760.06, 
#using original model, pred end = 664.76, 

#note 3 days at start before 10C or less, starts at 12C

glimpse(lmb_updated)

Warner_Tag_3<-Warner_data %>% 
  filter(id == "Tag_3") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_3 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#Tag_4 -----

#start weight = 1060, obs end = 952 

#using calculated MR, pred end = 762.37, pred end_red = 841.21
#using updated model, pred end = 786.24, 
#using original model, pred end = 692.56, 

#note 3 days at start before 10C or less, starts at 12C

glimpse(lmb_updated)

Warner_Tag_4<-Warner_data %>% 
  filter(id == "Tag_4") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_4 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#Tag_5 -----

#start weight = 1075, obs end = 903 

#using calculated MR, pred end = 811.56, pred end_red = 881.65
#using updated model, pred end = 811.91, 
#using original model, pred end = 715.38, 

#note 3 days at start before 10C or less, starts at 12C

glimpse(lmb_updated)

Warner_Tag_5<-Warner_data %>% 
  filter(id == "Tag_5") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_5 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#Tag_10 -----

#start weight = 1085, obs end = 954 

#using calculated MR, pred end = 938.13, pred end_red = 974.05
#using updated model, pred end = 937.23, 
#using original model, pred end = 881.02, 

#note 3 days at start before 10C or less, starts at 12C

glimpse(lmb_updated)

Warner_Tag_10<-Warner_data %>% 
  filter(id == "Tag_10") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_10 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

(IW_obs-FW)/IW_obs

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#Tag_11 -----

#start weight = 807, obs end = 762 

#using calculated MR, pred end = 672.84, pred end_red = 705.06
#using updated model, pred end = 692.46, 
#using original model, pred end = 645.27, 

#note 3 days at start before 10C or less, starts at 12C

glimpse(lmb_updated)

Warner_Tag_11<-Warner_data %>% 
  filter(id == "Tag_11") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_11 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

(IW_obs-FW)/IW_obs

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW


#Tag_12 -----

#start weight = 851, obs end = 837 

#using calculated MR, pred end = 698.62, pred end_red = 732.05
#using updated model, pred end = 717.88, 
#using original model, pred end = 669.35, 

#note 3 days at start before 10C or less, starts at 12C

glimpse(lmb_updated)

Warner_Tag_12<-Warner_data %>% 
  filter(id == "Tag_12") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_12 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#Tag_13 -----

#start weight = 905, obs end = 845 

#using calculated MR, pred end = 783.50, pred end_red = 812.20
#using updated model, pred end = 790.86, 
#using original model, pred end = 746.14, 

glimpse(lmb_updated)

Warner_Tag_13<-Warner_data %>% 
  filter(id == "Tag_13") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_13 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#Tag_14 -----

#start weight = 570, obs end = 571 

#using calculated MR, pred end = 477.24, pred end_red = 497.60
#using updated model, pred end = 504.30, 
#using original model, pred end = 475.80, 

glimpse(lmb_updated)

Warner_Tag_14<-Warner_data %>% 
  filter(id == "Tag_14") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations -

data<-Warner_Tag_14 %>% 
  dplyr::select(temp = temp_m, weight_start,weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#HR1010 -----

#start weight = 900, obs end = 777 

#using calculated MR, pred end = 596.85, pred end_red = 674.62
#using updated model, pred end = 604.18, 
#using original model, pred end = 503.38, 

glimpse(lmb_updated)

Warner_HR1010<-Warner_data %>% 
  filter(id == "HR1010") %>% 
  dplyr::select(id,date,numdays,start_date,end_recording,recap_date,weight_end = weight, weight_start, weight_change, temp_m,MR_oxycal,swim_day)

#using Updated Model

meta=lmb_updated

# Calculations

data<-Warner_HR1010 %>% 
  dplyr::select(temp = temp_m, weight_start, weight_end, Met.J = MR_oxycal,ndays = numdays, date)

data$ndays<-data$ndays - min(data$ndays)

glimpse(data)

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

data<-data %>% 
  mutate(Met.J = if_else(ndays > 50 & temp > 10, NA_real_, Met.J)) %>% 
  drop_na(Met.J)

glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

(IW_obs-FW)/IW_obs

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

#reduced metabolism due to fasting

data<-data %>% 
  mutate(Met.J_red = if_else(ndays > 14, Met.J*0.7, Met.J)) %>% 
  dplyr::select(temp, weight_start,weight_end, Met.J = Met.J_red,ndays, date)


glimpse(data)

ggplot(data, aes(x = date, y = temp))+
  geom_point()

data$weight[1] <- data$weight_start[1]
data$E <- NA
data$E[1] <- data$weight[1]*as.numeric(lmb_updated$ED)

head(data)

#grow function
FW_obs<-data$weight_end[1]
IW_obs<-data$weight_start[1]

#calculated

output_ind <- grow_ind(data=data, meta=lmb_updated, ndays=max(data$ndays),  p = 0) #assuming zero eating - how much is lost over the winter 
energetics_ind <- as.data.frame(output_ind[[1]])
FW <- output_ind[[2]]
head(energetics_ind)
FW

#updated

outputs_updated <- grow_ind_2(data=data, meta=lmb_updated, ndays=max(data$ndays),  p=0) 
energetics_lmb_updated <- as.data.frame(outputs_updated[[1]])
FW <- outputs_updated[[2]]
head(energetics_lmb_updated)
FW

#original

outputs_original <- grow_ind_2(data=data, meta=lmb, ndays=max(data$ndays),  p=0) 
energetics_lmb_original <- as.data.frame(outputs_original[[1]])
FW <- outputs_original[[2]]
head(energetics_lmb_original)
FW

####Assessing Relationship Between Mass, Weight Loss, and Consumption####

dat_morphs_combined<-fread("D:/2021 Recovered Tags - Manual Download/bioenergetic_morphometrics_updated2024.csv")

glimpse(dat_morphs_combined)

dat_morphs_combined$condition<-100*dat_morphs_combined$weight_start/(dat_morphs_combined$length/10)^3

#only assessing those individuals that were released in the fall and recapture in the spring
dat_morphs_close<-dat_morphs_combined %>% 
  filter(tagged_duration < 365 & release_day > 50 & recap_day < 200) 

condition_model<-lm(condition~prop_weight_change, data = dat_morphs_close)

summary(condition_model)
anova(condition_model)
confint(condition_model)

weight_model<-lm(weight_start~prop_weight_change, data = dat_morphs_close)

summary(weight_model)
anova(weight_model)
confint(weight_model)

condition_fig<-ggplot(dat_morphs_close, aes(x = condition, y = prop_weight_change))+
  geom_point()+
  geom_smooth(method = "lm",colour = "red", fill = "salmon", size = 0.8)+
  theme_classic()+
  coord_cartesian(ylim = c(-20, 25), xlim = c(1.18,1.93))+
  scale_y_continuous(breaks = seq(-20, 25, 10))+
  scale_x_continuous(breaks = seq(1.2,1.8,0.2), limits = c(1.18,1.93))+
  annotate(geom="text", x=(1.255), y=22, label = "italic(R^2) == 0.12", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(1.25), y=18, label = "italic(p) == 0.08", parse = TRUE, size = 4.5)+
  #annotate(geom="text", x=(1.85), y=3.2, label = "italic(R)^2[Marginal] == 0.72", parse = TRUE, size = 4)+
  #annotate(geom="text", x=(2), y=(-0.1), label = "n == 16", parse = TRUE, size = 8)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"))+
  labs(x = "Fall Condition", y = "Mass Change (% Body Mass)")

weight_fig<-ggplot(dat_morphs_close, aes(x = weight_start, y = prop_weight_change))+
  geom_point()+
  geom_smooth(method = "lm",colour = "red", fill = "salmon", size = 0.8)+
  theme_classic()+
  coord_cartesian(ylim = c(-20, 25))+
  scale_y_continuous(breaks = seq(-20, 25, 10))+
  # scale_x_continuous(breaks = seq(1.2,1.8,0.2), limits = c(1.18,1.93))+
  annotate(geom="text", x=(650), y=22, label = "italic(R^2) == 0.14", parse = TRUE, size = 4.5)+
  annotate(geom="text", x=(645), y=18, label = "italic(p) == 0.05", parse = TRUE, size = 4.5)+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "Fall Mass (g)", y = "Mass Change (% Body Mass)")

condition_fig|weight_fig

ggsave(("D:/2021 Recovered Tags - Manual Download/Figures-2024/Condition_Fig.png"),dpi = 600, width = 9, height = 5, units = "in")

#importing dataset produced from predictions of consumption generated from original, updated, and free-swimming models
bio_loss <- fread("D:/2021 Recovered Tags - Manual Download/Bioenergetics Analysis/Bioenergetics_Consumptions_Loss_updated.csv")

condition_df<-dat_morphs_close %>% 
  dplyr::select(Tag = id, condition)

bio_loss$cons_prop<-bio_loss$cons_winter/bio_loss$IW*100
bio_loss$condition<-100*bio_loss$IW/(bio_loss$length/10)^3


field_m1<-lm(cons_prop~IW, data = (bio_loss%>% filter(Model == "field")))
summary(field_m1)
anova(field_m1)
confint(field_m1)

lmb_updated_m1<-lm(cons_prop~IW, data = (bio_loss%>% filter(Model == "lmb_updated")))
summary(lmb_updated_m1)
confint(lmb_updated_m1)

lmb_original_m1<-lm(cons_prop~IW, data = (bio_loss%>% filter(Model == "lmb_original")))
summary(lmb_original_m1)
confint(lmb_original_m1)

field_m2<-lm(cons_prop~obs_loss, data = (bio_loss%>% filter(Model == "field")))
summary(field_m2)
anova(field_m2)
confint(field_m2)

lmb_updated_m2<-lm(cons_prop~obs_loss, data = (bio_loss%>% filter(Model == "lmb_updated")))
summary(lmb_updated_m2)
anova(lmb_updated_m2)
confint(lmb_updated_m2)

lmb_original_m2<-lm(cons_prop~obs_loss, data = (bio_loss%>% filter(Model == "lmb_original")))
summary(lmb_original_m2)
confint(lmb_original_m2)

field_m3<-lm(cons_prop~condition, data = (bio_loss%>% filter(Model == "field")))
summary(field_m3)
anova(field_m3)
confint(field_m3)

lmb_updated_m3<-lm(cons_prop~condition, data = (bio_loss%>% filter(Model == "lmb_updated")))
summary(lmb_updated_m3)
confint(lmb_updated_m3)

lmb_original_m3<-lm(cons_prop~condition, data = (bio_loss%>% filter(Model == "lmb_original")))
summary(lmb_original_m3)
confint(lmb_original_m3)

cons_IW<-ggplot()+
  geom_point(bio_loss, mapping = aes(x = IW, y = cons_prop, colour = factor(Model)))+
  geom_line(bio_loss, mapping = aes(x = IW, y = cons_prop, colour = factor(Model)), stat = "smooth", method = "lm", formula = y~x, alpha = 0.4, size = 1)+
  scale_y_continuous(limits = c(0,40), breaks = seq(0,40,5))+
  #scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"),
        legend.position = "none") +
  labs(x = "Fall Mass (g)", y = "Winter Consumption (% Body Mass)")+
  scale_colour_manual(values = c( "#0072B2","forestgreen", "#F21A00"), labels = c("Biologger","Original","Updated"))

cons_condition<-ggplot()+
  geom_point(bio_loss, mapping = aes(x = condition, y = cons_prop, colour = factor(Model)))+
  geom_line(bio_loss, mapping = aes(x = condition, y = cons_prop, colour = factor(Model)), stat = "smooth", method = "lm", formula = y~x, alpha = 0.4, size = 1)+
  scale_y_continuous(limits = c(0,40), breaks = seq(0,40,5))+
  scale_x_continuous(limits = c(1.3,2), breaks = seq(1.2,2,0.2))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(x = "Fall Condition", y = "Winter Consumption (% Body Mass)")+
  scale_colour_manual(values = c( "#0072B2","forestgreen", "#F21A00"), labels = c("Biologger","Original","Updated"))

cons_obsloss<-ggplot()+
  geom_point(bio_loss, mapping = aes(x = obs_loss, y = cons_prop, colour = factor(Model)))+
  geom_line(bio_loss, mapping = aes(x = obs_loss, y = cons_prop, colour = factor(Model)), stat = "smooth", method = "lm", formula = y~x, alpha = 0.4, size = 1)+
  # scale_y_continuous(limits = c(0.25,1.75), breaks = seq(0.25,1.75, 0.25))+
  # scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  theme_classic()+
  theme(axis.line = element_line(size = 0.7, linetype = "solid"),
        axis.ticks = element_line(colour = "black", size = 0.7),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()) +
  labs(x = "Mass Change (%)", y = "Winter Consumption (% Body Mass)")+
  scale_colour_manual(values = c( "#0072B2","forestgreen", "#F21A00"), labels = c("Biologger","Original","Updated"))

cons_IW|cons_condition|cons_obsloss

ggsave(("D:/2021 Recovered Tags - Manual Download/Figures-2024/Consumption_Fig_2.png"),dpi = 600, width = 11, height = 5, units = "in")


