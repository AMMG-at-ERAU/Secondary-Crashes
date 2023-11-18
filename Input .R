library(dplyr)
Book1 <- crashdata_for_2015_2016_and_2017_sorted

Book1 %>% select(date, time, long, lat)->Book1


##################### to select data only in the selected region############################################
#Tampa
ROI_data <- Book1[Book1$lat >= 27.96 & Book1$lat <= 27.97 & Book1$long >= -82.4464 & Book1$long <= -82.4129,]
crash.cat1 <- catalog(ROI_data, time.begin="2015/01/01",
                      study.start="2015/01/01", study.end="2018/01/01",
                      lat.range = c(27.96, 27.97),  long.range = c(-82.4464, -82.4129))

crash.cat1
#Plant city
ROI_data <- Book1[Book1$lat >= 28.0274 & Book1$lat <= 28.0366 & Book1$long >= -82.1479 & Book1$long <= -82.1163,]
crash.cat1 <- catalog(ROI_data, time.begin="2015/01/01",
                      study.start="2015/01/01", study.end="2018/01/01",
                      lat.range = c(28.0274, 28.0366),  long.range = c(-82.1479, -82.1163))

crash.cat1
#Kissimmee
ROI_data <- Book1[Book1$lat >= 28.3025 & Book1$lat <= 28.3253 & Book1$long >= -81.5735 & Book1$long <= -81.5520,]
crash.cat1 <- catalog(ROI_data, time.begin="2015/01/01",
                      study.start="2015/01/01", study.end="2018/01/01",
                      lat.range = c(28.3025, 28.3253),  long.range = c(-81.5735, -81.5520))

crash.cat1
#Orlando
ROI_data <- Book1[Book1$lat >= 28.5821 & Book1$lat <= 28.5969 & Book1$long >= -81.3857 & Book1$long <= -81.3747,]
crash.cat1 <- catalog(ROI_data, time.begin="2015/01/01",
                      study.start="2015/01/01", study.end="2018/01/01",
                      lat.range = c(28.5821, 28.5969),  long.range = c(-81.3857, -81.3747))

crash.cat1
#Sanford
ROI_data <- Book1[Book1$lat >= 28.8236 & Book1$lat <= 28.8458 & Book1$long >= -81.3308 & Book1$long <= -81.3088,]
crash.cat1 <- catalog(ROI_data, time.begin="2015/01/01",
                      study.start="2015/01/01", study.end="2018/01/01",
                      lat.range = c(28.8236, 28.8458),  long.range = c(-81.3308, -81.3088))

crash.cat1
#Daytona beach
ROI_data <- Book1[Book1$lat >= 29.11 & Book1$lat <= 29.145 & Book1$long >= -81.16 & Book1$long <= -81.10,]
crash.cat1 <- catalog(ROI_data, time.begin="2015/01/01",
                      study.start="2015/01/01", study.end="2018/01/01",
                      lat.range = c(29.11, 29.145),  long.range = c(-81.16, -81.10))

crash.cat1
plot(crash.cat1)
##########################################################################################

param01 <- c(0.05, 0.22, 2.7)
crash.fit1 <- etas(crash.cat1, param0 = param01)
pbd1 <- probs(crash.fit1)
write.table(pbd1, file="Daytona_beach_2miles_3years_cos_function_rush_hour.csv", row.names=F, sep=",")
event <- crash.fit1$object$revents
write.table(event, file="revent_Daytona_beach_2miles_3years_cos_function_rush_hour.csv", row.names=F, sep=",")
