install.packages("ggrepel")
install.packages("RColorbrewer")
install.packages("rstantools")
install.packages("pals")
install.packages("dplyr")
library(tidyverse)
library(ggstatsplot)
library(ggplot2)
library(rstatix)
library(readr)
library(dplyr)
install.packages(c('gapminder','ggplot2','gganimate','gifski'))
library(tidyverse)
require(scales)
library(gganimate)
library(gifski)
mycolors <- c("turquoise2","#80B1D3" , "#CD99D8" , "#FFED6F","#E495A5","turquoise2","#80B1D3" ,"#86B875", "#5CBD92", "#FDB462" ,"violetred1" ,"darkorchid1", "#ACA4E2",
              "#CD99D8" ,"tomato")

data1 <- read_csv(file="R_S1csv.csv")
data1$Global <- as.integer(data1$Global)
data2 <-filter(data1, Global%in% c(2,4,10,13))
data3 <-filter(data2, Country%in% c("China","EU","UK","Canada","Australia"))
p1 <- grouped_ggbetweenstats(
  data3,
  Country, 
  Rate, 
  grouping.var = Global,
  pairwise.comparisons = FALSE,
  color = Country,
  plotgrid.args = list(nrow = 4),
  results.subtitle = FALSE) + scale_colour_manual(values = mycolors)

plot(p1)

