install.packages("ggrepel")
install.packages("rstantools")
install.packages("pals")
install.packages("dplyr")
library(tidyverse)
library(ggstatsplot)
library(ggplot2)
library(rstatix)
library(readr)
library(RColorBrewer)
library(dplyr)
data1 <- read_csv(file="R_S1csv.csv")
data1$Global <- as.integer(data1$Global)
data1$Global <- as.factor(data1$Global)
data2 <-filter(data1, Global%in% c(2,6,10,13))

data3 <-filter(data2, Country == "US")

mycolors <- c("tomato","tomato","tomato","tomato")
p1 <- ggbetweenstats(
  data = data3,
  x = Global,
  y = Rate, 
  colour = Global,
  ylab = "Storage rate [Gt/yr]",
  xlab = "Global storage rate achieved in 2050 [Gt/yr]",
  title = "US",
  results.subtitle = FALSE,
  pairwise.comparisons = FALSE) + scale_colour_manual(values = mycolors)

plot(p1)


