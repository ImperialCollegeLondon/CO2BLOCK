library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(rstatix)
library(readr)
require(scales)
data1 <- read_csv(file="R_S1csv.csv")
data1$Country <- as.factor(data1$Country)
data1$Global <- as.integer(data1$Global)
mycolors <- c("#E495A5", "#FFED6F","turquoise2","#80B1D3" ,"#86B875",  "#FDB462" ,"violetred1" ,"darkorchid1", "#ACA4E2",
              "tomato")
data2 <-filter(data1, Global%in% c(10))
data2$Country <- as.factor(data2$Country)

p1 <- ggplot(data2, aes(x = Growth, color = Country, y = Storage,size = Rate, shape = Type)) 

p1+ geom_point(alpha = 0.5, aes(color = Country)) + scale_colour_manual(values = mycolors)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  scale_size(range = c(1, 13))+ labs(title = "Global storage rate achieved in 2050: 10 Gt/yr", x = "Growth Rate [%]", y = "Storage resource required [Gt]") +
  ggplot2::scale_x_continuous(limits = c(0, 25))

# 612, 537 resolution

