install.packages(c('gapminder','ggplot2','gganimate','gifski'))
install.packages('colorspace')
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(rstatix)
library(readr)
require(scales)
library(gganimate)
library(gifski)
data1 <- read_csv(file="R_animation.csv")
data1$Country <- as.factor(data1$Country)
data1$Global <- as.integer(data1$Global)
mycolors <- c("#E495A5", "#FFED6F","turquoise2","#80B1D3" ,"#86B875", "#5CBD92", "#FDB462" ,"violetred1" ,"darkorchid1", "#ACA4E2",
              "#CD99D8" ,"tomato")

p1 <- ggplot(data1, aes(x = Growth, color = Country, y = Storage, size = Rate, shape = Type)) +
  geom_point(alpha = 0.5) + scale_colour_manual(values = mycolors)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  scale_size(range = c(1, 10)) +
  labs(title = "Global storage rate achieved in 2050: {closest_state} GtCO2/yr", x = "Growth Rate [%]", y = "Storage resource required [Gt]") + transition_states(Global, transition_length = 1)  +
  ease_aes('linear') 
animate(p1)
anim_save("Country_distribution_compiled.gif")