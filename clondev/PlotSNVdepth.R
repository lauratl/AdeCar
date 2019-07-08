library(ggplot2)
library(tidyverse)

depth34 = read.table("07_WORKDIR/clondev/AC34.depth", stringsAsFactors = FALSE, header = FALSE)
colnames(depth34) = c("Adenoma", "Carcinoma")
depth34 = depth34 %>% gather(key = "sample", value = "depth", Adenoma, Carcinoma) %>% mutate(Patient="AC34")


depth33= read.table("07_WORKDIR/clondev/AC33.depth", stringsAsFactors = FALSE, header = FALSE)
colnames(depth33) = c("Adenoma", "Carcinoma")
depth33 = depth33 %>% gather(key = "sample", value = "depth", Adenoma, Carcinoma) %>% mutate(Patient="AC33")

ggplot(rbind(depth33,depth34)) + 
  geom_histogram(aes(x = depth, fill = sample), position = "dodge", binwidth = 5)+
  xlim(c(0,300))+
  scale_x_continuous(breaks = seq(0,300,20), limits = c(0,300))+
  facet_grid(sample~Patient)+
  ggtitle("Depth distribution (zoom at 0-300)")
