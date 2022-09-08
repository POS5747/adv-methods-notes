
library(tidyverse)

data <- read_csv("https://raw.githubusercontent.com/pos3713ri/3713ri-anes/main/data/anes-2016.csv") %>%
  glimpse()

  glimpse() %>%
  ggplot(aes(x = ft_scientists)) +
  geom_histogram() 


data1 <- filter(data, religion == "Not Religious", age < 30) %>%
  select(ft_scientists) %>%
  na.omit() %>%
  glimpse()

datapasta::df_paste(data.frame(ft_scientists = data1$ft_scientists))
data.frame(
  ft_scientists = c(70,70,60,100,85,70,100,70,100,
                    85,70,40,50,70,90,100,100,85,60,60,70,70,100,
                    60,100,60,100,100,100,100,100,100,100,100,100,
                    85,60,100,85,85,50,100,100,97,85,50,60,100,100,
                    100,50,100,85,50,85,100,85,50,84,100,85,85,100,
                    99,96,73,70,85,100,100,100,40,60,63,85,100,60,
                    50,100,100,100,30,95,85,100,96,100,81,90,90,
                    100,50,100,95,100,50,100,83,70,100,30,70,23,72,
                    70,85,100,0,98,90,97,100,100,99,83,100,99,75,
                    100,100,86,85,85,100,83,51,90,100,96,100,50,95,
                    100,100,30,99,70,85,40,85,100,100,85,99,100,100,
                    80,85,52,70,90,85,85,100,85,90)
)