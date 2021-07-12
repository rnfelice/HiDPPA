#plot results

library(tidyverse)
library(fs)
library(stringr)
library(ggthemes)
data_dir<-"./simulation_results"
resultsfiles <- dir(data_dir, full.names = TRUE)
my_results <- tibble(filename=resultsfiles)
my_results<- my_results %>%
  mutate(., n_Taxa = str_extract(filename, "(?<=results_).+(?=sp)")) %>%
  mutate(., p_Traits = str_extract(filename, "(?<=sp_).+(?=p.csv)"))
extractor_func<-function(x){
  result<-read.csv(x[1]) %>% mutate(n_Taxa = x[2], p_Traits = x[3])
  return(result)
}
compiled_results <- apply(my_results, 1, extractor_func)
compiled_results<- do.call(rbind, compiled_results) %>% as_tibble() %>% select(-X)

correct_answers<-compiled_results %>% filter(V1 == "A")

ggplot(correct_answers, aes(fill=V1, x=p_Traits))+
  #geom_bar(fill="grey30")+
  geom_bar(fill="#3C5488FF")+
  scale_y_continuous(name = "Percent Correct Model Selected", breaks = seq(0, 100, by = 10), minor_breaks = seq(0, 100, by = 5))+
  scale_x_discrete(name = "Number of Traits")+
  theme_bw()+
  theme(legend.position="none", panel.grid = element_line(size=.75), panel.grid.minor = element_line(size=.75))

ggsave("resultplot1.pdf", device = cairo_pdf,width =2, height = 3.149, units = "cm")
