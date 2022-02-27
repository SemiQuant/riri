# require(plotly)
# require(tidyverse)
# 
# files <- list.files("/Users/SemiQuant/Downloads/riboDelete/out_dir_test", full.names = T)
# 
# 
# data <- files %>%
#   map(read_tsv, col_names = F, show_col_types = FALSE) %>%
#   reduce(cbind) 
# colnames(data) <- gsub(".txt", "", basename(files))
# 
# data %>% 
#   mutate(n = 1:nrow(data)) %>% 
#   pivot_longer(!n, names_to = "Gene", values_to = "count") %>% 
#   plot_ly(x = ~n, y = ~count, type = 'scatter', mode = 'lines', color = ~Gene)
# 
# 
# 
# t() %>% 
#   data.frame() %>% 
#   
#   


