library(dplyr)
library(readr)
library(tidyr)

finalnet <- read_csv("finaltraits.csv", show_col_types = FALSE)

all_stats <- finalnet %>%
  summarise(across(where(is.numeric),
                   list(
                     n    = ~sum(!is.na(.x)),
                     mean = ~mean(.x, na.rm = TRUE),
                     sd   = ~sd(.x, na.rm = TRUE),
                     se   = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
                     min  = ~min(.x, na.rm = TRUE),
                     max  = ~max(.x, na.rm = TRUE)
                   ),
                   .names = "{.col}__{.fn}"
  )) %>%
  pivot_longer(everything(),
               names_to = c("trait", "stat"),
               names_sep = "__",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  arrange(trait)

all_stats
write_csv(all_stats, "supplementary_table_trait_descriptives.csv")
