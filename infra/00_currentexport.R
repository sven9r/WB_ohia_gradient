source(here::here("scripts/00_source.R"))

#saves current date into a variable
current_date <- format(Sys.Date(), "%Y-%m-%d")
current_df <- glue("output/{current_date}_zotu_df.csv")
current_tax <- glue("output/{current_date}_taxonomy_df.csv")
current_meta <- glue("output/{current_date}_meta_data.csv")

write.csv(zotu_df, current_zotu)
write.csv(taxonomy_df, current_tax)
write.csv(meta_data, current_meta)
