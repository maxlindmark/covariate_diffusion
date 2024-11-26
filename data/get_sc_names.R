library(stringr)
library(taxize)
library(readr)

root_dir <- here::here(".")

# Read species names files
names <- readxl::read_xlsx(paste0(root_dir, "/data/EBS_species.xlsx")) |>
  mutate(`Common name` = as.factor(`Common name`),
         CODE = as.factor(CODE)) |>
  # Fix typo
  mutate(CODE = fct_recode(CODE, "a_yfs" = "a_fys")) |>
  # Change some common names that do not get picked up by taxize
  mutate(`Common name` = fct_recode(`Common name`,
    "pacific halibut" = "halibut",
    "Adult yellowfin sole" = "Adult yellowfish sole"
  )) |>
  mutate(
    CODE = tolower(CODE),
    `Common name` = tolower(`Common name`)
  ) |>
  separate(CODE, sep = "_", into = c("life_stage", "species"), remove = FALSE) |>
  mutate(
    life_stage = ifelse(!life_stage %in% c("a", "ej", "j"), NA, life_stage),
    common_name2 = trimws(str_remove_all(`Common name`, "adult|juvenile|early juvenile|\\?")),
    sc_name = as.character(comm2sci(common_name2))
  ) |>
  # Some common names do not have a match, check these manually
  mutate(
    sc_name = ifelse(common_name2 == "greenland turbot", "Reinhardtius hippoglossoides", sc_name),
    sc_name = ifelse(common_name2 == "bering skate", "Bathyraja interrupta", sc_name),
    sc_name = ifelse(common_name2 == "bristol bay red king crab", "Paralithodes camtschaticus", sc_name),
    sc_name = ifelse(common_name2 == "bering sea flounder", "Hippoglossoides robustus", sc_name),
    sc_name = ifelse(common_name2 == "great skulpin", "Myoxocephalus polyacanthocephalus", sc_name),
    # https://www.fisheries.noaa.gov/science-blog/dutch-harbor-snow-and-tanner-crab-growth-study-post-1
    sc_name = ifelse(common_name2 == "tanner crab", "Chionoecetes bairdi", sc_name),
    sc_name = ifelse(common_name2 == "alaska plaice", "Pleuronectes quadrituberculatus", sc_name)
  ) |>
  separate(sc_name, into = c("family", "sp"), sep = " ", remove = FALSE) |>
  mutate(
    abb_name = substring(family, 1, 1),
    abb_name = paste0("<i>", paste(paste0(abb_name, "."), sp), "</i>"),
    abb_name = ifelse(!is.na(life_stage), paste0(abb_name, " (", life_stage, ")"), abb_name)
  ) |>
  # Distinguish the red king crabs
  mutate(abb_name = ifelse(CODE == "bking", paste0(abb_name, " (Bb)"), abb_name)) |>
  rename("X" = "CODE")

write_csv(names, paste0(root_dir, "/data/clean_EBS_species.csv"))

