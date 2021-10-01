library("tidyverse")
library("magrittr")

# Import drug list 
drugs.list <- read.delim("/DATA/usr/m.trauernicht/projects/EpiScreen/files_scripts/mt20190123_EpigeneticDrugList.txt")
drugs.list <- drugs.list %>%
  dplyr::select(-IC50..µM.) %>%
  rename(drug_well = Number) %>%
  mutate(drug_well = gsub("^([A-Z])(?<![0-9])([0-9])(?![0-9])", "\\10\\2", drug_well, perl = TRUE))

# Replicate 2 MiSeq
seq_folder = "/DATA/projects/DSBrepair/data/6421_DrugScreen"

file_list = list.files(paste0(seq_folder, "/bc_demux"), full.names = TRUE)
head(file_list)

full_meta_data_run1 = read.table("/DATA/projects/DSBrepair/git/EpiScreen/rs20210407_multiplexing_scheme_E1504.txt", 
                            header = F, sep = "\t",
                            col.names = c("no", "ID", "descr", NA, 
                                          "seq_barcode", "seq_index", NA)) %>%
  dplyr::select(no, ID, descr, seq_barcode, seq_index)  %>%
  as_tibble() %>%
  mutate(sample = substr(ID, 10, nchar(ID)),
         ID = paste0(ID, "_run1")) %>%
  distinct(ID, seq_barcode, seq_index, sample)

file_list_tib = tibble(file = file_list) %>% 
  mutate(sample = gsub(".*[CTGA]{8}_(.*).fq.gz", "\\1", file))

meta_rep2_run1_tib = left_join(full_meta_data_run1, file_list_tib) %>% 
  mutate(PCR_type = "indelPCR", index_length = "10", sample = NULL)

plate_no_tib = meta_rep2_run1_tib %>%
  mutate(plate_no = gsub("E1504_(..)_.*", "\\1", ID),
         ref =  gsub("E1504_.._(.*)_...", "\\1", ID)) %>% 
  distinct(plate_no, ref)

no_reads_rep2_tib = meta_rep2_run1_tib %>% 
  filter(is.na(file))

# Replicate 2 NextSeq
seq_folder = "/DATA/projects/DSBrepair/data/6475_DrugScreen_NextSeqMed"

file_list = list.files(seq_folder, full.names = TRUE)
head(file_list)

full_meta_data_run2 = full_meta_data_run1 %>% 
  mutate(ID = gsub("1$", "2", ID))

file_list_tib = tibble(file = file_list) %>% 
  mutate(sample = gsub(".*[CTGA]{8}_E1504_.._(.*).fq.gz", "\\1", file))

meta_rep2_run2_tib = left_join(full_meta_data_run2, file_list_tib) %>% 
  mutate(PCR_type = "indelPCR", index_length = "10", sample = NULL)

plate_no_tib = meta_rep2_run2_tib %>%
  mutate(plate_no = gsub("E1504_(..)_.*", "\\1", ID),
         ref =  gsub("E1504_.._(.*)_...", "\\1", ID)) %>% 
  distinct(plate_no, ref)

no_reads_rep2_tib = meta_rep2_run2_tib %>% 
  filter(is.na(file))



# Replicate 1 MiSeq
seq_folder = "/DATA/projects/DSBrepair/data/5280_DrugScreen"

file_list = list.files(paste0(seq_folder, "/raw"), full.names = TRUE)
head(file_list)

full_meta_data = read.table("/DATA/projects/DSBrepair/git/EpiScreen/rs20210526_EpiScreen1_GCF.txt", 
                            header = F, sep = "\t",
                            col.names = c("no", "name", "descr", NA, 
                                          "seq_barcode", "seq_index")) %>%
  dplyr::select(no, name, descr, seq_barcode, seq_index) %>%
  mutate(ref = gsub("(.*_._.)_.*$","\\1", name)) %>%
  left_join(., plate_no_tib, by = "ref") %>%
  mutate(name = gsub("(?<![0-9])([0-9])(?![0-9])$", "0\\1", name, perl = TRUE)) %>%
  mutate(ID = paste("E177", plate_no, name, "run1", sep = "_")) %>% 
  distinct(ID, seq_barcode, seq_index)

file_list_tib = tibble(file = file_list) %>% 
  mutate(ref = gsub(".*[CTGA]{8}_(.*)_R1.fq.gz", "\\1", file),
         well = gsub(".*_(.*)_[CTGA]{8}_.*", "\\1", file)) %>% 
  mutate(well = gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", well, perl = TRUE)) %>%
  left_join(plate_no_tib) %>%
  mutate(ID = paste("E177", plate_no, ref, well, "run1", sep = "_")) %>% 
  distinct(ID, file)
  

meta_rep1_tib = left_join(full_meta_data, file_list_tib) %>% 
  mutate(PCR_type = "indelPCR", index_length = "10")

no_reads_rep1_tib = meta_rep1_tib %>% 
  filter(is.na(file))

meta_tib = rbind(meta_rep2_run1_tib, meta_rep2_run2_tib, meta_rep1_tib)  %>%
  separate(ID, c("replicate", "plate", "concentration", "tech", "drug_plate", "well", "run"), "_", remove = FALSE) %>%
  mutate(drug_no = rep(sprintf('%0.3d', 1:192), times = 27),
         drug_well = paste(well, drug_no, sep = "_")) %>% 
  left_join(drugs.list)

meta_tib %<>% rename_with(tolower, -c("ID", "PCR_type"))


write.table(filter(meta_tib, !is.na(file)), file = "/DATA/projects/DSBrepair/config/rs20210526_E177_E1504_EpiScreen_metadata.txt", sep =  "\t", row.names = F, quote = F)
