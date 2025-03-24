
#!/usr/bin/env rscript
#2024 07 20
#this script takes in the file output from samtools coverage which calculated overall depth, and using the largest scaffolds (super scaffolds) calculates the proportion of reads
#to downsample by to take overall coverage to 5x



#read in command line argvs
inputs <- commandArgs(TRUE)
data_path<- inputs[1] #bed file with conservation scores
#workdir <- inputs[2] #working directory

library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)

# prepare files for downsampling high coverage individuals
options(scipen = 999)

# read in output of samtools coverage
#data_path <- "/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/markdups_bam_91123
files <- dir(data_path, pattern = "*_cov.gz")  # adjust pattern to only include the desired file type
print(files)



cov <- tibble(filename = files) %>%
  mutate(file_contents = map(filename, ~ fread(file.path(data_path, .), fill = TRUE))) %>%
  unnest(cols = c(file_contents))

# extract chr length scaffolds
scaffolds <- c(
  "Super-Scaffold_1", 
  "Super-Scaffold_102", 
  "Super-Scaffold_103", 
  "Super-Scaffold_104", 
  "Super-Scaffold_105", 
  "Super-Scaffold_107", 
  "Super-Scaffold_10", 
  "Super-Scaffold_111", 
  "Super-Scaffold_112", 
  "Super-Scaffold_114", 
  "Super-Scaffold_115", 
  "Super-Scaffold_118", 
  "Super-Scaffold_11", 
  "Super-Scaffold_120", 
  "Super-Scaffold_121", 
  "Super-Scaffold_122", 
  "Super-Scaffold_123", 
  "Super-Scaffold_124", 
  "Super-Scaffold_126", 
  "Super-Scaffold_128", 
  "Super-Scaffold_129", 
  "Super-Scaffold_12", 
  "Super-Scaffold_130", 
  "Super-Scaffold_132", 
  "Super-Scaffold_133", 
  "Super-Scaffold_134", 
  "Super-Scaffold_135", 
  "Super-Scaffold_136", 
  "Super-Scaffold_139", 
  "Super-Scaffold_13", 
  "Super-Scaffold_141", 
  "Super-Scaffold_146", 
  "Super-Scaffold_147", 
  "Super-Scaffold_148", 
  "Super-Scaffold_149", 
  "Super-Scaffold_14", 
  "Super-Scaffold_150", 
  "Super-Scaffold_151", 
  "Super-Scaffold_153", 
  "Super-Scaffold_154", 
  "Super-Scaffold_156", 
  "Super-Scaffold_158", 
  "Super-Scaffold_15", 
  "Super-Scaffold_160", 
  "Super-Scaffold_162", 
  "Super-Scaffold_163", 
  "Super-Scaffold_165", 
  "Super-Scaffold_166", 
  "Super-Scaffold_167", 
  "Super-Scaffold_168", 
  "Super-Scaffold_16", 
  "Super-Scaffold_170", 
  "Super-Scaffold_171", 
  "Super-Scaffold_172", 
  "Super-Scaffold_176", 
  "Super-Scaffold_179", 
  "Super-Scaffold_17", 
  "Super-Scaffold_181", 
  "Super-Scaffold_182", 
  "Super-Scaffold_186", 
  "Super-Scaffold_187", 
  "Super-Scaffold_188", 
  "Super-Scaffold_189", 
  "Super-Scaffold_18", 
  "Super-Scaffold_190", 
  "Super-Scaffold_194", 
  "Super-Scaffold_195", 
  "Super-Scaffold_196", 
  "Super-Scaffold_197", 
  "Super-Scaffold_199", 
  "Super-Scaffold_19", 
  "Super-Scaffold_201", 
  "Super-Scaffold_202", 
  "Super-Scaffold_203", 
  "Super-Scaffold_204", 
  "Super-Scaffold_206", 
  "Super-Scaffold_209", 
  "Super-Scaffold_20", 
  "Super-Scaffold_210", 
  "Super-Scaffold_211", 
  "Super-Scaffold_216", 
  "Super-Scaffold_217", 
  "Super-Scaffold_21", 
  "Super-Scaffold_221", 
  "Super-Scaffold_222", 
  "Super-Scaffold_225", 
  "Super-Scaffold_227", 
  "Super-Scaffold_228", 
  "Super-Scaffold_229", 
  "Super-Scaffold_22", 
  "Super-Scaffold_231", 
  "Super-Scaffold_234", 
  "Super-Scaffold_237", 
  "Super-Scaffold_238", 
  "Super-Scaffold_239", 
  "Super-Scaffold_23", 
  "Super-Scaffold_240", 
  "Super-Scaffold_241", 
  "Super-Scaffold_246", 
  "Super-Scaffold_247", 
  "Super-Scaffold_249", 
  "Super-Scaffold_24", 
  "Super-Scaffold_251", 
  "Super-Scaffold_252", 
  "Super-Scaffold_257", 
  "Super-Scaffold_25", 
  "Super-Scaffold_260", 
  "Super-Scaffold_261", 
  "Super-Scaffold_269", 
  "Super-Scaffold_26", 
  "Super-Scaffold_271", 
  "Super-Scaffold_272", 
  "Super-Scaffold_275", 
  "Super-Scaffold_276", 
  "Super-Scaffold_277", 
  "Super-Scaffold_279", 
  "Super-Scaffold_27", 
  "Super-Scaffold_281", 
  "Super-Scaffold_282", 
  "Super-Scaffold_284", 
  "Super-Scaffold_285", 
  "Super-Scaffold_286", 
  "Super-Scaffold_288", 
  "Super-Scaffold_289", 
  "Super-Scaffold_28", 
  "Super-Scaffold_291", 
  "Super-Scaffold_292", 
  "Super-Scaffold_293", 
  "Super-Scaffold_296", 
  "Super-Scaffold_297", 
  "Super-Scaffold_29", 
  "Super-Scaffold_2", 
  "Super-Scaffold_302", 
  "Super-Scaffold_303", 
  "Super-Scaffold_307", 
  "Super-Scaffold_30", 
  "Super-Scaffold_311", 
  "Super-Scaffold_312", 
  "Super-Scaffold_314", 
  "Super-Scaffold_318", 
  "Super-Scaffold_31", 
  "Super-Scaffold_323", 
  "Super-Scaffold_325", 
  "Super-Scaffold_326", 
  "Super-Scaffold_327", 
  "Super-Scaffold_328", 
  "Super-Scaffold_32", 
  "Super-Scaffold_330", 
  "Super-Scaffold_331", 
  "Super-Scaffold_333", 
  "Super-Scaffold_335", 
  "Super-Scaffold_336", 
  "Super-Scaffold_337", 
  "Super-Scaffold_33", 
  "Super-Scaffold_340", 
  "Super-Scaffold_341", 
  "Super-Scaffold_342", 
  "Super-Scaffold_344", 
  "Super-Scaffold_345", 
  "Super-Scaffold_348", 
  "Super-Scaffold_349", 
  "Super-Scaffold_352", 
  "Super-Scaffold_354", 
  "Super-Scaffold_355", 
  "Super-Scaffold_356", 
  "Super-Scaffold_358", 
  "Super-Scaffold_35", 
  "Super-Scaffold_362", 
  "Super-Scaffold_363", 
  "Super-Scaffold_364", 
  "Super-Scaffold_366", 
  "Super-Scaffold_367", 
  "Super-Scaffold_368", 
  "Super-Scaffold_369", 
  "Super-Scaffold_36", 
  "Super-Scaffold_370", 
  "Super-Scaffold_371", 
  "Super-Scaffold_372", 
  "Super-Scaffold_373", 
  "Super-Scaffold_374", 
  "Super-Scaffold_375", 
  "Super-Scaffold_376", 
  "Super-Scaffold_377", 
  "Super-Scaffold_378", 
  "Super-Scaffold_379", 
  "Super-Scaffold_37", 
  "Super-Scaffold_380", 
  "Super-Scaffold_382", 
  "Super-Scaffold_384", 
  "Super-Scaffold_385", 
  "Super-Scaffold_386", 
  "Super-Scaffold_387", 
  "Super-Scaffold_388", 
  "Super-Scaffold_389", 
  "Super-Scaffold_38", 
  "Super-Scaffold_390", 
  "Super-Scaffold_391", 
  "Super-Scaffold_393", 
  "Super-Scaffold_395", 
  "Super-Scaffold_397", 
  "Super-Scaffold_39", 
  "Super-Scaffold_3", 
  "Super-Scaffold_400", 
  "Super-Scaffold_401", 
  "Super-Scaffold_403", 
  "Super-Scaffold_405", 
  "Super-Scaffold_407", 
  "Super-Scaffold_408", 
  "Super-Scaffold_40", 
  "Super-Scaffold_412", 
  "Super-Scaffold_413", 
  "Super-Scaffold_416", 
  "Super-Scaffold_418", 
  "Super-Scaffold_41", 
  "Super-Scaffold_420", 
  "Super-Scaffold_427", 
  "Super-Scaffold_429", 
  "Super-Scaffold_42", 
  "Super-Scaffold_433", 
  "Super-Scaffold_434", 
  "Super-Scaffold_435", 
  "Super-Scaffold_43", 
  "Super-Scaffold_440", 
  "Super-Scaffold_445", 
  "Super-Scaffold_446", 
  "Super-Scaffold_447", 
  "Super-Scaffold_448", 
  "Super-Scaffold_449", 
  "Super-Scaffold_44", 
  "Super-Scaffold_450", 
  "Super-Scaffold_451", 
  "Super-Scaffold_452", 
  "Super-Scaffold_453", 
  "Super-Scaffold_458", 
  "Super-Scaffold_461", 
  "Super-Scaffold_462", 
  "Super-Scaffold_465", 
  "Super-Scaffold_467", 
  "Super-Scaffold_46", 
  "Super-Scaffold_470", 
  "Super-Scaffold_471", 
  "Super-Scaffold_472", 
  "Super-Scaffold_476", 
  "Super-Scaffold_477", 
  "Super-Scaffold_481", 
  "Super-Scaffold_487", 
  "Super-Scaffold_488", 
  "Super-Scaffold_489", 
  "Super-Scaffold_48", 
  "Super-Scaffold_498", 
  "Super-Scaffold_4", 
  "Super-Scaffold_502", 
  "Super-Scaffold_509", 
  "Super-Scaffold_51", 
  "Super-Scaffold_52", 
  "Super-Scaffold_55", 
  "Super-Scaffold_56", 
  "Super-Scaffold_57", 
  "Super-Scaffold_58", 
  "Super-Scaffold_59", 
  "Super-Scaffold_5", 
  "Super-Scaffold_60", 
  "Super-Scaffold_61", 
  "Super-Scaffold_62", 
  "Super-Scaffold_63", 
  "Super-Scaffold_65", 
  "Super-Scaffold_66", 
  "Super-Scaffold_67", 
  "Super-Scaffold_68", 
  "Super-Scaffold_69", 
  "Super-Scaffold_6", 
  "Super-Scaffold_70", 
  "Super-Scaffold_71", 
  "Super-Scaffold_72", 
  "Super-Scaffold_73", 
  "Super-Scaffold_74", 
  "Super-Scaffold_76", 
  "Super-Scaffold_77", 
  "Super-Scaffold_78", 
  "Super-Scaffold_79", 
  "Super-Scaffold_7", 
  "Super-Scaffold_80", 
  "Super-Scaffold_81", 
  "Super-Scaffold_82", 
  "Super-Scaffold_83", 
  "Super-Scaffold_84", 
  "Super-Scaffold_85", 
  "Super-Scaffold_86", 
  "Super-Scaffold_87", 
  "Super-Scaffold_88", 
  "Super-Scaffold_8", 
  "Super-Scaffold_91", 
  "Super-Scaffold_92", 
  "Super-Scaffold_95", 
  "Super-Scaffold_96", 
  "Super-Scaffold_98", 
  "Super-Scaffold_99", 
  "Super-Scaffold_9"
)


cov <- cov %>%
  filter(`#rname` %in% scaffolds) %>%
  mutate(bam = "merged")

# Visualize depth and coverage
hist(cov$meandepth) # depth
hist(cov$coverage)  # breadth

# read metadata
meta <- read.csv("/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/markdups_bam_91123/")

cov <- cov %>%
  mutate(filename = gsub("_cov.gz", "", filename)) %>%
  right_join(meta, by = c("filename" = "Sample_Id"))

# Mean depth per sample
cov %>%
  dplyr::group_by(filename, WGS) %>%
  dplyr::summarise(depth = mean(meandepth),
                   min = min(meandepth),
                   max = max(meandepth))

# Mean depth across WGS run
cov %>%
  dplyr::group_by(filename, WGS, bam) %>%
  dplyr::summarise(mean = mean(meandepth),
                   min = min(meandepth),
                   max = max(meandepth)) %>%
  ungroup() %>%
  dplyr::group_by(WGS, bam) %>%
  dplyr::summarise(meandepth = mean(mean),
                   min = min(min),
                   max = max(max))

# calculate proportion to subsample
subsample <- cov %>%
  dplyr::group_by(filename, WGS) %>%
  dplyr::summarise(depth = mean(meandepth),
                   min = min(meandepth),
                   max = max(meandepth)) %>%
  filter(WGS == 1) %>%
  mutate(prop_sub = 5 / depth) %>%
  mutate(prop_sub = round(prop_sub, 2)) %>%
  mutate(prop_sub = gsub("0.", "1.", prop_sub))  # adjust for samtools view -s syntax

write.table(subsample[c(1, 6)], "/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/markdups_bam_91123/proportion_to_downsample.txt",
            quote = f, col.names = f, row.names = f)
