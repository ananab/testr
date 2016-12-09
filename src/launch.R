

args <- commandArgs(trailingOnly=TRUE)

lkup_path <- args[1]
ref_path <- args[2]

rmarkdown::render("src/testdata.Rmd", output_dir = "results")
