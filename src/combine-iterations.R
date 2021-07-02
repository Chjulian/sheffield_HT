library(dplyr)

file.dir <- "./output/"

rds.files <- list.files(file.dir)

##Get file names
summary.files <- grep("summary", rds.files, value=T)
res.files <- grep("run_", rds.files, value=T)
step.files <- grep("steps.rds", rds.files, value=T)

##combine outcomes
summary.df <- lapply(summary.files, function(i){
        fun.df <- readRDS(paste0(file.dir, i))
        split.run <- strsplit(i, "_")
        run.name <- paste(split.run[[1]][3],
                          split.run[[1]][4], sep = "_")
        fun.df$run <- run.name
        print(paste0("Summary: " , run.name))
        return(fun.df)
}) %>% bind_rows()

saveRDS(summary.df, "./output/outcomes.rds")
rm(summary.df)

##combine outbreaker res objects
res.all <- readRDS(paste0(file.dir, res.files[1]))
split.run <- strsplit(res.files[1], "_")
run.name <- paste(split.run[[1]][3],
                 split.run[[1]][4], 
                 split.run[[1]][5],
                 split.run[[1]][6],
                 sep = "_")
res.all$run <- run.name

for (i in seq(2,length(res.files))){
  loop.df <- readRDS(paste0(file.dir, res.files[i]))
  split.run <- strsplit(res.files[i], "_")
  
  run.name <- paste(split.run[[1]][3],
                    split.run[[1]][4], 
                    split.run[[1]][5],
                    split.run[[1]][6],
                    sep = "_")
  loop.df$run <- run.name
  res.all <- rbind(res.all, loop.df)
  print(paste0("Res: " , i))
}

saveRDS(res.all, "./output/res-combined.rds")
rm(res.all)
      
##combine
steps.all <- readRDS(paste0(file.dir, step.files[1]))
split.run <- strsplit(step.files[1], "_")
run.name <- paste(split.run[[1]][3],
                  split.run[[1]][4], 
                  split.run[[1]][5],
                  sep = "_")
steps.all$run <- run.name

for (i in seq(2,length(step.files))){
  loop.df <- readRDS(paste0(file.dir, step.files[i]))
  split.run <- strsplit(step.files[i], "_")
  
  run.name <- paste(split.run[[1]][3],
                    split.run[[1]][4], 
                    split.run[[1]][5],
                    sep = "_")
  loop.df$run <- run.name
  steps.all <- rbind(steps.all, loop.df)
  print(paste0("Step: " , i))
}

saveRDS(steps.all, "./output/steps-combined.rds")
