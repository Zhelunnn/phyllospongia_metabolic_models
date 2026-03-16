library(Matrix)
library(lattice)
library(sybil)
library(tools)
library(optparse)
library(stringr)
option_list = list(
  optparse::make_option(c("-m", "--model"),help="model name "),
  optparse::make_option(c("-a", "--add_reaction"), help="reaction file")
)
opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);
React_id=read.delim(opt$add_reaction)

model_modified=readRDS(opt$model)

# The following commond should equal to:
# model_name <- addReact(
#   model_name,
#   id='Trans_NH4',
#   met=c('cpd00013[c0]','cpd00067[e0]','cpd00013[e0]'),
#   metName = c('NH3','H+[e]','NH4'),
#   Scoef=c(1,1,-1),
#   metComp = c(1,1,1,1)
#   reactName = 'Trans_NH4',
#   reversible= FALSE,
#   lb = 0,
#   up = 1000
# )


for (i in 1:nrow(React_id)) {
  model_modified <- addReact(
    model_modified,
    id         = as.character(React_id[i,]$Rid),
    met        = str_split(as.character(React_id[i,]$met_id),',')[[1]],
    Scoef      = as.integer(str_split(as.character(React_id[i,]$met_Scoef),',')[[1]]),
    reactName  = as.character(React_id[i,]$react_name),
    reversible = React_id[i,]$react_rev,
    metName    = str_split(as.character(React_id[i,]$met_name),',')[[1]],
    metComp    = as.integer(str_split(as.character(React_id[i,]$met_comp),',')[[1]]),
    lb         = as.integer(React_id[i,]$lowbnd),
    ub         = as.integer(React_id[i,]$uppbnd)
  )
  print(i)
  
}
#########generate the file name#########
first_part <-str_split(file_path_sans_ext(basename(opt$model)),'-')[[1]][1]####get the base name without suffix
second_part<-str_split(file_path_sans_ext(basename(opt$add_reaction)),'_',n=3)[[1]][3]####seperate the curated file name by '_' to 3 parts
final_name<-paste(first_part,'_draft_adapt_',second_part,'.RDS',sep='')


saveRDS(model_modified,file=final_name)
if ('cellwall_deficiency' %in% list.files(path = ".")){
  system(paste('mv ',final_name,' cellwall_deficiency/draft'))
}else{
  system(paste('mv ',final_name,' draft'))
}

# check the latest reaction added:
#model_new@lowbnd[model_new@react_num]
# [1] "R00148_Shan_OTU08"
#model_new@react_id[model_new@react_num]
