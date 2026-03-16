library(stringr)
library(sybil)
library(optparse)
library(glpkAPI)
getwd=getwd()
getwd


option_list = list(
  optparse::make_option(c("-f", "--file"),help="model name ")
  # optparse::make_option(c("-a", "--add_reaction"), help="reaction file")
)
opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);


model_nm=opt$file
model= readRDS(paste(model_nm,'.RDS',sep=''))

outfile_add_043= paste(getwd,'/043_model_',model_nm,'_model_list_react.tsv',sep = '')      # This output includes all reactions of your model.
outfile_add_045= paste(getwd,'/045_model_',model_nm,'_model_list_metabolites.tsv',sep = '')# This output includes all metabolites of your model.
outfile_add_046= paste(getwd,'/046_model_',model_nm,'_model_Smat.csv',sep = '')            # This output includes all stoichiometry of reactions.
infile_046     = outfile_add_046
outfile_050    = paste(getwd,'/050_model_',model_nm,'_model_Smat_fluxes.csv',sep = "")     # This output includes metabolic flux of each reaction.

########################## st4. Get reactions of a model (MM). #############################  
# 043.What are those reactions of the model?
#model@react_num
# Extract all reactions.
model_list<-list(model@react_id,model@react_name,model@react_rev,model@react_single, model@react_de,
                 model@lowbnd, model@uppbnd, model@obj_coef, model@gprRules, model@gpr,
                 model@react_attr)
df_model_list  <- as.data.frame(model_list, row.names = NULL)
colnames(df_model_list) <- c('react_id','react_name','react_rev','react_single','react_de','lowbnd','uppbnd','obj_coef','gprRules','gpr'
                                            ,'seed',"rxn","name","ec","tc","qseqid","pident","evalue","bitscore","qcovs","stitle","sstart"
                                            ,"send","pathway","blast.status","pathway.status","complex","exception","complex.status","gs.origin","annotation","MNX_ID","seedID","keggID","biggID","biocycID"
)


########################## st5. Get metabolites of a model (MM). #############################  
# 045.What are the metabolites of this model?

# number of reactions of ******model_Archaea2_model****: 
#model@met_num
model_met_name = model@met_name
model_met_id = model@met_id
# model_Archaea2_model: Generate a dataframe with metabolite info for each model
model_list_met <- list(model@met_id, model@met_name, model@met_comp, model@met_single, model@met_de)
df_model_list_met <- as.data.frame(model_list_met, row.names = NULL)
colnames(df_model_list_met) <- c('met_id','met_name','met_comp','met_single','met_de')
# R Sort a Data Frame using Order()
df_model_list_met_sort <- df_model_list_met[order(df_model_list_met$met_id),]


######################################################################################################
# 046.What are the stiochiometric matrix (model_xxx@S) in this model?


model_Smat <- model@S
model_Smat_mx = as.matrix(model@S)

colnames(model_Smat_mx)<-cbind(model@react_id)

row.names(model_Smat_mx)
row.names(model_Smat_mx)<-cbind(model@met_id)



write.table(df_model_list, file = outfile_add_043,quote = F,row.names = F,sep = '\t')

write.table(df_model_list_met_sort, file = outfile_add_045, quote = F,row.names = F,sep = '\t')

write.csv(model_Smat_mx, file = outfile_add_046, quote = F,sep = ',') 


write.csv(model_Smat_mx, file = outfile_add_046, quote = F,sep = ',');
  sol <- sybil::optimizeProb(model, retOptSol=F, algorithm = "fba")
  #sol<-simulation_loop[[1]]@mfluxlist[[2]][["STY_Merged_OTU04"]][1:(length(simulation_loop[[1]]@mfluxlist[[2]][["STY_Merged_OTU04"]])/3)]
  model_flux_list <- as.data.frame(sol$fluxes)
  #length(model_Archaea2_model_flux_list)
  
  
  #dim(model_Archaea2_model_flux_list)
  
  
  
  
  model_Smat <- read.csv(file = infile_046,row.names = 'X')
  t_model_Smat <- as.data.frame(t(model_Smat))
  #dim(t_model_Smat)
  
  cbind_t_model_Smat<- as.data.frame(cbind(model_flux_list, t_model_Smat))
write.csv(cbind_t_model_Smat, quote = FALSE,file = outfile_050,sep=',')
