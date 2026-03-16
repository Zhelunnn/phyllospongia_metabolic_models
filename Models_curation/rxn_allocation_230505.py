############################allocate the reactions to pwy based on different database. file=043######################
#######################database file was accepted ################
#####################pwys that found multiple rxns,single rxn and no rxn were added into -o,-s and -n #################
import re
import copy
import argparse
import sys
# def get_id_db(db,file):####correspond rxn_id with reaction_id in other database (search in the website column)
#     if db=='META':
#         label='META'
#     else:
#         label=db+'.reaction'
#     id_list=[]
#     for line in open(file):
#         if line.startswith('rxn'):
#             rxn_id=line.strip('\n').split('\t')[0][0:8]
#             result=re.search(label+':(.*?);',line)#############get the rxn_id in given database
#             if result:
#                 #print(result)
#                 db_id=result.group(1)
#                 #print(db_id)
#                 id=rxn_id+'_'+db_id
#                 id_list.append(id)
#             else:
#                 id=rxn_id+'_NA'
#                 id_list.append(id)
#     return id_list
#########search based on the database unique column in 043file
def correspond_rxnid(db,file):####db:db_type,file:043file
    if db=='META':
        ordinal=35
    elif db=='kegg':
        ordinal=33
    elif db=='metanetx':
        ordinal=31
    elif db=='bigg':########need to improve
        ordinal=34
    else:
        print('database type was not included')
        ordinal =0
    result_list=[]
    for line in open(file):
        if line.startswith('rxn'):
            rxn_id=line.strip('\n').split('\t')[0][0:8]
            correspond_id=line.strip('\n').split('\t')[ordinal]
            if db=='META':###############search the 'biocycID' column first then try to get it 'rxn' column
                if correspond_id=='NA':
                    correspond_id=line.strip('\n').split('\t')[11]
                    if line.strip('\n').split('\t')[12]=='NA':
                        correspond_id=rxn_id
            new_line=rxn_id+'_'+correspond_id
            result_list.append(new_line)
    return result_list

#x=get_id_db('META','043_model_cluster_CAR1_bin_2-draft_model_list_react.txt')
def get_info_pwy(pwy_line):########extract pwy information based on each pwy line from modelseed
    pwy_id=pwy_line.strip().split('\t')[0]######RXN ID
    #pwy_name=pwy_line.strip().split('\t')[1]+'\t'+pwy_line.strip().split('\t')[3]########RXN NAME
    pwy_name=pwy_line.strip().split('\t')[1]
    pwy_rxn=pwy_line.strip().split('\t')[5]
    #pwy_info=pwy_id+'_'+pwy_name+'_'+pwy_rxn
    pwy_info = pwy_name + '\t' + pwy_rxn
    return pwy_info
#########get the missed rxn in the pwy###########
def add_rxn_absen(key,pwyinfo_dict,rxn_input): 
    db_rxn=pwyinfo_dict[key].split('_')[-1] ###get all reactions belonging to this pwy
    
    rxn_list=db_rxn.split('\t')[1].split(',')### remove the pwy name
    result=''
    for i in rxn_list:
        if i not in rxn_input:
            result+=i+','
    return result.strip(',')
###############################allocate reaction into pwyid ########################
def pwy_rxn(datatype,reaction_list):#######reaction_list should be the result from get_id_db() function
    if datatype=='META':
        database ='/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/meta_pwy.tbl'
    elif datatype=='kegg':
        database = '/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/kegg_pwy.tbl'
    elif datatype=='biocyc':
        database = '/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/meta_pwy.tbl'
    # elif datatype=='metanetx':
    #     database =
    else:
        database='/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/seed_pwy.tbl'
        #database = 'seed_pwy.tbl'

    DB=open(database)
    db=DB.readlines()
    db_pwy={}########pwy:rxn
    rxn_no_pwy=[]#####reactions that are not found in pwy
    pwy_library={}########save the information of found pwy
    rxn_NA=[]
    for line in reaction_list:##############for loop in fill reaction file
        if line.startswith('rxn'):##########if this line refer to reaction
            element_list=line.strip().split('_')
            rxn_id=element_list[1]
            if rxn_id!='NA':
                dp_db = copy.deepcopy(db_pwy)##########deepcopy used to supervise the change of db_pwy
                for db_line in db:###########for loop search in databse
                    if rxn_id in db_line:
                       db_id=db_line.strip().split('\t')[0]
                       #print(db_id)
                       if db_id in list(db_pwy.keys()):
                           db_pwy[db_id]+=','+rxn_id
                       else:
                           db_pwy[db_id]=rxn_id
                           pwy_library[db_id]=get_info_pwy(db_line)#########save pwy information to a new dict
                if db_pwy==dp_db:
                    rxn_no_pwy.append(rxn_id)
            else:
                rxn_NA.append(element_list[0])
    return db_pwy,pwy_library,rxn_no_pwy,rxn_NA
def divide_pwy(db_pwy,pwy_library):
    rxn_single_pwy={}
    multi_pwy=[]
    for i in db_pwy.keys():
        if len(db_pwy[i].split(','))>1:
            rxn_col=''
            for rxn_str in db_pwy[i].split(','):
                rxn_col+=rxn_str+','
                #print(rxn_col)
            rest_rxn=add_rxn_absen(i,pwy_library,rxn_col)
            final=i+'\t'+pwy_library[i]+'\t'+rxn_col.strip(',')+'\t'+rest_rxn+'\n' #### second last column: found filled reactions. last column: the rest reactions of the pwy
            multi_pwy.append(final)
        else:
            rxn_sin=db_pwy[i]
            rxn_single_pwy[i]=rxn_sin
    return multi_pwy,rxn_single_pwy
# print(pwy_rxn('meta_pwy.tbl',x)[0])
# print(divide_pwy(pwy_rxn('meta_pwy.tbl',x)[0],pwy_rxn('meta_pwy.tbl',x)[1])[0])
if __name__=='__main__':
    allocate_pwy=argparse.ArgumentParser()
    allocate_pwy.add_argument('-dt',help='META,kegg, seed is default database, keep dt parameter empty if you want to use modelseed')
    allocate_pwy.add_argument('-f',required=True,help='reactions file 043')
    allocate_pwy.add_argument('-o',required=False,default='pwy_filled_rxn.txt',help='for pathway found with multi reactions')
    allocate_pwy.add_argument('-n',required=False,help='no pwy found reaactions')
    allocate_pwy.add_argument('-s', required=False, help=' pwy found with single reaaction')
    args=vars(allocate_pwy.parse_args())
    #result=pwy_rxn(args['dt'],get_id_db(args['dt'],args['f']))
    result=pwy_rxn(args['dt'],correspond_rxnid(args['dt'],args['f']))
    if args['o']!=None:
        file_output=open(args['o'],'w')
        file_output.write('%s\t%s\t%s\t%s\t%s\n' % ('pwy_id','pwy_name','pwy_all_reactions','filled_reactions','pwy_rest_reactions'))
        for line in divide_pwy(result[0],result[1])[0]:
            file_output.write(line)
        file_output.close()
    if args['n']!=None:
        file_output=open(args['n'],'w')
        for line in result[2]:
            file_output.write(line+'\n')
        file_output.close()
        file_NA=open('rxn_NA.txt','w')
        for NA in result[3]:
            file_NA.write(NA + '\n')
        file_NA.close()
    if args['s']!=None:
        file_output=open(args['s'],'w')
        single_pwy=divide_pwy(result[0], result[1])[1]
        for line in single_pwy.keys():
            file_output.write(line+'\t'+single_pwy[line]+'\n')
        file_output.close()

