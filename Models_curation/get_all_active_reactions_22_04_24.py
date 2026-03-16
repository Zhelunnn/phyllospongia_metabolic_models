import argparse
import re
def get_all(file_args):
    DB=open(file_args['db'])#######database file seed_reaction_corrected.tbl
    db=DB.readlines()


    input=file_args['i']
    output='original_active.txt'
    file = open(input,'r')
    new_file = open(output, 'w')
    file_list=file.readlines() ### 050file a matrix in form rxnid as row name and cpdid as column name each grid is coef, except the second column is flux
    header_list=file_list[0].strip('\n').split(',')##use the first line to correspond to the compound id
    indx = 0
    for line in file_list[1:]:#############read reactions line remove the first line
        element_list=line.strip('\n').split(',')
        if float(element_list[1])!=0:###########operate on fluxes!=0 line
           reaction=''##########extract compound in this rxn
           for index,i in enumerate(element_list[2:]):##############add coef before the compound id
               if float(i)!=0:
                   reaction+=i+header_list[index+2]+','

           IDX = element_list[0][0:8]######rxn_id in fluxes file
           if IDX.startswith('rxn'):
               for ref_line in db:###############extract the annotation of rxn from seed database
                   if IDX in ref_line:
                       annota=ref_line.strip('\n').split('\t')[2]
                       new_file.write('%s\t%s\t%s\t%s\n' %(element_list[0],annota,reaction.strip(","),element_list[1]))
                       break
           else:
               new_file.write('%s\t%s\n' %(element_list[0],element_list[1]))
    new_file.close()
    file.close()
###########################get the rxn equation from the active file##################
def get_rxn_eqa(active_file,metabolite_db):
    result=[]
    for line in open(active_file):
        if line.startswith('rxn'):
            rxn_eqa=line.strip('\n').split('\t')[2]
            flux=line.strip('\n').split('\t')[3]
            compound=rxn_eqa.split(',')
            consume=[]######consumption compound
            product=[]######product
            for cpd in compound:############allocate the reacting cpd and product to two cpd list
                cpd_id = re.search(f"{'cpd'}(.{{{5}}})", cpd).group(0) ### get cpd id by search for 'cpd' and followed 5 id number
                name=find_cpd_name(cpd_id,metabolite_db)
                if name!=None:
                    cpd_id_name = cpd.strip('-') + '<' + name + '>'  # write cpd id and name as 'coef+cpd_id+<cpd_name>'
                else:
                    cpd_id_name = cpd.strip('-') + '<' + cpd_id + '>'  # write cpd id and name as 'coef+cpd_id+<cpd_name>'

                if float(flux)<0:
                    if cpd.startswith('-'):
                        product.append(cpd_id_name)
                    else:
                        consume.append(cpd_id_name)
                else:
                    if cpd.startswith('-'):
                        consume.append(cpd_id_name)
                    else:
                        product.append(cpd_id_name)
################write the cpd from two cpd list into a equation
            consu='+'.join(consume)
            prod='+'.join(product)
            equation=consu+'   >>   '+prod
            new_line=line.strip('\n').split('\t')[0]+'\t'+line.strip('\n').split('\t')[1]+'\t'+equation.strip('+')+'\t'+flux.strip('-')+'\n'
            result.append(new_line)
    return result
def find_cpd_name(cpd_id,metabolite_db): ### input cpd_id, return the correspongding cpd name from the metabolites db file
    for line in metabolite_db:
        line_id=line.split('\t')[0]
        if cpd_id==line_id:
            cpd_name=line.split('\t')[3]
            if cpd_name==None:
                return cpd_id
            else:
                return cpd_name

if __name__=='__main__':
    active_reactions=argparse.ArgumentParser()
    active_reactions.add_argument('-i',required=True,help='input file should be 050 file')
    active_reactions.add_argument('-o',required=True,help='output file was required and the format is rxn,coeffi+cpdid,fluxes')
    active_reactions.add_argument('-db',default='/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/seed_reactions_corrected.tsv', help='the database for reactions')
    active_reactions.add_argument('-db_m', default='/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/seed_metabolites_edited.tsv',help='the metabolite db file ')
    args = vars(active_reactions.parse_args())
    get_all(args)
    output_file=open(args['o'],'w')
    metabolite_db=open(args['db_m']).readlines()
    for line in get_rxn_eqa('original_active.txt',metabolite_db):
        output_file.write(line)
    output_file.close()