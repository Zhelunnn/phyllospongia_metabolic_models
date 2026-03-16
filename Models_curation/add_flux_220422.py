import argparse
import sys
def add_fluxes(rxn_file,flux_file,output_fluxes):
    new_file=open(output_fluxes,'w')
    flux_file_open=open(flux_file).readlines()
    for line in open(rxn_file):
        rxn_id=line.strip('\n').split('\t')[0][0:8]
        for flux_line in flux_file_open:####open(flux_file):
            if rxn_id in flux_line:
                flux=flux_line.strip('\n').split(',')[1]
                new_line=line.strip('\n')+'\t'+flux+'\n'
                #print(new_line)
                new_file.write(new_line)
                #flux_file.close()
                break
    new_file.close()
##########add information of rxn##########
def add_info(db,rxn_file,output_info):
    new_file = open(output_info, 'w')
    db_file_open = open(db).readlines()
    for line in open(rxn_file):
        if line.startswith('rxn'):
            for db in db_file_open:
                rxn_id=line.strip('\n').split('\t')[0]
                db_id=db.strip('\n').split('\t')[0]
                if rxn_id in db_id:
                    info=db.strip('\n').split('\t')[2]+'\t'+db.strip('\n').split('\t')[4]
                    new_line=rxn_id+'\t'+info+'\n'
                    new_file.write(new_line)
    new_file.close()
#add_fluxes('final.txt','050_model_cluster_CAR1_bin_2_model_Smat_fluxes.csv')
if __name__=='__main__':
    add_flux=argparse.ArgumentParser()
    add_flux.add_argument('-i',required=True,help='input_file')
    add_flux.add_argument('-F', action='store_true',help='add the flux')
    add_flux.add_argument('-I',action='store_true',help='add the info')
    add_flux.add_argument('-f',default='050_model_cluster_CAR1_bin_16_model_Smat_fluxes.csv',help='050 flux file/no pwy found file')
    add_flux.add_argument('-o',default='fluxes_added.txt',help='flux_output file/info output file')
    add_flux.add_argument('-d',default='/srv/scratch/z5245780/software/gapseq/gapseq/dat/seed_reactions_corrected.tsv',help='database file')
    args=vars(add_flux.parse_args())
    if '-F' in sys.argv:
        add_fluxes(args['i'],args['f'],args['o'])
    if '-I' in sys.argv:
        add_info(args['d'],args['i'],args['o'])