import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('-E.g',
                    required=False,
                    help=' -indir -infile')
parser.add_argument('-mydir',default='/srv/scratch/z5245780/SpongeMAGs_Robbins_ISME/all_sponge_DB_bins_over50pctComp/CAR/genome/Bacarena_community',
                    required=False,
                    help='Full path of input files.')
parser.add_argument('-infile',
                    required=True,
                    help='Name of input files.(generated from RDS by the R script')

args = vars(parser.parse_args())

mydir = args['mydir']
infile_original = args['infile']

#################### Or run script within pycharm:
######################## user define starts #########################

######################## user define ends #########################
outfile_1_cpds_sum_in_species = mydir + '/' + infile_original.split('.csv')[0] + '_sum_in_species.csv'
outfile_2_network_format      = mydir + '/' + infile_original.split('.csv')[0] + '_network.csv'

#
# st1. get dicts of dicts: for the first dict, each species are keys,
# for the second dict, the cpd IDs are the keys and fluxes are values.
# At last, by providing species, will get you  its flux of each metabolites.
dict_of_dict ={}
header_list = []
#for each_line in open(mydir +'/'+ infile_original):
for each_line in open(mydir + '/'+infile_original):
    each_line_split = each_line.strip('\n').split(',')

    if each_line.startswith('"species",') or each_line.startswith('species,') :
        header_list = each_line_split
    else:
        species_id = each_line_split[0]

        if species_id not in dict_of_dict:
            dict_of_dict[species_id] = {}

        for (cpd, flx) in zip(header_list[1:], each_line_split[1:]):
            if flx != 'NA':

                if cpd not in dict_of_dict[species_id]:
                    dict_of_dict[species_id][cpd] = float(flx)
                else:
                    dict_of_dict[species_id][cpd] += float(flx)
                    
# st2. write the table of sum of different cpds in each species to an output file.
outfile_handle = open(outfile_1_cpds_sum_in_species, 'w')
outfile_handle.write(','.join(header_list) + '\n')
for each_species in dict_of_dict:
    flx_sum_list = []
    for each_cpd in header_list[1:]:
        flx_sum_list.append(dict_of_dict[each_species].get(each_cpd, 'NA'))

    flx_sum_list_str = [str(i) for i in flx_sum_list]


    outfile_handle.write('%s,%s\n' % (each_species, ','.join(flx_sum_list_str)))
outfile_handle.close()


# st3. write a table in the form for network graph in r.
outfile_handle = open(outfile_2_network_format, 'w')
outfile_handle.write('%s,%s,%s,%s,%s\n' % ('met', 'prod', 'cons', 'prod.flux', 'con.flux'))

for each_cpd in header_list:

    provide_list = []
    comsumer_list = []
    for each_species in dict_of_dict:
        # print('%s\t%s' % (each_species, dict_of_dict[each_species].get(each_cpd, 0)))
        flx_sum = dict_of_dict[each_species].get(each_cpd, 0)
        if flx_sum > 0:
            provide_list.append(each_species)
        if flx_sum < 0:
            comsumer_list.append(each_species)

    if (provide_list == []) and (comsumer_list != []):
        provide_list = ["Environment"]
        # print('%s\t%s\t%s' % (each_cpd, provide_list, comsumer_list))
    if (provide_list != []) and (comsumer_list == []):
        comsumer_list = ["Environment"]
        # print('%s\t%s\t%s' % (each_cpd, provide_list, comsumer_list))

    if (provide_list != []) and (comsumer_list != []):
        # pass
        #print('%s\t%s\t%s' % (each_cpd, provide_list, comsumer_list))

        for each_pro in provide_list:
            for each_con in comsumer_list:

                pro_flux = dict_of_dict.get(each_pro, {}).get(each_cpd, 0) # provide blanks if keys can not be found in the dictionaries.
                con_flux = dict_of_dict.get(each_con, {}).get(each_cpd, 0)

                outfile_handle.write('%s,%s,%s,%s,%s\n' % (each_cpd, each_pro, each_con, pro_flux, con_flux))

outfile_handle.close()

###########code from here to the top can be used to check the specific cycle crossfeeding, need to add the species at the beginning of first line and modify 'prod_flux' to 'prod.flux'...#######
##########change cpd_id to met, providers to prod, consumers to cons.#########
##################### Or run script within pycharm:
######################### user define starts #########################

dict_seed_metabolites_edited_gapseq = '/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/seed_metabolites_edited.tsv'

######################### user define ends #########################
infile_2  = mydir + '/' + infile_original.split('.csv')[0] + '_network.csv'
outfile_3 = mydir + '/' + infile_original.split('.csv')[0] + '_out_final_cf.csv'

# st1. create a dict of MS metabolite id to metabolite name.
dict_metid_2_metnm = {}
for each in open(dict_seed_metabolites_edited_gapseq):
    each_split = each.strip().split('\t')
    if not each.startswith('id'):
        met_id = each_split[0]
        met_nm = each_split[3]
        dict_metid_2_metnm[met_id] = met_nm
        # print(met_id)

# st2. use met_id to get met_nm from the input file.
outfile_handle = open(outfile_3, 'w')
for each in open(infile_2):
    each_split = each.strip().split(',')
    if each.startswith('met'):
        outfile_handle.write('%s,%s,%s\n' % (each_split[0], 'cpd_nm', ','.join(each_split[1:])))
    else:
        met_id = each_split[0].split('EX_')[1].split('_')[0]
        # print(met_id)
        met_nm = dict_metid_2_metnm[met_id]
        outfile_handle.write('%s,%s,%s\n' % (met_id, met_nm, ','.join(each_split[1:])))
outfile_handle.close()
os.system('rm '+infile_2)
