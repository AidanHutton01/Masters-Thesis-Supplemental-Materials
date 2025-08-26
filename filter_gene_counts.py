import pandas as pd
import argparse

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", type=str,
                        help="Enter the name or directory of the gene count file you want to split. File must be in a specific format; "
                             "each column must have a header, the header for the genes must be 'Geneid' (case sensitive), "
                             "the other headers should represent each sample")
    parser.add_argument("config", type=str,
                        help="Enter the name or directory of the config file containing the different sample names in each group. "
                             "Sample names should match the headers in the main data file.")
    parser.add_argument("output", type=str,
                        help="Enter the name or directory of the folder you want the data to be output into.")
    return parser.parse_args()


#Make pd dataframe
def read_config(fname):
    #reads each line of config file and prepares it for sorting
    with open(fname,'r') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]

    #opens list for each experimental group that we want to look at
    #(only works for space + ground groups with 2 other factors)
    space_group1 = []
    space_group2 = []
    ground_group1 = []
    ground_group2 = []

    group_flag = 0
    for line in lines:
        if line == 'S1:':
            group_flag = 0
            continue
        if line == 'S2:':
            group_flag = 1
            continue
        if line == 'G1:':
            group_flag = 2
            continue
        if line == 'G2:':
            group_flag = 3
            continue

        if group_flag == 0:
            space_group1.append(line)
        if group_flag == 1:
            space_group2.append(line)
        if group_flag == 2:
            ground_group1.append(line)
        if group_flag == 3:
            ground_group2.append(line)

        sgs = [space_group1, space_group2]
        ggs = [ground_group1, ground_group2]

    return sgs, ggs

#Truncate (rnd down) all values to int
def scale_data(data):
    for i in data:
        if i == 'Geneid':
            continue
        data[i] = pd.to_numeric(data[i]*100, errors='raise').astype(int)

    return data

def remove_spikeins(data):
    truth = ~data['Geneid'].str.contains('ERCC')
    filtered_data = data[truth].reset_index(drop=True)

    return filtered_data

#If there are any that are all 0s in one group
#and non 0s in the other extract them as significant
def extract_signif_gene(data, sg, gg):
    sig_rows = []
    for i in range(len(data)):
        check1 = (sum(data.loc[i][sg]) == 0) and (not 0 in data.loc[i][gg].values)
        check2 = (sum(data.loc[i][gg]) == 0) and (not 0 in data.loc[i][sg].values)
        if check1 or check2:
            sig_rows.append(i)

    sig_genes = data.loc[sig_rows].reset_index(drop=True)
    filtered_data = data.drop(sig_rows, axis=0).reset_index(drop=True)

    return filtered_data, sig_genes

#Go row by row, delete rows where 90% of counts are in one sample of one group
def remove_over_representation(data, sg, gg):
    rows_to_drop = []
    for i in range(len(data)):

        #for each line check to see if either the ground group or space group
        #is overrepresented by one sample, but ignores overrepresented samples
        #where the sum is 10 counts or fewer (e.g [0,0,0,0,10])
        gg_rep = (max(data.loc[i][gg]) >= sum(data.loc[i][gg])*0.9)
        sg_rep = (max(data.loc[i][sg]) >= sum(data.loc[i][sg])*0.9)

        #if either group is overrepresented, remove the row from the data
        if gg_rep or sg_rep:
            rows_to_drop.append(i)

    filtered_data = data.drop(rows_to_drop, axis=0).reset_index(drop=True)
    over_rep_data = data.loc[rows_to_drop].reset_index(drop=True)

    return filtered_data, over_rep_data




#make new filtered count file(s)
def main(dfile, config, output):
    sgs, ggs = read_config(config)
    for i in range(2):
        cols = ['Geneid'] + sgs[i] + ggs[i]
        data = pd.read_csv(dfile, usecols=cols)
        data = scale_data(data)
        data = remove_spikeins(data)
        data, sig_genes = extract_signif_gene(data, sgs[i], ggs[i])
        data, over_rep = remove_over_representation(data, sgs[i], ggs[i])
        data.to_csv(f'{output}/group{i+1}_filtered_gene_counts_scalex100.csv', index=False)
        sig_genes.to_csv(f'{output}/group{i+1}_automatic_significant_genes.csv', index=False)
        over_rep.to_csv(f'{output}/group{i+1}_over_represented_data.csv', index=False)


if __name__ == '__main__':
    args = arg_parser()
    main(args.fname, args.config, args.output)