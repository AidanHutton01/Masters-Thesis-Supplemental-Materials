import argparse
import pandas as pd
import glob


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="The path of the folder containing the files of differential expression data you want to filter")
    return parser.parse_args()


def annotate(path):
    dirs = glob.glob(f'{path}/*_filtered.csv')

    annotation = pd.read_csv(f'{path}/Gene_Annotation.txt')
    for f in dirs:
        DE_data = pd.read_csv(f"{f}")
        DE_data['Chromosome/scaffold name'] = 'N/A'
        DE_data['Gene start (bp)'] = 'N/A'
        DE_data['Gene end (bp)'] = 'N/A'
        DE_data['Strand'] = 'N/A'
        DE_data['Gene name'] = 'N/A'
        DE_data['Gene description'] = 'N/A'
        
        for i in range(len(DE_data)):
            try:
                index = annotation['Gene stable ID'][annotation['Gene stable ID'] == DE_data['GeneID'][i]].index[0]
                DE_data.loc[i, 'Gene start (bp)'] = annotation['Gene start (bp)'][index]
                DE_data.loc[i, 'Gene end (bp)'] = annotation['Gene end (bp)'][index]
                DE_data.loc[i, 'Strand'] = annotation['Strand'][index]
                DE_data.loc[i, 'Gene name'] = annotation['Gene name'][index]
                DE_data.loc[i, 'Chromosome/scaffold name'] = annotation['Chromosome/scaffold name'][index]
                DE_data.loc[i, 'Gene description'] = annotation['Gene description'][index]
            except:
                print('Error')
        DE_data.to_csv(f'{f}_annotated.csv', index=False)
        
if __name__ == '__main__':
    args = arg_parser()
    annotate(args.path)