'''
    This is script reformat result files from MSFragger for peptide elution analysis (McWhite et al., 2022).
'''
import os
import pandas as pd
import numpy as np
import glob

#need to pull out peptide.tsv file to make a wide format input table
#format of the peptide table is as follows:
#4 columns ExperimentID FractionID  Peptide PeptideCount
def format_files(root_folder, fractionation_name):
    
    #Each fraction has its own folder from Fragger. Compile all the folders.
    
    folder_list=glob.glob(f'{root_folder}/*/')
    folder_list.sort()
    peptide_filelist = []
    for folder in folder_list:
        os.chdir(folder)
        
        #format file 
        ExperimentID = fractionation_name
        FractionID = folder.split('/')[-2]
        peptide = pd.read_csv(os.path.join(folder, 'peptide.tsv'), 
                              sep ='\t', usecols=['Peptide', 'Spectral Count'])
        peptide = peptide.rename(columns={'Spectral Count':'PeptideCount'})
        peptide['ExperimentID'] = ExperimentID
        peptide['FractionID'] = FractionID
        cols = peptide.columns.tolist()
        cols = cols[2:4] + cols[0:2]
        peptide = peptide[cols]
        peptide_filelist.append(peptide)
        #concat all peptide files into one single dataframe
        combined_df = pd.concat(peptide_filelist)
    
    return combined_df

def pivot_file (input_file, fractionation_name):
    
    #ProtID, TotalCount, fraction1...n
    pivot = input_file.pivot(index= 'Peptide', columns='FractionID', values="PeptideCount")
    pivot = pivot.reset_index(level=0)
    pivot = pivot.fillna(0)
    pivot['ExperimentID'] = fractionation_name
    
    #pivot = pivot.drop(columns=['FractionID'])
    pivot['Peptide'] = pivot['Peptide'].str.replace('I', 'J')
    pivot['Peptide'] = pivot['Peptide'].str.replace('L', 'J')
    
    #create fraction order file
    fraction_order = pd.DataFrame({'FractionOrder':list(range(1,len(list(pivot.columns)[1:-1])+1)),
                               'FractionID':list(pivot.columns)[1:-1]})
    
    
    
    return pivot, fraction_order


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Reformat result files from MSFragger for peptide elution analysis')

    
    parser.add_argument('-r', '--root_folder', dest='root_folder', action="store", required=True, type=str, help= 'folder where MSFragger write result files and folders')

    parser.add_argument('-n', '--fractionation_name', dest='fractionation_name', action="store", required=True, type=str, help= 'name of fractionation experiment')

    parser.add_argument('-o', '--output_file', dest= 'output_file', action="store", type=str, help= 'name of output file in the format ready for peptide elution analysis')

    parser.add_argument('-f', '--fraction_order', dest= 'fraction_order', action="store", type=str, help= 'name of output fraction order file')
    
    args = parser.parse_args()
    combined_df = format_files(args.root_folder, args.fractionation_name)

    #pivot file to wide format and create fraction order file
    pivot, fraction_order = pivot_file (combined_df, args.fractionation_name)


    pivot.to_csv(args.output_file, index=False)
    fraction_order.to_csv(args.fraction_order, index=False)
    


