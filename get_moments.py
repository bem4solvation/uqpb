import os
import numpy 
import pandas
import sys
import glob
import argparse
from time import time
import gc

def check_parser(argv):
    """ 
    Reads the inputs from the command line
    """

    parser = argparse.ArgumentParser(description="Calculate moments for results in a csv file of a folder")
    parser.add_argument('-f','--folder', dest='folder', type=str, default=None, help='Folder with csv files')

    args = parser.parse_args(argv)

    if args.folder[-1] != "/":
        args.folder += "/"

    return args.folder


def get_moments(folder, columns_to_process =["solv_energy","elec_energy","cav_energy","disp_energy","nonpolar_energy"]):

    output_file = folder + "output_summary.csv"
    output_file_consolidated = folder + "output_consolidated.csv"

    #results_df = pandas.DataFrame(columns=['File','Mean','Standard deviation'])

    csv_files = glob.glob(folder + "**/*.csv")

    if output_file in csv_files:
        csv_files.remove(output_file)

    dfs = []
    dfs_cons = []
    for input_csv in csv_files:

        df = pandas.read_csv(input_csv)

        mean_values = df[columns_to_process].mean()
        std_dev_values = df[columns_to_process].std()

        local_result_df = pandas.DataFrame({
            'File': [input_csv],
            'Mean': [mean_values.to_dict()],
            'Standard Deviation': [std_dev_values.to_dict()]
        })

        dfs.append(local_result_df)
        dfs_cons.append(df)

    results_df = pandas.concat(dfs, ignore_index=True)
    results_df.to_csv(output_file, index=False)

    results_cons = pandas.concat(dfs_cons, ignore_index=True)
    
    # For avoiding warning on transform str to float, drop unnecesary columns
    results_cons = results_cons.drop(columns=['iter', 'pqr_file']) 
    
    mean_row = results_cons.mean()
    std_row = results_cons.std()
    
    results_cons.loc['Mean'] = mean_row
    results_cons.loc['Std_dev'] = std_row

    results_cons.to_csv(output_file_consolidated, index=False)
if __name__ == "__main__":
    
    folder  = check_parser(sys.argv[1:])
    get_moments(folder)
