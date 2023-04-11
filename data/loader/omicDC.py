import pandas as pd
import numpy as np
import dask.dataframe as dd
import stat
import time
import wget
import dask
from dask.distributed import Client
from dask.diagnostics import ProgressBar
from dask.distributed import progress
import gc

import multiprocessing
import os
import warnings


#defines
#TODO Worker num and cores num
PRIVATE_PATH = "private_omicON.txt"

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def create_matching_expirement_df(
            que, 
            filename, 
            options
        ):   
    """ Function to return expirement names list"""

    # match_exp_df - df for matching experiments
    match_exp_df = pd.read_csv(
                                FILE_PATH + filename,
                                sep = '\t', 
                                names = ['id', 'Genome assembly', 'Antigen class', 'Antigen', 'Cell type class', 'Cell type'],
                                usecols=range(6)
                            )
    
    if args.verbose:
        print("Find file " +  FILE_PATH + filename)

    for key in options.keys():
        if options[key]:
            tmp = options[key].split(',')
            match_exp_df = match_exp_df.loc[match_exp_df[key].isin(tmp)]

    return match_exp_df


def add_user_bed_markers(
        que,
        filename
    ):
    """ Merging sorted df with user`s .bed file as two additional cols
        Saving onli rows with several chr.
        Creating new 'intersect' column with booleans.
        Returning df with only intersected"""
    path_2_sorted_file_with_user_bed = FILE_PATH + "filtred_and_bed" + filename + ".csv"

    df = dd.read_csv(
        FILE_PATH + "filtred_" + filename + ".csv", 
        header=None, 
        sep=',',
        names = ['chr', 'begin', 'end', 'id', 'score'],
        blocksize = '10mb'
        )
    
    process_list = []
    
    df.set_index('chr')

    for part in range(df.npartitions):
        process_list.append(que.submit(
                                make_intersect,
                                df,
                                part,
                                path_2_sorted_file_with_user_bed
                                )) 
    
    a = [process.result() for process in process_list]

    return df.compute()


def add_sorted_bed_2_file( 
            filename,
            df,
            num,
            matching_experiments,
        ):
    """ Function to add lines to .csv file from part of sorted .bed files"""
    part = df.partitions[num]
    part = part.loc[part['id'].isin(matching_experiments)]
    part = part.compute()
    part.to_csv(filename, index=False, header=False, mode='a')
    return num


def create_sorted_bed_file(
        que,
        filename,
        match_exp_df
    ):
    """Create big .csv table with every finded match"""

    path_2_sorted_file = FILE_PATH + "filtred_" + filename + ".csv"

    process_list = []

    matching_experiments = list(match_exp_df.loc[:,'id'])

    df = dd.read_csv(   
                        FILE_PATH + filename,
                        sep = "\t", 
                        names = ['chr', 'begin', 'end', 'id', 'score'],
                        blocksize = '50mb'
                    )

    open(path_2_sorted_file, mode = 'w').close()  # Creating empty .csv for editing
    os.chmod(path_2_sorted_file, 33279)

    for part in range(df.npartitions):
        process_list.append(que.submit(
                                add_sorted_bed_2_file,
                                path_2_sorted_file,
                                df,
                                part,
                                matching_experiments
                                ))
    
    a = [process.result() for process in process_list]
    progress(a, notebook = False)


def parse_private():
    """Function to take data from private .txt file"""
    d = {}
    with open(PRIVATE_PATH) as f:
        for line in f:
            if str(line) == '___Doc_list___\n':
                if args.verbose:
                    print('___Doc_list___')
                continue
            (key, val) = line.split()
            d[str(key)] = val

            if args.verbose:
                print(key, ':', val)
    return(d)