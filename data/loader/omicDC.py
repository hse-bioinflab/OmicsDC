from pathlib import Path

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
            filepath, 
            options,
            verbose
        ):   
    """ Function to return expirement names df"""

    # match_exp_df - df for matching experiments
    match_exp_df = pd.read_csv(
                    filepath,
                    filepath,
                    sep = '\t', 
                    names = ['id', 'Genome assembly', 'Antigen class', 'Antigen', 'Cell type class', 'Cell type'],
                    usecols=range(6)
                )
    if verbose:
        print("Find file " +  filepath)

    for key in options.keys():
        if options[key]:
            tmp = options[key].split(',')
            match_exp_df = match_exp_df.loc[match_exp_df[key].isin(tmp)]

    return match_exp_df


def add_sorted_bed_2_file( 
            filename,
            df,
            num,
            matching_experiments,
        ):
    """ Function to add lines to .csv file from part of sorted .bed files"""
    return df.loc[df['id'].isin(matching_experiments)]


def create_sorted_bed_file(
        storage,
        filename,
        match_exp_df
    ):
    """Create big .csv table with every finded match"""

    path_2_sorted_file = storage + "filtred_" + filename + ".bed"

    process_list = []

    matching_experiments = list(match_exp_df.loc[:,'id'])

    df = dd.read_csv(   
                storage + filename,
                sep = "\t", 
                names = ['chr', 'begin', 'end', 'id', 'score'],
                blocksize = '100mb'
                )

    df = df.map_partitions(add_sorted_bed_2_file, matching_experiments)
    df.to_csv(path_2_sorted_file, index=False, header=False, single_file=True)


def del_slash(path):
    return path[:-1] if path[-1] == '/' else path


def omics(expid: str = None, assembly: str = 'hg38', assembly_threshold: str = '05' , antigen_class: str = None, antigen: str = None, cell_type: str = None, cell: str = None, storage: Path = 'storage', output_path: Path = './' , ncores: int = 2, nworkers: int = 4, verbose: bool = True):
    '''
    Function to create omics data from chip-atlas database.
    Arguments:
        expid: str of experiments Id
        assembly_threshold: str 
        antigen_class: str
        antigen: str
        cell_type: str
        cell: str
        storage: Path, default is './data/storage'
        output_path: Path, default is './'
        assembly: str, default is 'hg38'
        ncores and nworkers: int, multithread parametrs
        verbose: bool

    Outputs:
        In Path functions creates .csv file with omic data

    Returns:
        No return
    '''
    
    storage = del_slash(storage)
    output_path = del_slash(output_path)

    # move exp file
    if not os.path.isfile(storage + "/experimentList.tab"):
        os.system(f"gunzip {storage}/experimentList.tab.gz")

    if not os.path.isfile(storage + f"/allPeaks_light.{assembly}.{assembly_threshold}.bed"):
        print("new wget")
        wget.download(f"https://chip-atlas.dbcls.jp/data/hg38/allPeaks_light/allPeaks_light.{assembly}.{assembly_threshold}.bed.gz",
                      storage
                    )
        os.system(f"gunzip {storage}/allPeaks_light.{assembly}.{assembly_threshold}.bed.gz")
    
    options = {
        "id"                :   expid,
        "Genome assembly"   :   assembly,
        "Antigen class"     :   antigen_class,
        "Antigen"           :   antigen,
        "Cell type class"   :   cell_type,
        "Cell type"         :   cell
    }

    for key in options.keys():
        if options[key]:
            options[key] = options[key].replace('_', ' ')
    
    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")
        for warn in caught_warnings:
            if str(warn.message).find('Port 8787 is already in use') != -1:
                print(f"{bcolors.OKCYAN}U r not alone. Sorry but u have to w8.\nChill a bit!{bcolors.ENDC}") 
                exit()

    match_exp_df = create_matching_expirement_df(storage + "/experimentList.tab", options)

    create_sorted_bed_file(f"allPeaks_light.{assembly}.{assembly_threshold}.bed", match_exp_df)
    
    os.replace(storage + f"/filtred_allPeaks_light.{assembly}.{assembly_threshold}.bed", output_path)
    os.system(f"gzip {output_path}/filtred_allPeaks_light.{assembly}.{assembly_threshold}.bed.gz")
