from pathlib import Path
from loader.resources import cashed_ExperimentList
from typing import Union

import subprocess
import pandas as pd
import numpy as np
import dask.dataframe as dd
from dask.distributed import Client
import wget

import multiprocessing
import os
import warnings

def create_matching_expirement_df(
            filepath,
            options,
            verbose
        ):
    """ Function to return expirement names df"""

    # match_exp_df - df for matching experiments
    match_exp_df = pd.read_csv(
                    filepath,
                    sep = '\t', 
                    names = ['id', 'Genome assembly', 'Antigen class', 'Antigen', 'Cell type class', 'Cell type'],
                    usecols=range(6)
                )
    if verbose:
        print("Find file " +  filepath)

    for key in options.keys():
        if options[key]:
            match_exp_df = match_exp_df.loc[match_exp_df[key].isin(options[key])]

    return match_exp_df


def add_sorted_bed_2_file(
            df,
            matching_experiments
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


def omics(expid: Union[str,list] = None, assembly: Union[str,list] = 'hg38', assembly_threshold: Union[str,list] = '05' , antigen_class: Union[str,list] = None,
          antigen: Union[str,list] = None, cell_type: Union[str,list] = None, cell: Union[str,list] = None, storage: Path = 'storage', output_path: Path = './'):
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


    Outputs:
        In Path functions creates .csv file with omic data

    Returns:
        No return
    '''
    
    storage = del_slash(storage)
    output_path = del_slash(output_path)

    # move exp file
    if not os.path.isfile(storage + "/experimentList.tab"):
        subprocess.run(f"gunzip {storage}/experimentList.tab.gz")

    if not os.path.isfile(storage + f"/allPeaks_light.{assembly}.{assembly_threshold}.bed"):
        print("new wget")
        wget.download(f"https://chip-atlas.dbcls.jp/data/hg38/allPeaks_light/allPeaks_light.{assembly}.{assembly_threshold}.bed.gz",
                      storage
                    )
        subprocess.run(f"gunzip {storage}/allPeaks_light.{assembly}.{assembly_threshold}.bed.gz")
    
    options = {
        "id"                :   expid,
        "Genome assembly"   :   assembly,
        "Antigen class"     :   antigen_class,
        "Antigen"           :   antigen,
        "Cell type class"   :   cell_type,
        "Cell type"         :   cell
    }


# Dask multiuser checking.    
    # with warnings.catch_warnings(record=True) as caught_warnings:
    #     warnings.simplefilter("always")
    #     for warn in caught_warnings:
    #         if str(warn.message).find('Port 8787 is already in use') != -1:
    #             print(f"{bcolors.OKCYAN}U r not alone. Sorry but u have to w8.\nChill a bit!{bcolors.ENDC}") 
    #             exit()

    match_exp_df = create_matching_expirement_df(storage + "/experimentList.tab", options)

    create_sorted_bed_file(f"allPeaks_light.{assembly}.{assembly_threshold}.bed", match_exp_df)
    
    os.replace(storage + f"/filtred_allPeaks_light.{assembly}.{assembly_threshold}.bed", output_path)
    subprocess.run(f"gzip {output_path}/filtred_allPeaks_light.{assembly}.{assembly_threshold}.bed.gz")
