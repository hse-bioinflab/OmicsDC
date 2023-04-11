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


#defines
#TODO Worker num and cores num
PRIVATE_PATH = "/home/avoitetskii/private_omicON.txt"

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
        storage,
        filename,
        match_exp_df
    ):
    """Create big .csv table with every finded match"""

    path_2_sorted_file = FILE_PATH + "filtred_" + filename + ".csv"

    process_list = []

    matching_experiments = list(match_exp_df.loc[:,'id'])

    df = dd.read_csv(   
                storage + filename,
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


def omics(expid: str = None, assembly: str = 'hg38', assembly_threshold: str = None , antigen_class: str = None, antigen: str = None, cell_type: str = None, cell: str = None, storage: Path = './', output_path: Path = './' , ncores: int = 2, nworkers: int = 4):
    '''
    Function to create omics data from chip-atlas database.
    Arguments:
        expid: str of experiments Id
        assembly_threshold: str 
        antigen_class: str
        antigen: str
        cell_type: str
        cell: str
        storage: Path, default is './'
        assembly: str, default is 'hg38'
        ncores and nworkers: int, multithread parametrs

    Outputs:
        In Path functions creates .csv file with omic data

    Returns:
        No return
    '''
     
    NCORES    = ncores
    NWORKERS  = nworkers
    FILE_PATH = storage

    # move exp file
    if not os.path.isfile(FILE_PATH + f"/experimentList.tab"):
        os.replace("./data/storage/experimentList.tab.gz", FILE_PATH)
        os.system(f"gunzip {FILE_PATH}/experimentList.tab.gz")

    if not os.path.isfile(FILE_PATH + f"/allPeaks_light.{assembly}.{assembly_threshold}.bed"):
        wget.download(
                        f"https://chip-atlas.dbcls.jp/data/hg38/allPeaks_light/\
                        allPeaks_light.{assembly}.{assembly_threshold}.bed.gz",
                        FILE_PATH
                      )
        os.system(f"gunzip {FILE_PATH}/allPeaks_light.{assembly}.{assembly_threshold}.bed.gz")
        
    
    options = {
        #Parse arguments from cmd line to special dict
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
        que = Client(n_workers=NCORES, threads_per_worker=NWORKERS)
        for warn in caught_warnings:
            if str(warn.message).find('Port 8787 is already in use') != -1:
                print(f"{bcolors.OKCYAN}U r not alone. Sorry but u have to w8.\nChill a bit!{bcolors.ENDC}") 
                exit()

    match_exp_df = create_matching_expirement_df("experimentList.tab", options)

    create_sorted_bed_file(que, f"allPeaks_light.{assembly}.{assembly_threshold}.bed", match_exp_df)
    que.shutdown()
    
    match_exp_df.to_csv(
                        "result.csv.gz",
                        index=False, 
                        compression="gzip"
                        )
    
    os.remove(FILE_PATH + "filtred_" + f"allPeaks_light.{assembly}.{assembly_threshold}.bed" + ".csv")


def assembly(tag: str, saveto: Path, *_, force: bool = False):
    supported = {"GRCh38", "GRCm39"}
    assert tag in supported, f"Requested assembly ({tag}) is not among supported: {','.join(supported)}"

    # TODO: clearer error message
    assert saveto.name.endswith(".gz"), "Loaded assemblies must be saved as indexed & bgzip-ed files"

    # Use GENCODE:
    # After downloading, unzip -> bgzip -> samtools faidx files
    # GRCh38 - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
    # GRCm39 - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz
    raise NotImplementedError()


__all__ = ["omics", "assembly"]