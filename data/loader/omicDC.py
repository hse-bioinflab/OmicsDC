from pathlib import Path
from typing import Union
from time import gmtime, strftime

import subprocess
import pandas as pd
import numpy as np

import multiprocessing
import os
import subprocess

def create_matching_expirement_df(
            filepath,
            options
        ) -> pd.DataFrame:
    """ Function to return expirement names df"""

    # match_exp_df - df for matching experiments
    match_exp_df = pd.read_csv(
                    filepath,
                    sep = ',', 
                    names = ['id', 'Genome assembly', 'Antigen class', 'Antigen', 'Cell type class', 'Cell type']
                    #usecols=range(6)
                )

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



def omics(expid: list = None, assembly: list = ['hg38'], assembly_threshold: list = '05' , antigen_class: list = None,
          antigen: list = None, cell_type: list = None, cell: list = None, resources: Path = Path("./loader"), output_path: Path = Path("./storage/")):
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
    options = {
        "id"                :   expid,
        "Genome assembly"   :   assembly,
        "Antigen class"     :   antigen_class,
        "Antigen"           :   antigen,
        "Cell type class"   :   cell_type,
        "Cell type"         :   cell
    }

    #check files
    for a in assembly:
            SubporocessHub = subprocess.Popen(["python","-W","ignore","./loader/file_creator.py","-a", a, "-t", assembly_threshold, "-c", "1"], stdout= subprocess.PIPE)
            SubporocessHub.wait()
            Result = str(SubporocessHub.communicate())
            
            if "reload" in Result:
                print(f"Start reload {a}")
                #SubporocessHub = subprocess.Popen(["python","-W","ignore","./loader/file_creator.py","-a", a, "-t", assembly_threshold, "-r", "1", "-d", "1"],stdout=subprocess.PIPE)
                #SubporocessHub.wait()
            else:
                print(f"File check of {a} complete")




    match_exp_df = create_matching_expirement_df(Path("./loader/resources/experimentList.tab"), options)
    files_string = ''
    for a in assembly:
        files_dir = str(resources / f"resources/{a}/")
        match_exp_l = f'_*,'.join(match_exp_df["id"].loc[match_exp_df["Genome assembly"] == a].to_numpy()) + "_*"
        files_string+=match_exp_l
    print(f"Copying {files_string.count('*')}  files")
        #SubporocessHub = subprocess.Popen("cp ./"+files_dir+'/' + "{" +match_exp_l + "} ./storage/",shell = True)
    filename_d = strftime("%Y-%m-%d_%H:%M:%S", gmtime())
    if files_string.count('*') == 1:
        SubporocessHub = subprocess.Popen(f"tar -czf {output_path}/{filename_d}.tar.gz ./"+files_dir+'/' + match_exp_l,shell = True)
    else:
        SubporocessHub = subprocess.Popen(f"tar -czf {output_path}/{filename_d}.tar.gz ./"+files_dir+'/' + "{" +match_exp_l + "}",shell = True)
    
    SubporocessHub.wait()
    print("Process done")
