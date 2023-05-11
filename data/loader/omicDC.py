from pathlib import Path
from typing import Union
from time import gmtime, strftime
from loader import FILE_CREATOR,RESOURCES

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

def omics(expid: list = None, assembly: list = ['hg38'], assembly_threshold: str = '05' , antigen_class: list = None,
          antigen: list = None, cell_type: list = None, cell: list = None, output_path: Path = Path("./storage/")) -> Path:
    """
    Function for processing omics data.

    Args:
        expid (list): List of experiment IDs.
        assembly (list): List of genome assemblies.
        assembly_threshold (str): Assembly threshold.
        antigen_class (list): List of antigen classes.
        antigen (list): List of antigens.
        cell_type (list): List of cell type classes.
        cell (list): List of cell types.
        output_path (Path): Output path for the processed files.

    Returns:
        Path: Path to the processed files.

    """

    options = {
        "id"                :   expid,
        "Genome assembly"   :   assembly,
        "Antigen class"     :   antigen_class,
        "Antigen"           :   antigen,
        "Cell type class"   :   cell_type,
        "Cell type"         :   cell
    }

    # Check files for each assembly
    for a in assembly:
            Result = ""
            while not "done" in Result:
                # Execute the file check subprocess
                SubporocessHub = subprocess.Popen(["python",FILE_CREATOR,"-a", a, "-t", assembly_threshold, "-c", "1"], stdout= subprocess.PIPE)
                SubporocessHub.wait()
                Result = str(SubporocessHub.communicate())
                
                if "dir" in Result:
                    # Create working directory with all files download
                    print(f"Start creating working dir {a} with all files download? y/[N]")
                    user_answer = input()
                    if user_answer == 'y':
                        SubporocessHub = subprocess.Popen(["python","-W","ignore",FILE_CREATOR,"-a", a, "-t", assembly_threshold, "-r", "1", "-d", "1"],stdout=subprocess.PIPE)
                        SubporocessHub.wait()
                    else:
                        print("Process canceled")
                        return(0)

                elif "reload" in Result:
                    # Create local files
                    print(f"Start creating local files of {a}")
                    SubporocessHub = subprocess.Popen(["python","-W","ignore",FILE_CREATOR,"-a", a, "-t", assembly_threshold, "-r", "1"],stdout=subprocess.PIPE)
                    SubporocessHub.wait()
                    Result = str(SubporocessHub.communicate())
                    print(Result)

                elif "download" in Result:
                    # Download files
                    print(f"Need to download {a}. Start process? y/[N]")
                    user_answer = input()
                    if user_answer == 'y':
                        SubporocessHub = subprocess.Popen(["python","-W","ignore",FILE_CREATOR,"-a", a, "-t", assembly_threshold, "-r", "1", "-d", "1"],stdout=subprocess.PIPE)
                        SubporocessHub.wait()
                    else:
                        print("Download cancel")
                        return(0)
                else:
                    print(f"File check of {a} complete")



    # Create a dataframe with matching experiment information
    match_exp_df = create_matching_expirement_df(RESOURCES / "experimentList.tab", options)
    files_string = ''
    for a in assembly:
        files_dir = str(RESOURCES / f"resources/{a}/")
        match_exp_l = f'_*,'.join(match_exp_df["id"].loc[match_exp_df["Genome assembly"] == a].to_numpy()) + "_*"
        files_string+=match_exp_l
    print(f"Copying {files_string.count('*')}  files")

    filename_d = strftime("%Y-%m-%d_%H:%M:%S", gmtime())
    if files_string.count('*') == 1:
        # Compress files into a tar.gz archive
        SubporocessHub = subprocess.Popen(f"tar -czf {output_path}/{filename_d}.tar.gz ./"+files_dir+'/' + match_exp_l,shell = True)
    else:
        SubporocessHub = subprocess.Popen(f"tar -czf {output_path}/{filename_d}.tar.gz ./"+files_dir+'/' + "{" +match_exp_l + "}",shell = True)
    
    SubporocessHub.wait()
    print(f"Done. Created file ./{output_path}/{filename_d}.tar.gz ")

    return(output_path / f"{filename_d}.tar.gz")
