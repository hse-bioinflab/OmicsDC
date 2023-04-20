import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import multiprocessing as mp
from multiprocessing import Queue, Process, current_process,active_children, Manager, shared_memory
from tqdm import tqdm
import subprocess
import time
import wget

cmd_line = argparse.ArgumentParser(description='Script for preparing assembly file')

cmd_line.add_argument(
    '--assembly',
    '-i',
    type=str,
    default='hg19',
    help ="assembly name via chip atlas"
)
cmd_line.add_argument(
    '--check',
    '-c',
    type = int,
    default= 0,
    help = "checking files"
)
cmd_line.add_argument(
    '--nworkers',
    '-n',
    default=6
)
cmd_line.add_argument(
    "--chunk_size",
    '-c',
    type = int,
    default=100_000
)



#TODO Подумать про словарь
# exp = table.loc[table[1] == "hg38"]
# exp = np.array(exp.values.tolist())
# shm = shared_memory.SharedMemory(create=True, size=exp.nbytes)
# exp_copy = np.ndarray(exp.shape, dtype=exp.dtype, buffer=shm.buf)
# exp_copy[:] = exp[:]


# existing_shm = shared_memory.SharedMemory(name=shm.name)
# c = np.ndarray(exp.shape, dtype=exp.dtype, buffer=existing_shm.buf)
    
def SendPackage(que: list, parser_bufer: dict, pause: dict):
    for buf in parser_bufer:
        if len(parser_bufer[buf]) > 100:
            while pause[buf] > 10:
                if (que[buf].qsize() > 1000):
                    pause[buf] += 1
                else:
                    pause[buf] = 0
                    que[buf].put(parser_bufer[buf])
                    parser_bufer[buf] = []


def parser(chunk_que, que: dict, num_writers: int):
    """Process take chunk from chunk_que and line by line putting it into ParseLine funvtion"""
    parser_bufer = {key: [] for key in range(num_writers)}
    buf_trigger = 0
    pause = {key: 0 for key in range(num_writers)}
    while not chunk_que.empty():
        data = chunk_que.get()
        for i in range(data.shape[0]):
            line = data.iloc[i,:]
            q_num = int(line[3][3:])%num_writers
            parser_bufer[q_num].append(line)
            buf_trigger += 1
            if buf_trigger > 1000:
                buf_trigger = 0
                SendPackage(que, parser_bufer, pause)

def writer(que: Queue, file_names: dict):
    pause = 0
    verbose = 0
    is_working = True
    while is_working:
        if not que.empty():
            pause = 0
            lines = que.get()
            for line in lines:
                with open(f"./{current_process().name}.txt", "+a") as f:
                    f.write('\t'.join(map(str, line)) + '\n')
                    verbose += 1
            
            if verbose > 50000:
                print(f"{current_process().name} working on {line}")
                verbose = 0
        else:
            if pause > 20:
                print(f"{current_process().name} paused to {pause}")
            if pause > 100:
                #через секунду офается нафиг
                is_working = False
            time.sleep(0.01)
            pause += 1

    print(f"{current_process().name} stopped")
    

#TODO genome folder
def worker_file_creator(df,tasks, check:bool, file_dict: dict):
    while not tasks.empty():
        start, end = tasks.get()
        for f in range(start, end):
            file = list(df.iloc[f])
            path = str(Path("./resources"))+f'/{file[0]}_{file[1]}_{file[2]}_{file[3]}_{file[4]}_{file[5]}.bed'
            if check:
                if not Path.is_file(path):
                    file_dict["No_file"] += 1
            else:
                with open(path, 'w+') as f:
                    file_dict[file[0]] = path


def ExpListProcessing(data: pd.DataFrame, n_workers: int, check: bool, file_dict: dict):
    tasks = Queue()
    # if check ?
    for batch_start in range(0, len(data) - 1000, 1000):
        tasks.put((batch_start, batch_start + 1000))
    tasks.put((len(data) - 1000, len(data)))
    
    procs = [Process(target=worker_file_creator, args=(data, tasks, check, file_dict, )) for _ in range(n_workers) ]

    [p.start()  for p in procs]
    [p.join()   for p in procs]

    return 0
    


#TODO аргумент лист
if __name__ == '__main__':
#Как будто бы можно вьбать загрузку в этом модуле
    args = cmd_line.parse_args()
    NWORKERS = args.nworkers
    NWRITERS = NWORKERS * 9
    NPARSERS = NWORKERS * 3
    file_dick = Manager().dict()

    if not Path.is_file('./resources/experimentList.tab'):
        subprocess.run(f"gunzip {Path('./resources/experimentList.tab')}")
    df = pd.read_csv(Path("./resources/experimentList.tab"),
                        sep = '\t', 
                        usecols=range(6),
                        header = None
                        )
    df = df[df[1] == args.assembly]
    df .replace(
        [' ','/'], '_',
        inplace=True, regex=True
        )
    
    
    if not Path.is_dir(f"./resources/{args.assembly}"):
        print("new wget")
        wget.download(f"https://chip-atlas.dbcls.jp/data/hg38/allPeaks_light/allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed.gz",
                      Path(f"./resources")
                    )
        subprocess.run(f"gunzip ./resources/allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed.gz")
        subprocess.run(f"mkdir ./resources/{args.assembly}")
        
        ExpListProcessing(df, n_workers=8, check = False, file_dict = file_dick)
    
        parser_list = []
        writer_list = []
        workers_ques = [ Queue() for i in range(NWRITERS) ]
        chunk_que = Queue()

        iter_df = pd.read_csv(
                Path(f"./resources/{args.assembly}"),
                chunksize=1_000,
                sep = '\t',
                header= None
            )
        
        for chunk in iter_df:
            for line in chunk:
                    with open(file_dick[line[0]], "+a") as f:
                        f.write('\t'.join(map(str, line)) + '\n')
                        
        subprocess.run(f"rm ./resources/allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed.gz")
    
    elif args.check:
        file_dick["No_file"] = 0 
        result = ExpListProcessing(df, n_workers=8, check = False, file_dict = file_dick)
        print(file_dick["No_file"])
    
    else:
        print("Already satisfaied")    
