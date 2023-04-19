import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import multiprocessing as mp
from multiprocessing import Queue, Process, current_process,active_children, Manager, shared_memory
from tqdm import tqdm
import os
import queue
import time

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
                with open(file_names[lines[0]], "+a") as f:
                    f.write(str(line))
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
            if check:
                # проверка наличия файла с таким именем
                pass
            else:
                with open(str(Path("./resources"))+f'/{file[0]}_{file[1]}_{file[2]}_{file[3]}_{file[4]}_{file[5]}.bed', 'w+') as f:
                    file_dict[file[0]] = f'{file[0]}_{file[1]}_{file[2]}_{file[3]}_{file[4]}_{file[5]}.bed'


def ExpListProcessing(data: pd.DataFrame, n_workers: int, check: bool, file_dict: dict):
    tasks = Queue()
    # if check ?
    for batch_start in range(0,1000,1000):
    #range(0, len(data) - 1000, 1000):
        tasks.put((batch_start, batch_start + 1000))
    tasks.put((len(data) - 1000, len(data)))
    
    procs = [Process(target=worker_file_creator, args=(data, tasks,check, file_dict )) for _ in range(n_workers) ]

    [p.start()  for p in procs]
    [p.join()   for p in procs]
    


#TODO аргумент лист
if __name__ == '__main__':
#Как будто бы можно вьбать загрузку в этом модуле
    args = cmd_line.parse_args()
    NWORKERS = args.nworkers
    NWRITERS = NWORKERS * 9
    NPARSERS = NWORKERS * 3
    file_dick = Manager().dict()
    print(type(file_dick))

    # iterator = pd.read_csv(
    #         f'M:/Fast_Work/allPeaks_light.mm9.50.bed',
    #         chunksize=1_000,
    #         sep = '\t',
    #         header= None
    #     )

    df = pd.read_csv(Path("./resources/exp_edited.tab"),
                    sep = ',', 
                    #usecols=range(5),
                    header = None
                    )
    print(args.assembly)
    df = df[df[1] == args.assembly]
    df .replace(
        [' ','/'], '_',
        inplace=True, regex=True
    )
    print(df.head())

    ExpListProcessing(df, n_workers=8, check = args.check, file_dict = file_dick)
 
    parser_list = []
    writer_list = []
    workers_ques = [ Queue() for i in range(NWRITERS) ]
    chunk_que = Queue()


    for chunk in tqdm(iterator, total = 10):
        chunk_que.put(chunk)
        while chunk_que.qsize() > NPARSERS:
            if len(parser_list) == 0:
                for i in range(NWRITERS):
                    p = Process(target=writer, args=(workers_ques[i],))
                    writer_list.append(p)
                    p.start()
                for chunk,_ in zip(iterator,range(NPARSERS)):
                    p = Process(target=parser, args=(chunk_que, workers_ques,NWRITERS))
                    parser_list.append(p)
                    p.start()
                
            pass

    [par.join() for par in parser_list]
    [wr.join()  for wr  in writer_list]    
