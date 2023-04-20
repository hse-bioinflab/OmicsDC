import pandas as pd
import numpy as np
import multiprocessing as mp
from multiprocessing import Queue,shared_memory,Process, current_process,active_children
from tqdm import tqdm
import os
import queue
import time

table = pd.read_csv("experimentList.tab",
                    sep = '\t', 
                    usecols=range(6),
                    header = None
                    )

df = table[table[1] == 'hg38']
df = df.replace(' ', '_', regex=True)
df = df.replace('/', '_', regex=True)

exp = table.loc[table[1] == "hg38"]
exp = np.array(exp.values.tolist())
shm = shared_memory.SharedMemory(create=True, size=exp.nbytes)
exp_copy = np.ndarray(exp.shape, dtype=exp.dtype, buffer=shm.buf)
exp_copy[:] = exp[:]


existing_shm = shared_memory.SharedMemory(name=shm.name)
c = np.ndarray(exp.shape, dtype=exp.dtype, buffer=existing_shm.buf)

file_dict = Manager().dict()

iterator = pd.read_csv(
        f'M:/Fast_Work/allPeaks_light.mm9.50.bed',
        chunksize=1_000,
        sep = '\t',
        header= None
    )

def ParseLine(line, que: list, num_writers:int ):
    """Function to separate lines each into special writer"""
    que_num = int(line[3][3:])%num_writers
    while que[que_num].qsize() > 10:
        pass
    que[que_num].put(line)
    
def SendPackage(que: list, parser_bufer: dict):
    skipped = 0
    for buf in parser_bufer:
        if len(parser_bufer[buf]) > 100:
            if (que[buf].qsize() > 1000):
                skipped += 1
            que[buf].put(parser_bufer[buf])
            print(f"send {buf}")
            parser_bufer[buf] = []
    return skipped


def parser(chunk_que, que: dict, num_writers: int, shared_memory_name: str, SM_meta: dict):
    """Process take chunk from chunk_que and line by line putting it into ParseLine funvtion"""
    parser_bufer = {key: [] for key in range(num_writers)}
    buf_trigger = 0
    pause = {key: 0 for key in range(num_writers)}
    existing_shm = shared_memory.SharedMemory(name=shared_memory_name)
    list_copy = np.ndarray(SM_meta["shape"], dtype=SM_meta["dtype"], buffer=existing_shm.buf) # Это нужно для файла с именем
    while not chunk_que.empty():
        data = chunk_que.get()
        for i in range(data.shape[0]):
            line = data.iloc[i,:]
            q_num = int(line[3][3:])%num_writers
            parser_bufer[q_num].append(line)
            buf_trigger += 1
            if buf_trigger > 1000:
                buf_trigger = 0
                pause[q_num] += SendPackage(que, parser_bufer)
                if pause[q_num] > 100:
                    print(f"stacked at {pause[q_num]}")
            ParseLine(data.iloc[i,:], que, num_writers)

def writer(que):
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
            time.sleep(0.01)
            pause += 1

    print(f"{current_process().name} stopped")
    

def worker_file_creator(tasks):
    while not tasks.empty():
        start, end = tasks.get()
        for f in range(start, end):
            file = list(df.iloc[f])
            with open(f'./tmp/{file[0]}_{file[1]}_{file[2]}_{file[3]}_{file[4]}_{file[5]}.bed', 'w+') as f:
                file_dict[file[0]] = f'{file[0]}_{file[1]}_{file[2]}_{file[3]}_{file[4]}_{file[5]}.bed'


def create_files(n_workers):
    tasks = mp.Queue()

    for batch_start in range(0, len(df) - 1000, 1000):
        tasks.put((batch_start, batch_start + 1000))
    tasks.put((len(df) - 1000, len(df)))
    
    procs = [
        mp.Process(target=worker_file_creator, args=(tasks, file_dict,)) for _ in range(n_workers)
    ]

    for p in procs:
        p.start()
    
    for p in procs:
        p.join()




if __name__ == '__main__':
    create_files(n_workers=8)

    parser_list = []
    writer_list = []
    NWRITERS = 3
    NPARSERS = 3
    SM_list = {"shape":exp.shape, "dtype": exp.dtype}
    que_list = [ Queue() for i in range(NWRITERS) ]
    que_chunk = Queue()

    
    
    
    list_chunks = []
    parser_list = []
    verbose = 0

    for chunk in iterator:
        que_chunk.put(chunk)
        while que_chunk.qsize() > NPARSERS:
            if len(parser_list) == 0:
                for i in range(NWRITERS):
                    p = Process(target=writer, args=(que_list[i],))
                    writer_list.append(p)
                    p.start()
                for chunk,_ in zip(iterator,range(NPARSERS)):
                    p = Process(target=parser, args=(que_chunk,que_list,NWRITERS,shm.name,SM_list,))
                    parser_list.append(p)
                    p.start()
                
            pass

    
        

    for w in writer_list:
        w.terminate()