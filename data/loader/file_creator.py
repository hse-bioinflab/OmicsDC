import argparse,subprocess, time, urllib, os
import urllib.request

import multiprocessing as mp
from multiprocessing import Queue, Process, current_process,active_children, Manager, shared_memory

from pathlib import Path
from tqdm import tqdm
from collections import Counter
from multiprocessing.shared_memory import SharedMemory

import pandas as pd
import numpy as np

""" New class for checking progress in download"""
class TqdmUpTo(tqdm):
    """Provides `update_to(n)` which uses `tqdm.update(delta_n)`."""
    def update_to(self, b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            self.total = tsize
        return self.update(b * bsize - self.n)  # also sets self.n = b * bsize

"""Block for cmd line arguments"""
cmd_line = argparse.ArgumentParser(description='Script for preparing assembly file')

cmd_line.add_argument(
    '--assembly',
    '-a',
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
    default=6,
    type= int
)
cmd_line.add_argument(
    "--chunk_size",
    '-s',
    type = int,
    default=100000
)
cmd_line.add_argument(
    "--assembly_threshold",
    '-t',
    type=str,
    default='05'
)
"""End of block for cmd line arguments"""  

def writer(que: Queue, done:Queue, file_names: dict):
    """Process to write lines in separate files via numbers in que and id in file_name dict"""
    pause = 0
    is_working = True
    while is_working:
        if not que.empty(): 
            pause = 0   
            SM_meta, df, lines = que.get()
            try:
                ShMem = shared_memory.SharedMemory(name=df)
                list_copy = np.ndarray(SM_meta["shape"], dtype=SM_meta["dtype"], buffer=ShMem.buf)
                try:
                    [open(file_names[line[3]], "+a").write('\t'.join(line)) for line in lines]
                except:
                         pass
            except:
                if SM_meta == "Done":
                    is_working = False
                    done.put("Done")
                pass
                       #print(f"excepted {to_write[3]} ")
                
            ShMem.close()
            done.put(df)            
        else:
            pause +=1
            time.sleep(0.1)
            if pause > 100:
                print(f"{current_process().name} paused to {pause * 0.1}s")

    print(f"{current_process().name} stopped")
    

def worker_file_creator(df,tasks, path, check:bool, file_dict: dict):
    """Ortem`s subprocess"""
    while not tasks.empty():
        start, end = tasks.get()
        for f in range(start, end):
            file = list(df.iloc[f])
            actual_path = path / f'{file[0]}_{file[1]}_{file[2]}_{file[3]}_{file[4]}_{file[5]}.bed'
            if check:
                if not Path.is_file(actual_path):
                    file_dict["No_file"] += 1
            else:
                #with open(actual_path, 'w') as f:
                #    f.close()
                    file_dict[file[0]] = actual_path


def ExpListProcessing(data: pd.DataFrame, path: str, n_workers: int, check: bool, file_dict: dict):
    """ Ortem`s functions"""
    tasks = Queue()
    for batch_start in range(0, len(data) - 1000, 1000):
        tasks.put((batch_start, batch_start + 1000))
    tasks.put((len(data) - 1000, len(data)))
    
    procs = [Process(target=worker_file_creator, args=(data, tasks,path, check, file_dict )) for _ in range(n_workers) ]

    [p.start()  for p in procs]
    [p.join()   for p in procs]
    return 0
    
def Put2QandSM_Process(df:pd.DataFrame, lines_bufer: dict, writers_que: list):
    """Function for put and free buffer from __main__ process to child`s Q"""
    Df2List = np.array(df.values.tolist())
    ShMem = shared_memory.SharedMemory(create=True, size=Df2List.nbytes)
    SM_list = {"shape":Df2List.shape, "dtype": Df2List.dtype}
    l_copy = np.ndarray(Df2List.shape, dtype=Df2List.dtype, buffer=ShMem.buf)
    l_copy[:] = Df2List[:]
    pause = 0
    #print("w8")
    while (not all(q.empty() for q in writers_que)):
        pause += 1
        time.sleep(1)
        pass
    for buf in lines_bufer:
        writers_que[buf].put((SM_list,ShMem.name,lines_bufer[buf]))
        #print((SM_list,ShMem.name,lines_bufer[buf]))
    ShMem.close()
    #print(f"this iter was waiting to {pause*0.01}s")

def FreeTheMemmory_Process(DoneQ : Queue, n_writers: int):
    """Process to free shared memory blocks"""
    is_working = True
    done_counter = Counter()
    while is_working:
        while not DoneQ.empty():
            done_counter[done_que.get()] +=1
        for mem_buf in done_counter:
                if done_counter[mem_buf] == n_writers:
                    ShMem = shared_memory.SharedMemory(name=mem_buf)
                    ShMem.unlink()
                    done_counter[mem_buf] = 0
                    #print(f"free {mem_buf}")
        if  done_counter["Done"] != 0:
            for mem_buf in done_counter:
                 if done_counter[mem_buf] != 0:
                    ShMem = shared_memory.SharedMemory(name=mem_buf)
                    ShMem.unlink()
            is_working = False
        time.sleep(1)
    print(f"{current_process().name}(GC) stopped")
        

if __name__ == '__main__':
    args = cmd_line.parse_args()
    NWORKERS = args.nworkers
    NWRITERS = 20
    file_dick = Manager().dict()
    WORKING_DIR = Path(f"./resources/{args.assembly}")
    chunk_size = 500000
    

    if not Path('./resources/experimentList.tab').exists(): # if user has no ungz list tab - ungz it
        proc = subprocess.Popen(["gunzip",Path('./resources/experimentList.tab.gz')])
        proc.wait()

    # read ExpList df
    ExpList = pd.read_csv(Path("./resources/experimentList.tab"),
                        sep = ',', 
                        #usecols=range(5),
                        header = None
                        )
    ExpList = ExpList[ExpList[1] == "mm9"].replace(
        [' ','/'], '_',
        inplace=False, regex=True
        )
    

    if WORKING_DIR.exists():
        #os.mkdir(WORKING_DIR) 
        """This condition for first setup with download assembly file from CA DB"""
        #urllib = getattr(urllib, 'request', urllib)
        Path2File = WORKING_DIR / f"allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed.gz"
        # eg_link = f"https://chip-atlas.dbcls.jp/data/{args.assembly}/allPeaks_light/allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed.gz"
        # with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1,
        #             desc=eg_link.split('/')[-1]) as t:  # all optional kwargs
        #     urllib.urlretrieve(eg_link, filename= Path2File,
        #                     reporthook=t.update_to, data=None)
        #     t.total = t.n

        #proc = subprocess.Popen(["gunzip",Path2File])     

        """Create files for each id experiment"""
        ExpListProcessing(ExpList, n_workers=12,path = WORKING_DIR, check = False, file_dict = file_dick)

        #proc.wait()
        Path2File = WORKING_DIR / f"allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed"
        LinesInDock = str(subprocess.check_output(['wc','-l',Path2File]))
        LinesInDock = int(LinesInDock[2:LinesInDock.find(' ')])

        print(LinesInDock)

        iter_df = pd.read_csv(
                Path2File,
                chunksize=chunk_size,
                sep = '\t',
                header= None
            )

        lines_buf = {key:[] for key in range(NWRITERS)}# buffer for each chunk
        writer_list = []# buffer for writers

        """Creating queues"""
        workers_ques = [ Queue() for i in range(NWRITERS) ]
        done_que = Queue()
        
        # setup and run checker to free memory
        Checker = Process(target = FreeTheMemmory_Process, args = (done_que,NWRITERS,))
        Checker.start()

        #print("exala")
        for chunk in tqdm(iter_df, total = int(LinesInDock/chunk_size)):
            """W8 till each process take his part from writers_Q"""
            for index,id_exp in enumerate(chunk[3]):
                """ Putting lines for each writer to work in buffer for  writersQ"""
                q_num = int(id_exp[3:])%NWRITERS
                lines_buf[q_num].append(index)
            if not writer_list:
                 # Setup and run writers        
                for i in range(NWRITERS):
                    p = Process(target=writer, args=(workers_ques[i], done_que, file_dick,))
                    writer_list.append(p)
                    p.start()
            """without sending package to writerQ, writer cannot take smth from chunkQ"""
            Put2QandSM_Process(chunk,lines_buf,workers_ques)
            lines_buf = {key:[] for key in range(NWRITERS)}# refresh buffer
           
        
        [q.put(("Done",0,0)) for q in workers_ques]
        [wr.join()  for wr  in writer_list]
        Checker.join()
        #subprocess.run(f"rm ./resources/allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed.gz")
    
    elif args.check:
        file_dick["No_file"] = 0 
        result = ExpListProcessing(df, n_workers=8, check = False, file_dict = file_dick)
        print(file_dick["No_file"])
    
    else:
        print("Already satisfaied")    