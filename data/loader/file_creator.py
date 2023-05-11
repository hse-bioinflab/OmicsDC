import argparse,subprocess, time, urllib, os, warnings
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
    default=8,
    type= int
)
cmd_line.add_argument(
    "--chunksize",
    '-s',
    type = int,
    default=500000
)
cmd_line.add_argument(
    "--assembly_threshold",
    '-t',
    type=str,
    default='05'
)
cmd_line.add_argument(
    "--reload",
    '-r',
    type=int,
    default=0
)
cmd_line.add_argument(
    "--download",
    "-d",
    type=int,
    default=0
)
"""End of block for cmd line arguments"""  

def writer(que: Queue, done:Queue, file_names: dict):
    """
    Process to write lines in separate files via numbers in que and id in file_name dict
    """
    pause = 0
    is_working = True
    while is_working:
        if not que.empty():
            # Working cycle 
            pause = 0   
            SM_meta, df, lines = que.get()
            try: 
            # Writing only experiment list tab experiments
                ShMem = shared_memory.SharedMemory(name=df)
                list_copy = np.ndarray(SM_meta["shape"], dtype=SM_meta["dtype"], buffer=ShMem.buf)
                for id,values in lines.items():
                    try: 
                        # Writing only experiment list tab experiments
                        open(file_names[id],'+a').writelines(['\t'.join(line)+'\n' for line in list_copy[[values]][0]])
                    except KeyError:
                        pass 
                ShMem.close()
                done.put(df)   
            except Exception as e:
                if SM_meta == "Done":
                    is_working = False
                    done.put("Done")
                else:
                    raise e
        else:
            pause +=1
            time.sleep(0.1)
            if pause > 100:
                print(f"{current_process().name} paused to {pause * 0.1}s")
    #print(f"{current_process().name} stopped")
    

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
                open(actual_path, 'w+').close()
                file_dict[file[0]] = actual_path


def ExpListProcessing(data: pd.DataFrame, path: Path, n_workers: int, check: bool, file_dict: dict):
    """ Ortem`s functions"""
    tasks = Queue()
    for batch_start in range(0, len(data) - 1000, 1000):
        tasks.put((batch_start, batch_start + 1000))
    tasks.put((len(data) - 1000, len(data)))
    if check:
        file_dict["No_file"] = 0
    
    procs = [Process(target=worker_file_creator, args=(data, tasks,path, check, file_dict )) for _ in range(n_workers) ]

    [p.start()  for p in procs]
    [p.join()   for p in procs]
    return 0
    
def Put2QandSM_Process(df:pd.DataFrame, lines_bufer: list, writers_que: list):
    """
    Function for put and free buffer from __main__ process to child`s Q
    """
    Df2List = df.to_numpy(dtype=str)
    ShMem = shared_memory.SharedMemory(create=True, size=Df2List.nbytes)
    SM_list = {"shape":Df2List.shape, "dtype": Df2List.dtype}
    l_copy = np.ndarray(Df2List.shape, dtype=Df2List.dtype, buffer=ShMem.buf)
    l_copy[:] = Df2List[:]
    pause = 0
    while (not all(q.empty() for q in writers_que)):
        pause += 1
        time.sleep(1)
        pass
    for buf,lines in enumerate(lines_bufer):
        writers_que[buf].put((SM_list,ShMem.name,lines))
    ShMem.close()
    #print(f"this iter was waiting to {pause*0.01}s")

def FreeTheMemmory_Process(DoneQ : Queue, n_writers: int):
    """
    Process to free shared memory blocks
    """
    is_working = True
    done_counter = Counter()
    while is_working:
        while not DoneQ.empty():
            done_counter[done_que.get()] +=1
        for mem_buf in done_counter:
            if done_counter[mem_buf] == n_writers and not ("Done" in mem_buf):
                ShMem = shared_memory.SharedMemory(name=mem_buf)
                ShMem.unlink()
                done_counter[mem_buf] = 0
                #print(f"free {mem_buf}")
            elif done_counter[mem_buf] == n_writers:
                is_working = False
        time.sleep(1)
    print(f"{current_process().name}(GC) stopped")

def DownloadChipAtlasFile(Dir: Path, assembly: str, threshold: str) -> subprocess.Popen:
    """
    Downloads ChipAtlas file for a specific assembly and threshold
    """
    os.mkdir(Dir) 
    import urllib.request,urllib
    urllib = getattr(urllib, 'request', urllib)
    Path2File = Dir / f"allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed.gz"
    eg_link = f"https://chip-atlas.dbcls.jp/data/{args.assembly}/allPeaks_light/allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed.gz"
    with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1,
                  desc=eg_link.split('/')[-1]) as t:  # all optional kwargs
        urllib.urlretrieve(eg_link, filename= Path2File,
                           reporthook=t.update_to, data=None)
        t.total = t.n 
    return subprocess.Popen(["gunzip",Path2File])    

        

if __name__ == '__main__':
    args = cmd_line.parse_args()
    SubporocessHub = None
    NWORKERS = args.nworkers
    if NWORKERS == -1:
        NWORKERS = mp.cpu_count()
    Files = Manager().dict()
    WORKING_DIR = Path(f"./loader/resources/{args.assembly}")
    Path2File = WORKING_DIR / f"allPeaks_light.{args.assembly}.{args.assembly_threshold}.bed"
    chunk_size = args.chunksize
    

    if not Path('./loader/resources/experimentList.tab').exists(): # if user has no ungz list tab - ungz it
        SubporocessHub = subprocess.Popen(["gunzip","-k",Path('./loader/resources/experimentList.tab.gz')])
        SubporocessHub.wait()
        SubporocessHub = None

    # read ExpList df
    ExpList = pd.read_csv(Path("./loader/resources/experimentList.tab"),
                        sep = ',', 
                        #usecols=range(5),
                        header = None
                        )
    ExpList = ExpList[ExpList[1] == args.assembly].replace(
        [' ','/'], '_',
        inplace=False, regex=True
        )
    
    if args.reload and args.download:
        SubporocessHub = subprocess.Popen(["rm","-rf", WORKING_DIR])
        SubporocessHub.wait()
        SubporocessHub = None

    if args.download: SubporocessHub = DownloadChipAtlasFile(WORKING_DIR,args.assembly, args.assembly_threshold)

    if (not WORKING_DIR.exists()):
        print("Need to create working dir and load files")
        exit()
   
    if args.reload or args.check:
        """Create files for each id experiment"""
        ExpListProcessing(ExpList, n_workers=NWORKERS,path = WORKING_DIR, check = args.check, file_dict = Files)
    
    if args.check and Files["No_file"]:
        print("Need reload")
        exit()
    
        
    if SubporocessHub: SubporocessHub.wait()

    if args.reload:
        SubporocessHub = subprocess.Popen(['wc','-l',Path2File], stdout= subprocess.PIPE)
        if args.check and not Path2File.exists():
            print(f"No {Path2File} file exist. Place it there or download")
            exit()
        iter_df = pd.read_csv(
                Path2File,
                chunksize=chunk_size,
                sep = '\t',
                header= None
            )

        lines_buf = [{} for key in range(NWORKERS)]# buffer for each chunk
        writer_list = []# buffer for writers

        """Creating queues"""
        workers_ques = [ Queue() for i in range(NWORKERS) ]
        done_que = Queue()
        
        # setup and run checker to free memory
        Checker = Process(target = FreeTheMemmory_Process, args = (done_que,NWORKERS,))
        Checker.start()

 
        SubporocessHub.wait()
        LinesInDock = str(SubporocessHub.communicate())
        LinesInDock = int(LinesInDock[3:LinesInDock.find(' ')])
        print("Cashing files\n")
        for chunk in tqdm(iter_df, total = int(LinesInDock/chunk_size)):
            """W8 till each process take his part from writers_Q"""
            for index,id_exp in enumerate(chunk[3]):
                """ Putting lines for each writer to work in buffer for  writersQ"""
                q_num = int(id_exp[3:])%NWORKERS
                try:
                    lines_buf[q_num][id_exp] += index,
                except KeyError:
                    lines_buf[q_num][id_exp] = index,
            if not writer_list:
                 # Setup and run writers        
                for i in range(NWORKERS):
                    p = Process(target=writer, args=(workers_ques[i], done_que, Files,))
                    writer_list.append(p)
                    p.start()
            """without sending package to writerQ, writer cannot take smth from chunkQ"""
            Put2QandSM_Process(chunk,lines_buf,workers_ques)
            lines_buf = [{} for key in range(NWORKERS)]# refresh buffer

           
        time.sleep(2)
        [q.put(("Done",0,0)) for q in workers_ques]
        [wr.join()  for wr  in writer_list]
        Checker.join()
        subprocess.Popen(["rm", Path2File])
        subprocess.Popen([])
    if args.check:
        print("Check done. Ready to go!")
    exit()