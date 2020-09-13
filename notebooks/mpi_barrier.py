import sys
import time
sys.path.insert(0, '../build/lib')
import intelqs_py as iqs
from intelqs_py import MPIEnvironment


if __name__ == '__main__':
    iqs.init()
    rank = MPIEnvironment.GetRank()
    if not rank:
        print(flush=True)
    print('Process {} before the barrier'.format(rank), flush=True)
    if not rank:
        time.sleep(1)
        print(flush=True)
    MPIEnvironment.Barrier()
    print('Process {} after the barrier'.format(rank), flush=True)
    iqs.finalize()
