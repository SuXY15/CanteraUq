from utils import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for m,name,mech in mechs:
    for props['phi'],props['P'],props['T'] in conditions:
        file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                UF,props['phi'],props['P'],props['T'],samplingSize)
        if checkexists(file_name,delete=False):
            print("File %s exists. changing..."%(file_name))
            data_arr = json_reader(file_name)
            os.remove(file_name)
            for data_dict in data_arr:
                data_dict['props']['idt']/=1000
                json_writer(file_name,data_dict)
        else:
            print("File %s not found"%(file_name))

