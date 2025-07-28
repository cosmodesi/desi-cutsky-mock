import glob
import numpy as np
import h5py
import time

def generate_shell(file_):
	print(file_)
	
	f = h5py.File(file_, 'r+')
	data = f['galaxy']
			
	start = time.time()

	if "STATUS" in data.keys():
		print("STATUS EXISTS")	
		data["STATUS"][:] = np.zeros(len(data["STATUS"]))

	else:
		print("ERROR: STATUS DOES NOT EXIST")

	if "RAN_NUM_0_1" in data.keys():
		print("RAN_NUM_0_1 EXISTS")
		data['RAN_NUM_0_1'][:] = np.zeros(len(data["STATUS"]))
	else:
		print("ERROR: RAN_NUM_0_1 DOES NOT EXIST")

	f.close()
	print("TIME: It took {} seconds to insert the column the bits.".format(time.time()-start), flush=True)

    
def main():
    folders = glob.glob("/global/cscratch1/sd/avariu/desi/FirstGenMocks/AbacusSummit/CubicBox/ELG/ELG_ran_S*/")
    print(len(folders))

    for folder_ in folders:
        print(folder_)
        files = glob.glob(folder_ + "/*h5py")
        print(len(files))
        for file_ in files:
            generate_shell(file_)

if __name__== '__main__':
    main()