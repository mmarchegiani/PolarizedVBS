import os
import sys
import time
import numpy as np
import pandas as pd
import uproot
from lib.lhe import save_tree
import multiprocessing

print("Starting ", end='')
print(time.ctime())
start = time.time()
argc = len(sys.argv)
if argc < 4:
	print("Error: Missing arguments: [filename] [label] [Nchunks] [MC_generator]")
	sys.exit(1)

for arg in sys.argv:
	print(arg, end=" ")
print("")

filename = sys.argv[1]
label = sys.argv [2]
Nchunks = int(sys.argv[3])
MC_generator = "MadGraph"
phantom = False
if argc > 5:
	MC_generator = sys.argv[4]
	if MC_generator == "phantom":
		phantom = True

home = os.path.expanduser('~')
path = home + '/dataframes/' + label + '/'

varset = ['Particle.PT', 'Particle.Eta', 'Particle.Phi', 'Particle.PID']
print("varset: ", end='')
print(varset)

file = uproot.open(filename)
tree = file[b'Delphes;1/Particle']
Nentries = len(tree)

def df_chunk(entrystart, entrystop):
	return save_tree (filename, varset, save=False, label=label, path=path, entrystart=entrystart, entrystop=entrystop, phantom=phantom)

delimiters = np.linspace(0, Nentries, Nchunks + 1).astype(int)
chunks = [(delimiters[i], delimiters[i+1]) for i in range(len(delimiters[:-1]))]
starttime = time.time()
pool = multiprocessing.Pool()
df_list = pool.starmap(df_chunk, chunks)
pool.close()
df = pd.concat(df_list, ignore_index=True)
out_file = path + 'df_' + label + '.csv'
if phantom == True:
	df['Nevts_gen'] = [Nentries]*df.shape[0]
	df['overall_efficiency'] = [df['cut_efficiency'].mean()]*df.shape[0]
	out_file = out_file.replace('df_', 'df_phantom_')
print("Saving dataframe into " + out_file)
df.to_csv(out_file)
print("Finishing ", end='')
print(time.ctime())
print("Processed tree in %d s" % (time.time()-start))
