import sys
import time
import numpy as np
import pandas as pd
import uproot
from lhe import save_tree
import multiprocessing

print("Starting ", end='')
print(time.ctime())
start = time.time()
argc = len(sys.argv)
if argc < 4:
	print("Error: Missing arguments: [filename] [part_type] [label] [Nchunks]")
	#print("particle types: 'mu', 'l', 'j', 'w', 'final'")
	print("particle types: 'l', 'j', 'final'")
	print("particle types: 'final'")
	sys.exit(1)

for arg in sys.argv:
	print(arg, end=" ")
print("")

filename = sys.argv[1]
part_type = sys.argv[2]
label = sys.argv [3]
Nchunks = int(sys.argv[4])

#varset = ['Particle.PT', 'Particle.Eta', 'Particle.Phi', 'Particle.PID', 'Event.ScalePDF']
varset = ['Particle.PT', 'Particle.Eta', 'Particle.Phi', 'Particle.PID']
print("varset: ", end='')
print(varset)

file = uproot.open(filename)
tree = file[b'Delphes;1/Particle']
Nentries = len(tree)
path = 'dataframes/' + label + '/'

def df_chunk(entrystart, entrystop):
	return save_tree (filename, varset, part_type, save=False, label=label, path=path, entrystart=entrystart, entrystop=entrystop)

delimiters = np.linspace(0, Nentries, Nchunks + 1).astype(int)
chunks = [(delimiters[i], delimiters[i+1]) for i in range(len(delimiters[:-1]))]
starttime = time.time()
pool = multiprocessing.Pool()
df_list = pool.starmap(df_chunk, chunks)
pool.close()
df = pd.concat(df_list, ignore_index=True)
out_file = path + 'df_' + part_type + '_' + label + '.csv'
for dataframe in df_list:
	print(dataframe.index[0], end=", ")
print("")
print("Saving dataframe into " + out_file)
df.to_csv(out_file)
print("Finishing ", end='')
print(time.ctime())
print("Processed tree in %d s" % (time.time()-start))
