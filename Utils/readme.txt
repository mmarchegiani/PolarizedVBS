Readme about the use of the Python3 script root2pandas.py:
1) Four arguments are required:
     python root2pandas.py [filename] [label] [mode] [plot_card]
2) Save the MadGraph output TTree into a Pandas dataframe and produce plots
     python root2pandas.py zz_events.root zz_mumuee save cards/plot_card.dat
3) Read the Pandas dataframe and produce plots
     python root2pandas.py df_zz_events.csv zz_mumuee read cards/plot_card.dat
4) lib/lhe.py: functions to read and process the ROOT TTree 
5) lib/plotclass.py: Python class to manage the plots of kinematical variables 
6) Information about binning and cuts for the plots is passed through the input file [plot_card]

root2pandas_multiproc.py is a Python script to read the .root file exploiting to allow for multiprocessing:
1) Four arguments are required:
     python root2pandas_multiproc.py [filename] [label] [Nchunks] [MC_generator]
2) Save the MadGraph output TTree into a Pandas dataframe using 20 processors:
     python root2pandas_multiproc.py zz_events.root zz_mumuee 20 madgraph

pandas2plots_multiproc.py is a Python script to read the Pandas dataframe and produce plots exploiting multiprocessing:
1) Five arguments are required:
     python pandas2plots_multiproc.py [filename] [label] [plot_card] [Nchunks] [MC_generator]
2) Read the Pandas dataframe and produce plots using 20 processors:
     python pandas2plots_multiproc.py df_zz_events.csv zz_mumuee read cards/plot_card.dat 20
