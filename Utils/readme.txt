Readme about the use of the Python3 script root2pandas.py:
1) Five arguments are required:
python root2pandas.py [filename] [part_type] [label] [mode] [plot_card]
2) Save the MadGraph output TTree into a Pandas dataframe and produce plots
python root2pandas.py zz_events.root final zz_mumuee save cards/plot_card.dat
3) Read the Pandas dataframe and produce plots
python root2pandas.py df_zz_events.csv final zz_mumuee read cards/plot_card.dat
4) lhe.py: functions to read and process the ROOT TTree 
5) plotclass.py: Python class to manage the plots of kinematical variables 
6) Information about binning and cuts for the plots is passed through the input file [plot_card]
