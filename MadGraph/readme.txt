Readme about all commands to generate MadGraph events:
1) Move to the MadGraph dir
2) Execute MadGraph5_aMC
     bin/mg5_aMC
3) Generate diagrams
     generate p p > j j w+{0} w- QCD=0 QED<=4, w+ > mu+ vm, w- > j j
4) Generate output directory for the process
     output wp0wm_muvmjj
     quit
5) Prepare the cards for the events generation
     Cards/run_card.dat   --->  nevents, ren. scale, cuts
     Cards/param_card.dat --->  physical parameters
6) Renormalization and factorization for the PDFs is defined in
     SubProcesses/setscales.f
7) Generate events
     bin/generate_events     
