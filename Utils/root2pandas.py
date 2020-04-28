import numpy as np
import pandas as pd
import uproot
#from ROOT import TLorentzVector
import matplotlib.pyplot as plt
import os
import sys
import time
from lhe import *
from plotclass import *

print("Starting ", end='')
print(time.ctime())
start = time.time()
argc = len(sys.argv)
if argc < 5:
	print("Error: Missing arguments: [filename] [part_type] [label] [mode] [plot_card]")
	#print("particle types: 'mu', 'l', 'j', 'w', 'final'")
	print("particle types: 'l', 'j', 'final'")
	print("particle types: 'final'")
	sys.exit(1)

filename = sys.argv[1]
part_type = sys.argv[2]
label = sys.argv [3]
mode = sys.argv [4]
plot_card = sys.argv [5]
plot_dir = "plots/" + label + "/"
if not os.path.exists(plot_dir):
	os.makedirs(plot_dir)

#varset = ['Particle.PT', 'Particle.Eta', 'Particle.Phi', 'Particle.PID', 'Event.ScalePDF']
varset = ['Particle.PT', 'Particle.Eta', 'Particle.Phi', 'Particle.PID']
print("varset: ", end='')
print(varset)

modes = ["save", "read"]
df = None
df_events = None
if mode == "save":
	df = save_tree(filename, varset, part_type, save=True, label=label, path='dataframes/' + label + '/')
if mode == "read":
	if label == "wpwm0_jjmuvm":
		dtype = {'PT_l' : object, 'Eta_l' : object, 'Phi_l' : object, 'PID_l' : object, 'mll' : float, 'PT_miss' : float, 'mww' : float,
				 'PT_j' : object, 'Eta_j' : object, 'Phi_j' : object, 'PID_j' : object, 'mjj' : object, 'DeltaEta_jj' : object, 'vbs_tag' : object}
		df = pd.read_csv(filename, dtype=dtype)
		varset = ['PT_l', 'Eta_l', 'PT_j', 'Eta_j', 'mjj', 'DeltaEta_jj']
		create_lists(df, varset)
	if label in ["zz_mumuee", "zz0_mumuee", "zzT_mumuee"]:
		dtype = {'PT_l' : object, 'Eta_l' : object, 'Phi_l' : object, 'PID_l' : object, 'Mothup_l' : object, 'mll' : object, 'm4l' : float,
				 'PT_j' : object, 'Eta_j' : object, 'Phi_j' : object, 'PID_j' : object, 'mjj' : float, 'DeltaEta_jj' : float,
				 'PT_z' : object, 'Eta_z' : object, 'Phi_z' : object, 'PID_z' : object, 'fUniqueID_z' : object, 'mz' : object, 'mzz' : float, 'Theta_e' : float, 'Theta_mu' : float}
		df = pd.read_csv(filename, usecols=['PT_l', 'Eta_l', 'mll', 'm4l', 'PT_j', 'Eta_j', 'mjj', 'DeltaEta_jj', 'PT_z', 'Eta_z', 'mz', 'mzz', 'Theta_e','Theta_mu'], dtype=dtype)
		#df = pd.read_csv(filename, dtype=dtype)
		varset = ['PT_l', 'Eta_l', 'mll', 'PT_j', 'Eta_j']
		create_lists(df, varset)

# Saving plot_card configuration
with open(plot_card, 'r') as file:
	first = True
	for line in file:
		if line.startswith('#'):
			continue
		line = line.split('!')[0]
		varname = line.split('=')[1].strip()
		var = line.split('=')[0].strip()
		if ',' in var:
			var = [float(x) for x in var.split(',')]
			if "bin" in varname:
				var = np.linspace(var[0], var[1], var[2])
		else:
			if first == False:
				var = float(var)
		
		globals()[varname] = var
		first = False

if "zz" in label:
	fileroot = "trees/unweighted_events_zz_mumuee_10k.root"		# HARDCODED
	print("Opening %s" % fileroot)
	file = uproot.open(fileroot)
	tree_delphes = file[b'Delphes;1']
	df_events = tree_delphes.pandas.df([b'Event.ScalePDF'])
	print(str(tree_delphes.name) + " contains " + str(len(tree_delphes)) + " entries")
	file.close()
	
	PT_l = inclusive(df['PT_l'])
	PT1_l = column_leading(df['PT_l'], 1)
	Eta_l = inclusive(df['Eta_l'])
	PT_j = inclusive(df['PT_j'])
	Eta_j = inclusive(df['Eta_j'])
	mll = inclusive(df['mll'])
	mjj = df['mjj'].values.tolist()
	DeltaEta_jj = df['DeltaEta_jj']
	m4l = df['m4l'].values.tolist()
	mzz = df['mzz'].values.tolist()
	cosTheta_e = np.cos(df['Theta_e'].values).tolist()
	cosTheta_mu = np.cos(df['Theta_mu'].values).tolist()
	df_highmzz = df.query('mzz > 750')
	cosTheta_e_highmzz = np.cos(df_highmzz['Theta_e'].values).tolist()
	cosTheta_mu_highmzz = np.cos(df_highmzz['Theta_mu'].values).tolist()

	pt_leptons = plot("pt_leptons", [PT_l, PT1_l], bin_pt_l, ["leptons", "leading"], cut=[pt_l_min, pt1_l_min], color=['blue', 'red'], type='hist', plot_dir=plot_dir)
	pt_jets = plot("pt_jets", PT_j, bin_pt_j, "jets", cut=pt_j_min, color='blue', type='hist', plot_dir=plot_dir)
	eta_particles = plot("eta_particles", [Eta_l, Eta_j], bin_eta, ["leptons", "jets"], cut=[eta_l_max, eta_j_max], color=['yellow', 'green'], type='hist', plot_dir=plot_dir)
	dilepton = plot("mll", mll, bin_mll, ["", ""], cut=[mmll, mmllmax], color='blue', type='hist', plot_dir=plot_dir)
	dijet = plot("mjj", mjj, bin_mjj, "jets", cut=mjj_min, color='blue', type='hist', plot_dir=plot_dir)
	deltaetajj = plot("deltaetajj", DeltaEta_jj, bin_deta, "jets", cut=deta_jj_min, color='yellow', type='hist', plot_dir=plot_dir)
	fourlepton = plot("m4l", m4l, bin_mzz, "", cut=m4l_min, color='cyan', type='hist', plot_dir=plot_dir)
	diboson = plot("mzz", mzz, bin_mzz, "", cut=m4l_min, color='cyan', type='hist', plot_dir=plot_dir)
	costheta = plot("costheta", [cosTheta_e, cosTheta_mu], bin_ctheta, ["cos$θ_e$", "cos$θ_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir)
	costheta_highmzz = plot("costheta", [cosTheta_e_highmzz, cosTheta_mu_highmzz], bin_ctheta, ["cos$θ_e$", "cos$θ_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='highmzz')
	m4l_vs_mzz = plot("m4l_vs_mzz", [m4l, mzz], [bin_mzz[0], bin_mzz[-1]], type='scatter', plot_dir=plot_dir)

	for graph in [pt_leptons, pt_jets, eta_particles, dilepton, dijet, deltaetajj, fourlepton, diboson, costheta, costheta_highmzz, m4l_vs_mzz]:
		graph.draw()
		del graph

if "wp" in label:
	fileroot = "trees/unweighted_events_wpwm0_jjmuvm_5k.root"
	print("Opening %s" % fileroot)
	file = uproot.open(fileroot)
	tree_delphes = file[b'Delphes;1']
	df_events = tree_delphes.pandas.df([b'Event.ScalePDF'])
	print(str(tree_delphes.name) + " contains " + str(len(tree_delphes)) + " entries")

	PT_mu = column_leading(df['PT_l'], 1)
	PT_vm = column_leading(df['PT_l'], 2)
	Eta_mu = column_leading(df['Eta_l'], 1)
	PT_miss = df['PT_miss'].values.tolist()
	PT_j = inclusive(df['PT_j'])
	Eta_j = inclusive(df['Eta_j'])
	mjj = inclusive(df['mjj'])
	mjj_vbs = column_leading(df['mjj'], 1)
	mjj_w = column_leading(df['mjj'], 2)
	DeltaEta_jj = inclusive(df['DeltaEta_jj'])
	DeltaEta_jj_vbs = column_leading(df['DeltaEta_jj'], 1)
	DeltaEta_jj_w = column_leading(df['DeltaEta_jj'], 2)
	mww = df['mww'].values.tolist()

	pt_leptons = plot("pt_leptons", [PT_mu, PT_vm], bin_pt_l, ["muon", "neutrino"], cut=[pt_l_min, pt1_l_min], color=['blue', 'red'], type='hist', plot_dir=plot_dir)
	pt_jets = plot("pt_jets", PT_j, bin_pt_j, "jets", cut=pt_j_min, color='blue', type='hist', plot_dir=plot_dir)
	eta_particles = plot("eta_particles", [Eta_mu, Eta_j], bin_eta, ["leptons", "jets"], cut=[eta_l_max, eta_j_max], color=['yellow', 'green'], type='hist', plot_dir=plot_dir)
	dijet = plot("mjj", mjj, bin_mjj, "jets", cut=mjj_min, color='blue', type='hist', plot_dir=plot_dir)
	dijet_vbs = plot("mjj", [mjj_vbs, mjj_w], bin_mjj, ["forward jets", "central jets"], cut=mjj_min, color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='vbs')
	deltaetajj = plot("deltaetajj", DeltaEta_jj, bin_deta, "jets", cut=deta_jj_min, color='yellow', type='hist', plot_dir=plot_dir)
	deltaetajj_vbs = plot("deltaetajj", [DeltaEta_jj_vbs, DeltaEta_jj_w], bin_deta, "jets", cut=deta_jj_min, color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='vbs')
	diboson = plot("mww", mww, bin_mww, "", cut=m4l_min, color='cyan', type='hist', plot_dir=plot_dir)

	for graph in [pt_leptons, pt_jets, eta_particles, dijet, dijet_vbs, deltaetajj, deltaetajj_vbs, diboson]:
		graph.draw()
		del graph

	cov = np.correlate(mww, df_events['Event.ScalePDF'])
	corrcoeff = cov/(np.array(mww).std()*np.array(df_events['Event.ScalePDF']).std())
	print("Correlation between Mww and ScalePDF = %f" % corrcoeff)

end = time.time()
print("Finishing ", end='')
print(time.ctime())
print("Processed tree in %d s" % (end-start))

"""
	plt.figure(figsize=[6, 6])
	plt.scatter((1./np.sqrt(2))*np.array(mww), df_events['Event.ScalePDF'], s=1, color="blue")
	plt.xlim(0, 1650)
	plt.ylim(0, 1650)
	plt.title("Correlation between $M_{WW}/\sqrt{2}$ and ScalePDF")
	plt.xlabel("$M_{WW}/\sqrt{2}$ [GeV]")
	plt.ylabel("ScalePDF [GeV]")
	plt.savefig(plot_dir + "scalepdf.png", format="png")
	print("Saving " + plot_dir + "scalepdf.png")
	plt.show()
	plt.close()
"""
