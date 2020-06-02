import numpy as np
import pandas as pd
import uproot
#from ROOT import TLorentzVector
import matplotlib.pyplot as plt
import os
import sys
import time
from lib.lhe import *
from lib.plotclass import *
import multiprocessing

print("Starting ", end='')
print(time.ctime())
start = time.time()
argc = len(sys.argv)
if argc < 5:
	print("Error: Missing arguments: [filename] [label] [plot_card] [Nchunks] [MC_generator]")
	print("labels: zz_mumuee, z0z0_mumuee, zTzT_mumuee, z0zT_mumuee, zTz0_mumuee")
	sys.exit(1)

for arg in sys.argv:
	print(arg, end=" ")
print("")

filename = sys.argv[1]
label = sys.argv [2]
plot_card = sys.argv [3]
Nchunks = int(sys.argv [4])
MC_generator = sys.argv[5].lower()
home = os.path.expanduser('~')
phantom = False
if MC_generator == "phantom":
	phantom = True
plot_dir = home + "/plots/" + label + "/"
if phantom == True:
	plot_dir = plot_dir + "Phantom/"
else:
	if MC_generator == "madgraph":
		plot_dir = plot_dir + "Madgraph/"
	elif MC_generator == "delphes":
		plot_dir = plot_dir + "Delphes/"

if not os.path.exists(plot_dir):
	os.makedirs(plot_dir)

df = pd.DataFrame()
df_list = []
varset = None

if label == "wpwm0_jjmuvm":
	dtype = {'PT_l' : object, 'Eta_l' : object, 'Phi_l' : object, 'PID_l' : object, 'mll' : float, 'PT_miss' : float, 'mww' : float,
			 'PT_j' : object, 'Eta_j' : object, 'Phi_j' : object, 'PID_j' : object, 'mjj' : object, 'DeltaEta_jj' : object, 'vbs_tag' : object}
	df = pd.read_csv(filename, dtype=dtype)
	varset = ['PT_l', 'Eta_l', 'PT_j', 'Eta_j', 'mjj', 'DeltaEta_jj', 'DeltaPhi_jj', 'DeltaR_jj']
if label in ["zz_mumuee", "z0z0_mumuee", "zTzT_mumuee", "z0zT_mumuee", "zTz0_mumuee"]:
	dtype = {'PT_l' : object, 'Eta_l' : object, 'Phi_l' : object, 'PID_l' : object, 'Mothup_l' : object, 'mll' : object, 'm4l' : float,
			 'PT_j' : object, 'Eta_j' : object, 'Phi_j' : object, 'PID_j' : object, 'mjj' : float, 'DeltaEta_jj' : float,
			 'PT_z' : object, 'Eta_z' : object, 'Phi_z' : object, 'PID_z' : object, 'fUniqueID_z' : object, 'mz' : object, 'mzz' : float, 'PT_ZZ' : float,
			 'Theta_e' : float, 'Theta_mu' : float, 'PT_Ze': float, 'PT_Zmu': float, 'Eta_Ze': float, 'Eta_Zmu': float, 'ScalePDF' : float}
	# Reading of multiple files and concat in an unique dataframe (e.g. filename = "dataframe0*" for 'dataframe01.csv', 'dataframe02.csv'...)
	if '*' in filename:
		csv_file = filename.split('/')[-1]
		data_dir = filename.strip(csv_file)
		for file in os.listdir(data_dir):
			if csv_file.strip('*') in file:
				df = pd.read_csv(data_dir + file, usecols=['PT_l', 'Eta_l', 'mll', 'm4l',
														   'PT_j', 'Eta_j', 'Phi_j', 'mjj', 'DeltaEta_jj',
														   'PT_z', 'Eta_z', 'Phi_z', 'mz', 'mzz', 'PT_ZZ',
														   'Theta_e', 'Theta_mu', 'PT_Ze', 'PT_Zmu', 'Eta_Ze', 'Eta_Zmu', 'ScalePDF'], dtype=dtype)
				df_list.append(df)
		df = pd.concat(df_list, ignore_index=True)
		df_list = []
	else:
		df = pd.read_csv(filename, usecols=['PT_l', 'Eta_l', 'mll', 'm4l',
											'PT_j', 'Eta_j', 'Phi_j', 'mjj', 'DeltaEta_jj',
											'PT_z', 'Eta_z', 'Phi_z', 'mz', 'mzz', 'PT_ZZ',
											'Theta_e', 'Theta_mu', 'PT_Ze', 'PT_Zmu', 'Eta_Ze', 'Eta_Zmu', 'ScalePDF'], dtype=dtype)
	varset = ['PT_l', 'Eta_l', 'mll', 'PT_j', 'Eta_j', 'Phi_j', 'Phi_z']

Nentries = df.shape[0]
print(df.shape)
print("Nentries = %d" % Nentries)
def df_chunk(entrystart, entrystop):
	return create_lists(df, varset, entrystart, entrystop)

# Create the lists for each chunk of the dataframe launching a subprocess for each chunk
delimiters = np.linspace(0, Nentries, Nchunks + 1).astype(int)
chunks = [(delimiters[i], delimiters[i+1]) for i in range(len(delimiters[:-1]))]
print(chunks)
starttime = time.time()
pool = multiprocessing.Pool()
df_list = pool.starmap(df_chunk, chunks)
pool.close()

# Concat all the dataframe chunks
df = pd.concat(df_list, ignore_index=True)

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

if label in ["zz_mumuee", "z0z0_mumuee", "zTzT_mumuee", "z0zT_mumuee", "zTz0_mumuee"]:
	PT_l = inclusive(df['PT_l'])
	PT1_l = column_leading(df['PT_l'], 1)
	PT2_l = column_leading(df['PT_l'], 2)
	PT3_l = column_leading(df['PT_l'], 3)
	PT4_l = column_leading(df['PT_l'], 4)
	Eta_l = inclusive(df['Eta_l'])
	Eta_l_sorted = sort_column(df['Eta_l'], key=abs)
	Eta1_l = column_leading(Eta_l_sorted, 1)
	Eta2_l = column_leading(Eta_l_sorted, 2)
	Eta3_l = column_leading(Eta_l_sorted, 3)
	Eta4_l = column_leading(Eta_l_sorted, 4)
	PT_j = inclusive(df['PT_j'])
	Eta_j = inclusive(df['Eta_j'])
	Phi1_j = column_leading(df['Phi_j'], 1)
	Phi2_j = column_leading(df['Phi_j'], 2)
	DeltaPhi_jj = np.minimum(np.abs(Phi1_j - Phi2_j), 2*np.pi - np.abs(Phi1_j - Phi2_j))
	mll = inclusive(df['mll'])
	mjj = df['mjj'].values.tolist()
	DeltaEta_jj = df['DeltaEta_jj']
	m4l = df['m4l'].values.tolist()
	Phi1_z = column_leading(df['Phi_z'], 1)
	Phi2_z = column_leading(df['Phi_z'], 2)
	DeltaPhi_zz = np.minimum(np.abs(Phi1_z - Phi2_z), 2*np.pi - np.abs(Phi1_z - Phi2_z))
	mzz = df['mzz'].values.tolist()
	cosTheta_e = np.cos(df['Theta_e'].values).tolist()
	cosTheta_mu = np.cos(df['Theta_mu'].values).tolist()
	df_highptz_e = df.query('mzz > 200')
	df_highptz_e = df_highptz_e.query('200 < PT_Ze < 400')
	df_highptz_mu = df.query('mzz > 200')
	df_highptz_mu = df_highptz_mu.query('200 < PT_Zmu < 400')
	df_highmzz = df.query('mzz > 750')
	cosTheta_e_highmzz = np.cos(df_highmzz['Theta_e'].values).tolist()
	cosTheta_mu_highmzz = np.cos(df_highmzz['Theta_mu'].values).tolist()
	cosTheta_e_highptz = np.cos(df_highptz_e['Theta_e'].values).tolist()
	cosTheta_mu_highptz = np.cos(df_highptz_mu['Theta_mu'].values).tolist()
	PT_Ze = df['PT_Ze'].values.tolist()
	PT_Zmu = df['PT_Zmu'].values.tolist()
	PT_ZZ = df['PT_ZZ'].values.tolist()
	df_lowmzz = df.query('mzz < 250')
	df_midmzz = df.query('250 < mzz < 350')
	df_highmzz = df.query('mzz > 350')
	Eta1_l_lowmzz = column_leading(df_lowmzz['Eta_l'], 1)
	Eta2_l_lowmzz = column_leading(df_lowmzz['Eta_l'], 2)
	Eta3_l_lowmzz = column_leading(df_lowmzz['Eta_l'], 3)
	Eta4_l_lowmzz = column_leading(df_lowmzz['Eta_l'], 4)
	Eta1_l_midmzz = column_leading(df_midmzz['Eta_l'], 1)
	Eta2_l_midmzz = column_leading(df_midmzz['Eta_l'], 2)
	Eta3_l_midmzz = column_leading(df_midmzz['Eta_l'], 3)
	Eta4_l_midmzz = column_leading(df_midmzz['Eta_l'], 4)
	Eta1_l_highmzz = column_leading(df_highmzz['Eta_l'], 1)
	Eta2_l_highmzz = column_leading(df_highmzz['Eta_l'], 2)
	Eta3_l_highmzz = column_leading(df_highmzz['Eta_l'], 3)
	Eta4_l_highmzz = column_leading(df_highmzz['Eta_l'], 4)
	Eta_Ze = df['Eta_Ze'].values.tolist()
	Eta_Zmu = df['Eta_Zmu'].values.tolist()
	scalepdf = df['ScalePDF'].values.tolist()
	mzz_sqrt2 = ((df['mzz'].values)/np.sqrt(2)).tolist()

	pt_leptons_cuts = plot("pt_leptons", [PT_l, PT1_l], bin_pt_l, ["leptons", "leading"], cut=[pt_l_min, pt1_l_min], color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='cuts')
	pt_leptons = plot("pt_leptons", [PT1_l, PT2_l, PT3_l, PT4_l], bin_pt_l, ["$p_{T,1}^{\ell}$", "$p_{T,2}^{\ell}$", "$p_{T,3}^{\ell}$", "$p_{T,4}^{\ell}$"], color=['red', 'green', 'cyan', 'blue'], type='histstep', plot_dir=plot_dir, save=True)
	pt_jets = plot("pt_jets", PT_j, bin_pt_j, "jets", cut=pt_j_min, color='blue', type='hist', plot_dir=plot_dir, save=True)
	eta_particles = plot("eta_particles", [Eta_l, Eta_j], bin_eta, ["leptons", "jets"], cut=[eta_l_max, eta_j_max], color=['yellow', 'green'], type='hist', plot_dir=plot_dir)
	eta_leptons = plot("eta_leptons", [Eta1_l, Eta2_l, Eta3_l, Eta4_l], bin_eta, ["$\eta_1^{\ell}$", "$\eta_2^{\ell}$", "$\eta_3^{\ell}$", "$\eta_4^{\ell}$"], color=['blue', 'cyan', 'green', 'red'], type='histstep', plot_dir=plot_dir, save=True)
	eta_leptons_lowmzz = plot("eta_leptons", [Eta1_l_lowmzz, Eta2_l_lowmzz, Eta3_l_lowmzz, Eta4_l_lowmzz], bin_eta, ["$\eta_1^{\ell}$", "$\eta_2^{\ell}$", "$\eta_3^{\ell}$", "$\eta_4^{\ell}$"], color=['blue', 'cyan', 'green', 'red'], type='histstep', plot_dir=plot_dir, filelabel='lowmzz', save=True)
	eta_leptons_midmzz = plot("eta_leptons", [Eta1_l_midmzz, Eta2_l_midmzz, Eta3_l_midmzz, Eta4_l_midmzz], bin_eta, ["$\eta_1^{\ell}$", "$\eta_2^{\ell}$", "$\eta_3^{\ell}$", "$\eta_4^{\ell}$"], color=['blue', 'cyan', 'green', 'red'], type='histstep', plot_dir=plot_dir, filelabel='midmzz', save=True)
	eta_leptons_highmzz = plot("eta_leptons", [Eta1_l_highmzz, Eta2_l_highmzz, Eta3_l_highmzz, Eta4_l_highmzz], bin_eta, ["$\eta_1^{\ell}$", "$\eta_2^{\ell}$", "$\eta_3^{\ell}$", "$\eta_4^{\ell}$"], color=['blue', 'cyan', 'green', 'red'], type='histstep', plot_dir=plot_dir, filelabel='highmzz', save=True)
	deltaphijj = plot("deltaphijj", DeltaPhi_jj, bin_dphi, "jets", color='cyan', type='hist', plot_dir=plot_dir, save=True)
	dilepton = plot("mll", mll, bin_mll, ["", ""], cut=[mmll, mmllmax], color='blue', type='hist', plot_dir=plot_dir)
	dijet = plot("mjj", mjj, bin_mjj, "jets", cut=mjj_min, color='blue', type='hist', plot_dir=plot_dir)
	deltaetajj = plot("deltaetajj", DeltaEta_jj, bin_deta, "jets", cut=deta_jj_min, color='yellow', type='hist', plot_dir=plot_dir, save=True)
	fourlepton = plot("m4l", m4l, bin_mzz, "", cut=m4l_min, color='cyan', type='hist', plot_dir=plot_dir)
	diboson = plot("mzz", mzz, bin_mzz, "", cut=m4l_min, color='cyan', type='hist', plot_dir=plot_dir, save=True)
	costheta = plot("costheta", [cosTheta_e, cosTheta_mu], bin_ctheta, ["cos$θ_e$", "cos$θ_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, save=True)
	costheta_highmzz = plot("costheta", [cosTheta_e_highmzz, cosTheta_mu_highmzz], bin_ctheta, ["cos$θ_e$", "cos$θ_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='highmzz', save=True)
	costheta_highptz = plot("costheta", [cosTheta_e_highptz, cosTheta_mu_highptz], bin_ctheta, ["cos$θ_e$", "cos$θ_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='highptz', save=True)
	deltaphizz = plot("deltaphizz", DeltaPhi_zz, bin_dphi, "Z", color='cyan', type='hist', plot_dir=plot_dir, save=True)
	pt_z = plot("pt_z", [PT_Ze, PT_Zmu], bin_pt_j, ["$Z_e$","$Z_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, save=True)
	pt_zz = plot("pt_zz", PT_ZZ, bin_pt_j, color='blue', type='hist', plot_dir=plot_dir, save=True)
	eta_z = plot("eta_z", [Eta_Ze, Eta_Zmu], bin_eta, ["$Z_e$","$Z_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, save=True)
	m4l_vs_mzz = plot("m4l_vs_mzz", [m4l, mzz], [bin_mzz[0], bin_mzz[-1]], type='scatter', plot_dir=plot_dir)
	pdfscale_vs_mzz_sqrt2 = plot("pdfscale_vs_mzz_sqrt2", [scalepdf, mzz_sqrt2], [lim_pdf[0], lim_pdf[-1]], type='scatter', plot_dir=plot_dir)

	for graph in [pt_leptons, pt_jets, eta_particles, eta_leptons, eta_leptons_lowmzz, eta_leptons_midmzz, eta_leptons_highmzz,
				  deltaphijj, dilepton, dijet, deltaetajj, fourlepton, diboson, costheta, costheta_highmzz, costheta_highptz, deltaphizz,
				  pt_z, pt_zz, eta_z, m4l_vs_mzz, pdfscale_vs_mzz_sqrt2]:
		graph.draw()
		del graph

if "wp" in label:
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
	DeltaPhi_jj = inclusive(df['DeltaPhi_jj'])
	DeltaPhi_jj_vbs = column_leading(df['DeltaPhi_jj'], 1)
	DeltaPhi_jj_w = column_leading(df['DeltaPhi_jj'], 2)
	DeltaR_jj = inclusive(df['DeltaR_jj'])
	DeltaR_jj_vbs = column_leading(df['DeltaR_jj'], 1)
	DeltaR_jj_w = column_leading(df['DeltaR_jj'], 2)
	mww = df['mww'].values.tolist()
	mww_sqrt2 = ((df['mww'].values)/np.sqrt(2)).tolist()
	scalepdf = df['ScalePDF'].values.tolist()

	pt_leptons = plot("pt_leptons", [PT_mu, PT_vm], bin_pt_l, ["muon", "neutrino"], cut=[pt_l_min, pt1_l_min], color=['blue', 'red'], type='hist', plot_dir=plot_dir)
	pt_jets = plot("pt_jets", PT_j, bin_pt_j, "jets", cut=pt_j_min, color='blue', type='hist', plot_dir=plot_dir)
	eta_particles = plot("eta_particles", [Eta_mu, Eta_j], bin_eta, ["leptons", "jets"], cut=[eta_l_max, eta_j_max], color=['yellow', 'green'], type='hist', plot_dir=plot_dir)
	dijet = plot("mjj", mjj, bin_mjj, "jets", cut=mjj_min, color='blue', type='hist', plot_dir=plot_dir)
	dijet_vbs = plot("mjj", [mjj_vbs, mjj_w], bin_mjj, label=["forward jets", "central jets"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='vbs')
	deltaetajj = plot("deltaetajj", DeltaEta_jj, bin_deta, "jets", cut=deta_jj_min, color='yellow', type='hist', plot_dir=plot_dir)
	deltaetajj_vbs = plot("deltaetajj", [DeltaEta_jj_vbs, DeltaEta_jj_w], bin_deta, ["forward jets", "central jets"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='vbs')
	deltaphijj = plot("deltaphijj", DeltaPhi_jj, bin_dphi, "jets", color='cyan', type='hist', plot_dir=plot_dir)
	deltaphijj_vbs = plot("deltaphijj", [DeltaPhi_jj_vbs, DeltaPhi_jj_w], bin_dphi, ["forward jets", "central jets"], color=['cyan', 'purple'], type='hist', plot_dir=plot_dir, filelabel='vbs')
	deltarjj = plot("deltarjj", DeltaR_jj, bin_dr, "jets", color='cyan', type='hist', plot_dir=plot_dir)
	deltarjj_vbs = plot("deltarjj", [DeltaR_jj_vbs, DeltaR_jj_w], bin_dr, ["forward jets", "central jets"], color=['cyan', 'purple'], type='hist', plot_dir=plot_dir, filelabel='vbs')
	diboson = plot("mww", mww, bin_mww, "", cut=m4l_min, color='cyan', type='hist', plot_dir=plot_dir)
	pdfscale_vs_mww_sqrt2 = plot("pdfscale_vs_mww_sqrt2", [scalepdf, mww_sqrt2], [lim_pdf[0], lim_pdf[-1]], type='scatter', plot_dir=plot_dir)

	for graph in [pt_leptons, pt_jets, eta_particles, dijet, dijet_vbs, deltaetajj, deltaetajj_vbs, deltaphijj, deltaphijj_vbs, deltarjj, deltarjj_vbs, diboson, pdfscale_vs_mww_sqrt2]:
		graph.draw()
		del graph

	cov = np.correlate(mww, scalepdf)
	corrcoeff = cov/(np.array(mww).std()*np.array(scalepdf).std())
	print("Correlation between Mww and ScalePDF = %f" % corrcoeff)

end = time.time()
print("Finishing ", end='')
print(time.ctime())
print("Processed tree in %d s" % (end-start))

