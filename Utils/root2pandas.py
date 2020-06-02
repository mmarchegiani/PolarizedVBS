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

print("Starting ", end='')
print(time.ctime())
start = time.time()
argc = len(sys.argv)
if argc < 4:
	print("Error: Missing arguments: [filename] [label] [mode] [plot_card]")
	sys.exit(1)

for arg in sys.argv:
	print(arg, end=" ")
print("")

filename = sys.argv[1]
label = sys.argv [2]
mode = sys.argv [3]
plot_card = sys.argv [4]
plot_dir = "plots/" + label + "/"
if not os.path.exists(plot_dir):
	os.makedirs(plot_dir)

#varset = ['Particle.PT', 'Particle.Eta', 'Particle.Phi', 'Particle.PID', 'Event.ScalePDF']
varset = ['Particle.PT', 'Particle.Eta', 'Particle.Phi', 'Particle.PID']
print("varset: ", end='')
print(varset)

modes = ["save", "read"]
df = None
if mode == "save":
	df = save_tree(filename, varset, save=True, label=label, path='dataframes/' + label + '/')
	sys.exit()
if mode == "read":
	if label in ["wpwm0_jjmuvm", "wp0wm_muvmjj"]:
		dtype = {'PT_l' : object, 'Eta_l' : object, 'Phi_l' : object, 'PID_l' : object, 'mll' : float, 'PT_miss' : float, 'mww' : float,
				 'PT_j' : object, 'Eta_j' : object, 'Phi_j' : object, 'PID_j' : object, 'mjj' : object, 'DeltaEta_jj' : object, 'vbs_tag' : object}
		df = pd.read_csv(filename, dtype=dtype)
		varset = ['PT_l', 'Eta_l', 'PT_j', 'Eta_j', 'mjj', 'DeltaEta_jj', 'DeltaPhi_jj', 'DeltaR_jj']
		create_lists(df, varset)
	if label in ["zz_mumuee", "z0z0_mumuee", "zTzT_mumuee", "z0zT_mumuee", "zTz0_mumuee"]:
		dtype = {'PT_l' : object, 'Eta_l' : object, 'Phi_l' : object, 'PID_l' : object, 'Mothup_l' : object, 'mll' : object, 'm4l' : float,
				 'PT_j' : object, 'Eta_j' : object, 'Phi_j' : object, 'PID_j' : object, 'mjj' : float, 'DeltaEta_jj' : float,
				 'PT_z' : object, 'Eta_z' : object, 'Phi_z' : object, 'PID_z' : object, 'fUniqueID_z' : object, 'mz' : object, 'mzz' : float,
				 'Theta_e' : float, 'Theta_mu' : float, 'PT_Ze': float, 'PT_Zmu': float, 'Eta_Ze': float, 'Eta_Zmu': float, 'ScalePDF' : float}
		df = pd.read_csv(filename, usecols=['PT_l', 'Eta_l', 'mll', 'm4l',
											'PT_j', 'Eta_j', 'Phi_j', 'mjj', 'DeltaEta_jj',
											'PT_z', 'Eta_z', 'mz', 'mzz',
											'Theta_e', 'Theta_mu', 'PT_Ze', 'PT_Zmu', 'Eta_Ze', 'Eta_Zmu', 'ScalePDF'], dtype=dtype)
		#df = pd.read_csv(filename, dtype=dtype)
		varset = ['PT_l', 'Eta_l', 'mll', 'PT_j', 'Eta_j', 'Phi_j']
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

if label in ["zz_mumuee", "z0z0_mumuee", "zTzT_mumuee", "z0zT_mumuee", "zTz0_mumuee"]:
	PT_l = inclusive(df['PT_l'])
	PT1_l = column_leading(df['PT_l'], 1)
	Eta_l = inclusive(df['Eta_l'])
	PT_j = inclusive(df['PT_j'])
	Eta_j = inclusive(df['Eta_j'])
	Phi1_j = column_leading(df['Phi_j'], 1)
	Phi2_j = column_leading(df['Phi_j'], 2)
	DeltaPhi_jj = np.minimum(np.abs(Phi1_j - Phi2_j), 2*np.pi - np.abs(Phi1_j - Phi2_j))
	mll = inclusive(df['mll'])
	mjj = df['mjj'].values.tolist()
	DeltaEta_jj = df['DeltaEta_jj']
	m4l = df['m4l'].values.tolist()
	mzz = df['mzz'].values.tolist()
	cosTheta_e = np.cos(df['Theta_e'].values).tolist()
	cosTheta_mu = np.cos(df['Theta_mu'].values).tolist()
	df_highptz = df.query('mzz > 200')
	df_highptz = df_highptz.query('200 < PT_Ze < 400')
	df_highptz = df_highptz.query('200 < PT_Zmu < 400')
	df_highmzz = df.query('mzz > 750')
	cosTheta_e_highmzz = np.cos(df_highmzz['Theta_e'].values).tolist()
	cosTheta_mu_highmzz = np.cos(df_highmzz['Theta_mu'].values).tolist()
	cosTheta_e_highptz = np.cos(df_highptz['Theta_e'].values).tolist()
	cosTheta_mu_highptz = np.cos(df_highptz['Theta_mu'].values).tolist()
	PT_Ze = df['PT_Ze'].values.tolist()
	PT_Zmu = df['PT_Zmu'].values.tolist()
	Eta_Ze = df['Eta_Ze'].values.tolist()
	Eta_Zmu = df['Eta_Zmu'].values.tolist()
	scalepdf = df['ScalePDF'].values.tolist()
	mzz_sqrt2 = ((df['mzz'].values)/np.sqrt(2)).tolist()

	pt_leptons = plot("pt_leptons", [PT_l, PT1_l], bin_pt_l, ["leptons", "leading"], cut=[pt_l_min, pt1_l_min], color=['blue', 'red'], type='hist', plot_dir=plot_dir)
	pt_jets = plot("pt_jets", PT_j, bin_pt_j, "jets", cut=pt_j_min, color='blue', type='hist', plot_dir=plot_dir)
	eta_particles = plot("eta_particles", [Eta_l, Eta_j], bin_eta, ["leptons", "jets"], cut=[eta_l_max, eta_j_max], color=['yellow', 'green'], type='hist', plot_dir=plot_dir)
	deltaphijj = plot("deltaphijj", DeltaPhi_jj, bin_dphi, "jets", color='cyan', type='hist', plot_dir=plot_dir, save=True)
	dilepton = plot("mll", mll, bin_mll, ["", ""], cut=[mmll, mmllmax], color='blue', type='hist', plot_dir=plot_dir)
	dijet = plot("mjj", mjj, bin_mjj, "jets", cut=mjj_min, color='blue', type='hist', plot_dir=plot_dir)
	deltaetajj = plot("deltaetajj", DeltaEta_jj, bin_deta, "jets", cut=deta_jj_min, color='yellow', type='hist', plot_dir=plot_dir)
	fourlepton = plot("m4l", m4l, bin_mzz, "", cut=m4l_min, color='cyan', type='hist', plot_dir=plot_dir)
	diboson = plot("mzz", mzz, bin_mzz, "", cut=m4l_min, color='cyan', type='hist', plot_dir=plot_dir, save=True)
	costheta = plot("costheta", [cosTheta_e, cosTheta_mu], bin_ctheta, ["cos$θ_e$", "cos$θ_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, save=True)
	costheta_highmzz = plot("costheta", [cosTheta_e_highmzz, cosTheta_mu_highmzz], bin_ctheta, ["cos$θ_e$", "cos$θ_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='highmzz', save=True)
	costheta_highptz = plot("costheta", [cosTheta_e_highptz, cosTheta_mu_highptz], bin_ctheta, ["cos$θ_e$", "cos$θ_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, filelabel='highptz', save=True)
	pt_z = plot("pt_z", [PT_Ze, PT_Zmu], bin_pt_j, ["$Z_e$","$Z_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, save=True)
	eta_z = plot("eta_z", [Eta_Ze, Eta_Zmu], bin_eta, ["$Z_e$","$Z_{\mu}$"], color=['blue', 'red'], type='hist', plot_dir=plot_dir, save=True)
	m4l_vs_mzz = plot("m4l_vs_mzz", [m4l, mzz], [bin_mzz[0], bin_mzz[-1]], type='scatter', plot_dir=plot_dir)
	pdfscale_vs_mzz_sqrt2 = plot("pdfscale_vs_mzz_sqrt2", [scalepdf, mzz_sqrt2], [lim_pdf[0], lim_pdf[-1]], type='scatter', plot_dir=plot_dir)

	for graph in [pt_leptons, pt_jets, eta_particles, deltaphijj, dilepton, dijet, deltaetajj, fourlepton, diboson, costheta, costheta_highmzz, costheta_highptz, m4l_vs_mzz, pdfscale_vs_mzz_sqrt2, pt_z, eta_z]:
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
