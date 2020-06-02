import numpy as np
import pandas as pd
import uproot
from ROOT import TLorentzVector, TVector3
#import matplotlib.pyplot as plt
import os
import sys

Mz = 91.188		# Z boson mass (PDG) [GeV]
Mw = 80.419		# W boson mass (PDG) [GeV]
Gf = 1.16639e-5	# Fermi coupling constant (PDG) [GeV^2]

def variable(tree, branchname):		# Inclusive variable from Tree
	var = tree[branchname]
	var_inclusive = np.array(var.array().flatten())
	return var_inclusive

def df_particle(df, part_type):
	if(part_type == "e"):
		df_p = df.loc[np.abs(df['Particle.PID']) == 11]
		print("%d electrons selected" % len(df_p))
	if(part_type == "mu"):
		df_p = df.loc[np.abs(df['Particle.PID']) == 13]
		print("%d muons selected" % len(df_p))
	if(part_type == "l"):
		df_p = df.loc[(np.abs(df['Particle.PID']) >= 11) & (np.abs(df['Particle.PID']) <= 16) ]
		print("%d leptons selected" % len(df_p))
	if(part_type == "j"):
		df_p = df.loc[(np.abs(df['Particle.PID']) >= 1) & (np.abs(df['Particle.PID']) <= 6)]
		df_p = df_p.loc[df_p['Particle.PT'] > 0]       # By requiring PT>0 we remove the unphysical partons with PT=0
		print("%d partons selected" % len(df_p))
	if(part_type == "z"):
		df_p = df.loc[np.abs(df['Particle.PID']) == 23]
		print("%d Z bosons selected" % len(df_p))
	if(part_type == "w"):
		df_p = df.loc[np.abs(df['Particle.PID']) == 24]
		print("%d W bosons selected" % len(df_p))
	return df_p

def argmaxs(l):			# Function which returns the ranked index list
	arg_list = []
	buffer = l.copy()
	for i in range(len(l)):
		i_max = np.argmax(buffer)
		arg_list.append(i_max)
		buffer[i_max] = float('-inf')
	return arg_list

def metric(mxx, M):
	metric_list = []
	for m in mxx:
		metric_list.append(-abs(m - M))		# This choice for the metric as negative is such that the Z,W candidates have the highest metric
	return metric_list

def build_col(df, varset):   		# Leading and sub-leading particles are defined with respect to
	varname = varset[0]      		# the first variable of the list (varset[0]) (normally PT)
	entrystart = df.index[0][0]		# The function generates lists of the type PT=[40.5, 36, 20, 10]
	x = []
	list_var = []
	for (j, var) in enumerate(varset):
		list_var.append([])
	length = len(df.loc[entrystart][varname].values)
	i = 0
	for index, row in df.iterrows():
		x.append(row[varname])
		
		if(i == length - 1):
			i_max = argmaxs(x)      # N.B.: i_max is the list with indices of maxima
			for (j, var) in enumerate(varset):
				list_ordered = []
				for k in range(length):
					list_ordered.append(df.loc[index[0]][var].values[i_max[k]])
				list_var[j].append(list_ordered)
			i = -1
			x = []
		i += 1
	return list_var

def build_dataframe(df, varset):		# Function which builds the dataframe with the ordered lists
	print("Building dataframe")
	list_ordered = build_col(df, varset)
	varlist = []
	for var in varset:
		varlist.append(var.split('.')[-1])   # Remove the Particle. prefix in the column name
	d = { varlist[i] : list_ordered[i] for i in range(0, len(varlist) ) }
	df_leading = pd.DataFrame(data=d)
	return df_leading

def drop_singleZ(df, Nevts):
	print("Dropping single-Z events")
	singleZ_indices=[]
	for index, row in df.iterrows():
		nPart = len(df.loc[index[0]])
		if nPart < 10:		# Condition for missing Z boson in the event
			if index[0] in singleZ_indices:
				continue
			else:
				singleZ_indices.append(index[0])
	df.drop(singleZ_indices, inplace=True)
	efficiency = 1.0 - float(len(singleZ_indices))/float(Nevts)
	print("len(singleZ_indices) = %d" % len(singleZ_indices))
	print("Nevts = %d" % Nevts)
	return singleZ_indices, efficiency

def eemumu(df, mll_ordering='flavor'):
	boolean_list = []
	c = 0
	for index, row in df.iterrows():		# Remember that 'row' is just a copy of the row of the dataframe! Need to access to the element with df.loc[index, 'variable']
		doubleflavor = True
		pid = row['PID']
		first = pid[0]
		pair_1 = [0]
		pair_2 = [0, 1, 2, 3]
		pair_list = []
		l = [0]
		l_bar = []
		for i in range(1,len(pid)):
			if pid[i] == first:
				doubleflavor = False
				boolean_list.append(doubleflavor)
				continue
		if doubleflavor == True:
			c += 1
			boolean_list.append(doubleflavor)
			for i in range(1,len(pid)):
				if (pid[i] == -first):
					pair_1.append(pair_2.pop(i))
					pair_2.pop(0)
			pair_list = [pair_1, pair_2]
			mll = []
			for pair in pair_list:
				p = []
				for i in range(2):
					p.append(TLorentzVector())
					p[i].SetPtEtaPhiM(row['PT'][pair[i]], row['Eta'][pair[i]], row['Phi'][pair[i]], 0.)
				mll.append((p[0]+p[1]).M())
			#z_pole_distance = metric(mll, Mz)
			#i_max = argmaxs(z_pole_distance)
			#df['mll'][index] = [mll[i_max[0]], mll[i_max[1]]]
		else:
			for i in range(1,len(pid)):
				if (pid[i] == first):
					l.append(i)
				else:
					l_bar.append(i)
			pairs = ((x,y) for x in l for y in l_bar)
			pair_list = []
			#pair_with_leading_lepton = []
			mll = []
			for u, v in pairs:
				pair_list.append([u, v])
			
			for pair in pair_list:
				p = []
				leading_lepton = False
				for i in range(2):
					if pair[i] == 0:
						leading_lepton = True
					p.append(TLorentzVector())
					p[i].SetPtEtaPhiM(row['PT'][pair[i]], row['Eta'][pair[i]], row['Phi'][pair[i]], 0.)
				pair_with_leading_lepton.append(leading_lepton)
				mll.append((p[0]+p[1]).M())

		z_pole_distance = metric(mll, Mz)
		i_max = argmaxs(z_pole_distance)
		(u,v) = pair_list[i_max[0]]
		if mll_ordering == "flavor":
			if np.abs(pid[u]) == 11:								# mll ordered by [ee pair, mumu pair]
				df['mll'][index] = [mll[i_max[0]], mll[i_max[1]]]
			else:
				df['mll'][index] = [mll[i_max[1]], mll[i_max[0]]]
		if mll_ordering == "pt":
			if pair_with_leading_lepton[i_max[0]] == True:			# mll ordered by [pair with leading lepton, pair without leading lepton]
				df['mll'][index] = [mll[i_max[0]], mll[i_max[1]]]
			else:
				df['mll'][index] = [mll[i_max[1]], mll[i_max[0]]]
		
		p = []
		for i in range(len(row['PT'])):
			p.append(TLorentzVector())
			p[i].SetPtEtaPhiM(row['PT'][i], row['Eta'][i], row['Phi'][i], 0.)
		m4l = (p[0] + p[1] + p[2] + p[3]).M()
		df['m4l'][index] = m4l
	print("%d 2e2mu events found" % c)
	df['2e2mu'] = boolean_list
	return df

def jjmuvm(df_j, df_l):
	varset = df_j.columns
	print("Computing 'vbs_tag', 'mjj', 'DeltaEtajj', 'DeltaRjj'")
	for index, row in df_j.iterrows():
		pid = row['PID']
		mjj = []
		DeltaEta_jj = []
		DeltaPhi_jj = []
		DeltaR_jj = []
		pairs = []
		p = []
		for i in range(len(pid)):
			p.append(TLorentzVector())
			p[i].SetPtEtaPhiM(row['PT'][i], row['Eta'][i], row['Phi'][i], 0.)
		for i in range(len(pid)-1):
			for j in range(i+1, len(pid)):
				pairs.append((i,j))
				mjj.append((p[i]+p[j]).M())
				deta = np.abs(p[i].Eta()-p[j].Eta())
				dphi = np.minimum(np.abs(p[i].Phi()-p[j].Phi()), 2*np.pi - np.abs(p[i].Phi()-p[j].Phi()))
				DeltaEta_jj.append(deta)
				DeltaPhi_jj.append(dphi)
				DeltaR_jj.append(np.sqrt(deta**2 + dphi**2))

		w_pole_distance = metric(mjj, Mw)
		i_w = argmaxs(w_pole_distance)[0]
		pair_w = pairs[i_w]
		pair_vbs = (-1,-1)
		mjj_ordered = []
		DeltaEta_jj_ordered = []
		DeltaPhi_jj_ordered = []
		DeltaR_jj_ordered = []
		for u,v in pairs:
			if not (u in pair_w) and not (v in pair_w):
				pair_vbs = (u,v)
		i_vbs = pairs.index(pair_vbs)
		mjj_ordered = [mjj[i_vbs], mjj[i_w]]
		DeltaEta_jj_ordered = [DeltaEta_jj[i_vbs], DeltaEta_jj[i_w]]
		DeltaPhi_jj_ordered = [DeltaPhi_jj[i_vbs], DeltaPhi_jj[i_w]]
		DeltaR_jj_ordered = [DeltaR_jj[i_vbs], DeltaR_jj[i_w]]
		pop_list = []
		for u,v in pairs:
			if ((u,v) == pair_w) or ((u,v) == pair_vbs):
				pop_list.append(pairs.index((u,v)))

		pop_list.reverse()
		for i in pop_list:
			mjj.pop(i)
			DeltaEta_jj.pop(i)
			DeltaPhi_jj.pop(i)
			DeltaR_jj.pop(i)

		i_max = argmaxs(mjj)	# 'other' jets ordered by mjj
		for i in i_max:
			mjj_ordered.append(mjj[i])
			DeltaEta_jj_ordered.append(DeltaEta_jj[i])
			DeltaPhi_jj_ordered.append(DeltaPhi_jj[i])
			DeltaR_jj_ordered.append(DeltaR_jj[i])

		df_j['vbs_tag'][index] = [ i in pair_vbs for i in range(4) ]
		df_j['mjj'][index] = mjj_ordered					# Order: VBS jet pair, W jet pair, others (ordered by mjj)
		df_j['DeltaEta_jj'][index] = DeltaEta_jj_ordered	# Order: VBS jet pair, W jet pair, others (ordered by mjj)
		df_j['DeltaPhi_jj'][index] = DeltaPhi_jj_ordered	# Order: VBS jet pair, W jet pair, others (ordered by mjj)
		df_j['DeltaR_jj'][index] = DeltaR_jj_ordered		# Order: VBS jet pair, W jet pair, others (ordered by mjj)

	print("Computing 'mll', 'PT_miss'")
	varset = df_l.columns
	for index, row in df_l.iterrows():
		pid = row['PID']
		PT_miss = 0.0
		p = []
		first = False
		for i in range(len(pid)):
			p.append(TLorentzVector())
			p[i].SetPtEtaPhiM(row['PT'][i], row['Eta'][i], row['Phi'][i], 0.)
			if np.abs(pid[i]) in [12, 14, 16]:		# if particle is a neutrino, save PT_miss
				PT_miss = p[i].Pt()
				if i == 0:
					first = True
		if first == True:
			for var in varset:
				if type(df_l[var][index]) is type(None):
					continue
				df_l[var][index].reverse()		# Revert lists in order to have [charged lepton, neutrino]
		df_l['mll'][index] = (p[0]+p[1]).M()
		df_l['PT_miss'][index] = PT_miss

	df_j.rename(columns={varset[i] : varset[i] + "_j" for i in range(len(varset))}, inplace=True)
	df_l.rename(columns={varset[i] : varset[i] + "_l" for i in range(len(varset))}, inplace=True)
	df_l.rename(columns={'mll_l' : 'mll', 'PT_miss_l' : 'PT_miss'}, inplace=True)
	df = df_j.copy()
	df[df_l.columns] = df_l
	df['mww'] = [None]*df.shape[0]

	print("Computing 'mww'")
	for index, row in df.iterrows():
		vbs_tag = row['vbs_tag']
		pid = row['PID_l']
		p_j = []
		p_l = []
		i_j = 0
		for i in range(len(vbs_tag)):
			if vbs_tag[i] == False:
				p_j.append(TLorentzVector())
				p_j[i_j].SetPtEtaPhiM(row['PT_j'][i], row['Eta_j'][i], row['Phi_j'][i], 0.)
				i_j += 1
		for i in range(len(pid)):
			p_l.append(TLorentzVector())
			p_l[i].SetPtEtaPhiM(row['PT_l'][i], row['Eta_l'][i], row['Phi_l'][i], 0.)
		df['mww'][index] = (p_j[0]+p_j[1]+p_l[0]+p_l[1]).M()

	return df

def jj(df):
	print("Computing 'mjj', 'DeltaEtajj'")
	for index, row in df.iterrows():
		pid = row['PID']
		p = []
		for i in range(len(pid)):
			p.append(TLorentzVector())
			p[i].SetPtEtaPhiM(row['PT'][i], row['Eta'][i], row['Phi'][i], 0.)

		df['mjj'][index] = (p[0]+p[1]).M()
		df['DeltaEta_jj'][index] = np.abs(p[0].Eta()-p[1].Eta())

	return df

def zz(df_j, df_l, df_z):
	varset = list(df_j.columns)
	pop_list = []
	for (i,var) in enumerate(varset):
		if 'jj' in var:
			pop_list.append(i)
	pop_list.reverse()
	for i in pop_list:
		varset.pop(i)

	print("Computing 'mzz'")
	for index, row in df_z.iterrows():
		pid = row['PID']
		p = []
		for i in range(len(pid)):
			p.append(TLorentzVector())
			p[i].SetPtEtaPhiM(row['PT'][i], row['Eta'][i], row['Phi'][i], row['Mass'][i])

		df_z['mzz'][index] = (p[0]+p[1]).M()

	df_l.rename(columns={'M2' : 'Mothup_l', 'Mass' : 'mz'}, inplace=True)
	df_z.rename(columns={'fUniqueID' : 'fUniqueID_z', 'Mass' : 'mz'}, inplace=True)
	df_j.rename(columns={varset[i] : varset[i] + "_j" for i in range(len(varset))}, inplace=True)
	df_l.rename(columns={varset[i] : varset[i] + "_l" for i in range(len(varset))}, inplace=True)
	df_z.rename(columns={varset[i] : varset[i] + "_z" for i in range(len(varset))}, inplace=True)
	df = df_j.copy()
	df[df_l.columns] = df_l
	df[df_z.columns] = df_z

	print("Computing 'Theta_e', 'Theta_mu'")
	for index, row in df_l.iterrows():
		pid = row['PID_l']
		i_e = None
		i_mu = None
		p_e = TLorentzVector()
		p_mu = TLorentzVector()
		p_zee = TLorentzVector()
		p_zmumu = TLorentzVector()
		for i in range(len(pid)):
			if pid[i] == 11:
				i_e = i
				p_e.SetPtEtaPhiM(row['PT_l'][i], row['Eta_l'][i], row['Phi_l'][i], 0.)
			if pid[i] == 13:
				i_mu = i
				p_mu.SetPtEtaPhiM(row['PT_l'][i], row['Eta_l'][i], row['Phi_l'][i], 0.)

		pid = df['PID_z'][index]
		for i in range(len(pid)):	# This function works if mll[] is ordered by flavor: [ee pair, mumu pair]
			index_z = df['fUniqueID_z'][index][i] - 1
			if index_z == row['Mothup_l'][i_e]:
				p_zee.SetPtEtaPhiM(df['PT_z'][index][i], df['Eta_z'][index][i], df['Phi_z'][index][i], df['mz'][index][i])
			else:
				if index_z == row['Mothup_l'][i_mu]:
					p_zmumu.SetPtEtaPhiM(df['PT_z'][index][i], df['Eta_z'][index][i], df['Phi_z'][index][i], df['mz'][index][i])
				else:
					print("Z not linked to any lepton pair.")

			"""
			if np.abs(df['mll'][index][0] - df['mz'][index][i]) < 0.01:
				p_zee.SetPtEtaPhiM(df['PT_z'][index][i], df['Eta_z'][index][i], df['Phi_z'][index][i], df['mz'][index][i])
			else:
				if np.abs(df['mll'][index][0] - df['mz'][index][i]) < 0.01:
					p_zmumu.SetPtEtaPhiM(df['PT_z'][index][i], df['Eta_z'][index][i], df['Phi_z'][index][i], df['mz_z'][index][i])
				else:
					sys.exit("Z not linked to any lepton pair. Exit.")
			"""
		boost_e = -p_zee.BoostVector()		# boost vector to boost in the Ze c.m. frame
		boost_mu = -p_zmumu.BoostVector()	# boost vector to boost in the Zmu c.m. frame
		p_e.Boost(boost_e)
		p_mu.Boost(boost_mu)

		df['Theta_e'][index] = p_zee.Angle(p_e.Vect())
		df['Theta_mu'][index] = p_zmumu.Angle(p_mu.Vect())
		df['PT_Ze'][index] = p_zee.Pt()
		df['PT_Zmu'][index] = p_zmumu.Pt()
		df['PT_ZZ'][index] = (p_zee + p_zmumu).Pt()
		df['Eta_Ze'][index] = p_zee.Eta()
		df['Eta_Zmu'][index] = p_zmumu.Eta()

	return df	

def save_tree(filename, varset, save=False, label='', path='dataframes/', entrystart=None, entrystop=None, phantom=False):
	print("Opening %s" % filename)
	file = uproot.open(filename)

	tree = file[b'Delphes;1/Particle']
	Nentries = len(tree)
	Nentries_chunk = -1
	print(str(tree.name) + " contains " + str(Nentries) + " entries")

	tree_delphes = file[b'Delphes;1']
	df_part = pd.DataFrame()
	df_events = pd.DataFrame()
	efficiency = 1.0
	if entrystop != None:
		df_part = tree_delphes.pandas.df([b'Particle.fUniqueID', b'Particle.PID', b'Particle.M2', b'Particle.PT', b'Particle.E', b'Particle.Eta', b'Particle.Phi', b'Particle.Mass'], entrystart=entrystart, entrystop=entrystop)
		df_events = tree_delphes.pandas.df([b'Event.ScalePDF'], entrystart=entrystart, entrystop=entrystop)
	else:
		df_part = tree_delphes.pandas.df([b'Particle.fUniqueID', b'Particle.PID', b'Particle.M2', b'Particle.PT', b'Particle.E', b'Particle.Eta', b'Particle.Phi', b'Particle.Mass'])
		df_events = tree_delphes.pandas.df([b'Event.ScalePDF'])

	Nentries_chunk = df_events.shape[0]
	print("Reading " + str(df_events.shape[0]) + " entries")
	if phantom == True:
		singleZ_indices, efficiency = drop_singleZ(df_part, Nentries_chunk)
		df_events.drop(singleZ_indices, inplace=True)
		print("cut_efficiency = %f" % efficiency)
	df_events.rename(columns={'Event.ScalePDF' : 'ScalePDF'}, inplace=True)
	df_events.reset_index(drop=True, inplace=True)
	df_leading = None
	if label in ["zz_mumuee", "z0z0_mumuee", "zTzT_mumuee", "z0zT_mumuee", "zTz0_mumuee"]:
		df_j = df_particle(df_part, 'j')
		df_l = df_particle(df_part, 'l')
		df_z = df_particle(df_part, 'z')
		df_j_leading = build_dataframe(df_j, varset)
		df_l_leading = build_dataframe(df_l, varset + ['Particle.M2'])
		varset_z = varset
		if 'Particle.Mass' not in varset:
			varset_z = varset + ['Particle.fUniqueID', 'Particle.Mass']
		df_z_leading = build_dataframe(df_z, varset_z)
		df_l_leading['mll'] = [[]]*df_l_leading.shape[0]
		df_l_leading['m4l'] = [None]*df_l_leading.shape[0]
		df_l_leading['2e2mu'] = [None]*df_l_leading.shape[0]
		df_l_leading = eemumu(df_l_leading)
		df_j_leading['mjj'] = [None]*df_j_leading.shape[0]
		df_j_leading['DeltaEta_jj'] = [None]*df_j_leading.shape[0]
		df_j_leading = jj(df_j_leading)
		df_z_leading['mzz'] = [None]*df_z_leading.shape[0]
		df_z_leading['PT_ZZ'] = [None]*df_z_leading.shape[0]
		df_z_leading['Theta_e'] = [None]*df_z_leading.shape[0]
		df_z_leading['Theta_mu'] = [None]*df_z_leading.shape[0]
		df_z_leading['PT_Ze'] = [None]*df_z_leading.shape[0]
		df_z_leading['PT_Zmu'] = [None]*df_z_leading.shape[0]
		df_z_leading['Eta_Ze'] = [None]*df_z_leading.shape[0]
		df_z_leading['Eta_Zmu'] = [None]*df_z_leading.shape[0]
		df_leading = zz(df_j_leading, df_l_leading, df_z_leading)
		df_leading['ScalePDF'] = df_events['ScalePDF']
		if phantom == True:
			df_leading['cut_efficiency'] = [efficiency]*df_leading.shape[0]
			df_leading['Nentries_chunk'] = [Nentries_chunk]*df_leading.shape[0]

	if "wp" in label:
		df_j = df_particle(df_part, 'j')
		df_l = df_particle(df_part, 'l')
		df_j_leading = build_dataframe(df_j, varset)
		df_l_leading = build_dataframe(df_l, varset)
		#varset = df_j_leading.columns
		df_j_leading['mjj'] = [[]]*df_j_leading.shape[0]
		df_l_leading['mll'] = [None]*df_l_leading.shape[0]
		df_j_leading['DeltaEta_jj'] = [[]]*df_j_leading.shape[0]
		df_j_leading['DeltaPhi_jj'] = [[]]*df_j_leading.shape[0]
		df_j_leading['DeltaR_jj'] = [[]]*df_j_leading.shape[0]
		df_l_leading['PT_miss'] = [None]*df_l_leading.shape[0]
		df_j_leading['vbs_tag'] = [[]]*df_j_leading.shape[0]
		df_leading = jjmuvm(df_j_leading, df_l_leading)
		df_leading['ScalePDF'] = df_events['ScalePDF']

	if not label in ["zz_mumuee", "z0z0_mumuee", "zTzT_mumuee", "z0zT_mumuee", "zTz0_mumuee", "wp0wm_muvmjj", "wpwm0_jjmuvm"]:
		print("Label name not available.")
		return

	print(df_leading[:10])
	if save == True:
		if entrystop != None:
			path = path + "chunks/"
		if not os.path.exists(path):
			os.makedirs(path)
		out_file = path + 'df_' + label + '_' + str(entrystop) + '.csv'
		if phantom == True:
			out_file = out_file.replace('df_', 'df_phantom_')
		print("Saving dataframe into " + out_file)
		df_leading.to_csv(out_file)
	return df_leading

def selection(df_list, df, var, cfr, val, abs=False):
	values = None
	indices = []
	if type(var) == str:
		values = df[var]
	elif type(var) in [list, np.ndarray, tuple]:
		values = np.array(var)
	if abs == True:
		values = np.abs(values)
	if cfr == '>':
		if type(var) == str:
			indices = [i for (i,j) in df.index[values > val].tolist()]
		else:
			indices = np.where(values > val)[0].tolist()
	elif cfr == '<':
		if type(var) == str:
			indices = [i for (i,j) in df.index[values < val].tolist()]
		else:
			indices = np.where(values < val)[0].tolist()
	elif cfr == '==':
		if type(var) == str:
			indices = [i for (i,j) in df.index[values == val].tolist()]
		else:
			indices = np.where(values == val)[0].tolist()
	elif cfr == '>=':
		if type(var) == str:
			indices = [i for (i,j) in df.index[values >= val].tolist()]
		else:
			indices = np.where(values >= val)[0].tolist()
	elif cfr == '<=':
		if type(var) == str:
			indices = [i for (i,j) in df.index[values <= val].tolist()]
		else:
			indices = np.where(values <= val)[0].tolist()
	df_list_output = []
	for (i, df_part) in enumerate(df_list):
		df_list_output.append(df_part.loc[indices])
	return df_list_output

def n_particles(df_part):
	n_part = []
	key = ''
	if 'Electron' in df_part.columns[0]:
		key = "electrons"
	elif 'Muon' in df_part.columns[0]:
		key = "muons"
	elif 'Jet' in df_part.columns[0]:
		key = "jets"
	else:
		key = "ZZ candidates"
	print("Counting " + key + " per event")
	for i in range(1 + max(df_part.index)[0]):
		if i in df_part.index:
			n_part.append(df_part.loc[i].shape[0])
		else:
			n_part.append(-1)	# This is a 'nan' value
	return np.array(n_part)

def total_charge(df_part):
	var = None
	for col in df_part.columns:
		if 'Charge' in col:
			var = col
	print("Computing total " + var.split('.')[0] + " charge per event")
	tot_charge = []
	for i in range(1 + max(df_part.index)[0]):
		if i in df_part.index:
			tot_charge.append(df_part.loc[i][var].sum())
		else:
			tot_charge.append(100)	# This is a 'nan' value
	return np.abs(tot_charge)

def candidate_pairs(charges):
	index_positive = []
	index_negative = []
	for i in charges.index.tolist():
		if charges[i] > 0:
			index_positive.append(i)
		elif charges[i] < 0:
			index_negative.append(i)
	pairs = [(u,v) for u in index_positive for v in index_negative]	# indices: (e+, e-), (mu+, mu-)
	return pairs

def save_delphes(filename, save=False, label='', path='dataframes/', entrystart=None, entrystop=None):
	print("Opening %s" % filename)
	file = uproot.open(filename)
	filename = filename.split("/")[-1] + ": "	# Filename formatting for standard output
	tree_delphes = file[b'Delphes;1']
	Nevts_gen = len(tree_delphes)
	print(str(tree_delphes.name) + " contains " + str(Nevts_gen) + " events")
	df_Electron = tree_delphes.pandas.df([b'Electron.PT', b'Electron.Eta', b'Electron.Phi', b'Electron.Charge'])
	df_MuonLoose = tree_delphes.pandas.df([b'MuonLoose.PT',b'MuonLoose.Eta',b'MuonLoose.Phi', b'MuonLoose.Charge'])
	df_Jet = tree_delphes.pandas.df([b'Jet.PT', b'Jet.Eta', b'Jet.Phi', b'Jet.Mass'])
	df_Event = tree_delphes.pandas.df([b'Event.ScalePDF'])
	df_Event.index = df_Event.index.droplevel("subentry")

	# Select electrons (muons) with a pT > 7 (5) GeV and jets with a pT > 30 GeV
	print("Selecting electrons (muons) with a pT > 7 (5) GeV")
	df_Electron = df_Electron[df_Electron['Electron.PT'] > 7]
	df_MuonLoose = df_MuonLoose[df_MuonLoose['MuonLoose.PT'] > 5]
	df_Electron = df_Electron[np.abs(df_Electron['Electron.Eta']) < 2.5]
	df_MuonLoose = df_MuonLoose[np.abs(df_MuonLoose['MuonLoose.Eta']) < 2.5]
	print("Selecting jets with a pT > 30 GeV")
	df_Jet = df_Jet[df_Jet['Jet.PT'] > 30]
	df_Jet = df_Jet[np.abs(df_Jet['Jet.Eta']) < 4.7]

	# Select events with at least 2 electrons, 2 muons and 2 jets
	df_list = [df_Electron, df_MuonLoose, df_Jet]
	print(filename, end="")
	n_Electron = n_particles(df_list[0])
	df_list = selection(df_list, df_list[0], n_Electron, '>=', 2)
	print(filename, end="")
	n_MuonLoose = n_particles(df_list[1])
	df_list = selection(df_list, df_list[1], n_MuonLoose, '>=', 2)
	print(filename, end="")
	n_Jet = n_particles(df_list[2])
	df_list = selection(df_list, df_list[2], n_Jet, '>=', 2)

	# Select events with at least 2 Z candidates
	print(filename, end="")
	n_Electron = n_particles(df_list[0])
	print(filename, end="")
	tot_charge_Electron = total_charge(df_list[0])
	df_list = selection(df_list, df_list[0], n_Electron + tot_charge_Electron, '<', 2*n_Electron)
	print(filename, end="")
	n_MuonLoose = n_particles(df_list[1])
	print(filename, end="")
	tot_charge_MuonLoose = total_charge(df_list[1])
	df_list = selection(df_list, df_list[1], n_MuonLoose + tot_charge_MuonLoose, '<', 2*n_MuonLoose)

	# Select the two leading jets
	df_Jet = df_list[2]
	selected_entries = df_Jet.index.get_level_values('entry').tolist()
	selected_entries = list( dict.fromkeys(selected_entries) )		# Remove duplicates from entries list
	indices_to_drop = []
	for i in selected_entries:
		row = df_Jet.loc[i].sort_values(by='Jet.PT', ascending=False)
		index = row.index.tolist()
		for j in index[2:]:
			indices_to_drop.append((i,j))
	df_Jet = df_Jet.drop(indices_to_drop)
	df_list[2] = df_Jet

	# Construct dataframe of ZZ candidates
	print(filename, end="")
	print("Constructing ZZ candidates")
	selected_entries = df_list[0].index.get_level_values('entry').tolist()
	selected_entries = list( dict.fromkeys(selected_entries) )		# Remove duplicates from entries list
	max_Z_candidates = max(n_Electron)*max(n_MuonLoose)
	iterables = [selected_entries, range(max_Z_candidates)]
	index = pd.MultiIndex.from_product(iterables, names=['entry', 'subentry'])
	columns = ['Z2.ScalarPTsum', 'Lepton.PT1', 'Lepton.PT2', 'Lepton.PT3', 'Lepton.PT4',
			   'Z1.Mass', 'Z2.Mass', 'e+mu-.Mass', 'e-mu+.Mass', 'Z1.PoleDistance', 'mzz',
			   'mjj', 'DeltaEta_jj',
			   'PT_l', 'Eta_l', 'mll',
			   'PT_j', 'Eta_j', 'Phi_j',
			   'PT_z', 'Eta_z', 'Phi_z', 'mz', 'PT_Ze', 'PT_Zmu', 'Eta_Ze','Eta_Zmu', 'PT_ZZ', 'm4l',
			   'Theta_e', 'Theta_mu', 'ScalePDF']
	DeltaR_columns = ['Jet.DeltaR' + str(i) for i in range(1,9)]
	columns = columns + DeltaR_columns
	df_Candidates = pd.DataFrame(index=index, columns=columns)

	for i in selected_entries:
		row_jet = df_list[2].loc[i]
		jet_index = row_jet.index
		p_Jet = []
		for j in jet_index:
			p = TLorentzVector()
			p.SetPtEtaPhiM(row_jet.loc[j]['Jet.PT'], row_jet.loc[j]['Jet.Eta'], row_jet.loc[j]['Jet.Phi'], row_jet.loc[j]['Jet.Mass'])
			p_Jet.append(p)
		mjj = (p_Jet[0] + p_Jet[1]).M()
		DeltaEta_jj = np.abs(p_Jet[0].Eta() - p_Jet[1].Eta())
		Jet_PT = [p_Jet[0].Pt(), p_Jet[1].Pt()]
		Jet_PT.sort(reverse=True)
		Electron_pairs = candidate_pairs(df_list[0].loc[i]['Electron.Charge'])
		MuonLoose_pairs = candidate_pairs(df_list[1].loc[i]['MuonLoose.Charge'])
		n_Candidates = 0
		n_LeptonPairs = 0
		for pair_e in Electron_pairs:
			Electron_PT = []
			p_Electron = []
			for j in pair_e:
				Electron_PT.append(df_list[0].loc[(i,j)]['Electron.PT'])
				p = TLorentzVector()
				p.SetPtEtaPhiM(df_list[0].loc[(i,j)]['Electron.PT'], df_list[0].loc[(i,j)]['Electron.Eta'], df_list[0].loc[(i,j)]['Electron.Phi'], 0.)
				p_Electron.append(p)
			p_zee = p_Electron[0] + p_Electron[1]
			mz_e = p_zee.M()
			PTsum_e = p_Electron[0].Pt() + p_Electron[1].Pt()
			PT_Ze = p_zee.Pt()
			Eta_Ze = p_zee.Eta()
			Phi_Ze = p_zee.Phi()
			boost_e = -p_zee.BoostVector()		# boost vector to boost in the Ze c.m. frame
			p_e = TLorentzVector(p_Electron[1])
			p_e.Boost(boost_e)
			Theta_e = p_zee.Angle(p_e.Vect())		# p_Electron[1] corresponds to e-
			for pair_mu in MuonLoose_pairs:
				MuonLoose_PT = []
				p_MuonLoose = []
				for k in pair_mu:
					MuonLoose_PT.append(df_list[1].loc[(i,k)]['MuonLoose.PT'])
					p = TLorentzVector()
					p.SetPtEtaPhiM(df_list[1].loc[(i,k)]['MuonLoose.PT'], df_list[1].loc[(i,k)]['MuonLoose.Eta'], df_list[1].loc[(i,k)]['MuonLoose.Phi'], 0.)
					p_MuonLoose.append(p)
				p_zmumu = p_MuonLoose[0] + p_MuonLoose[1]
				mz_mu = p_zmumu.M()
				PTsum_mu = p_MuonLoose[0].Pt() + p_MuonLoose[1].Pt()
				PT_Zmu = p_zmumu.Pt()
				Eta_Zmu = p_zmumu.Eta()
				Phi_Zmu = p_zmumu.Phi()
				boost_mu = -p_zmumu.BoostVector()	# boost vector to boost in the Zmu c.m. frame
				p_mu = TLorentzVector(p_MuonLoose[1])
				p_mu.Boost(boost_mu)
				Theta_mu = p_zmumu.Angle(p_mu.Vect())		# p_MuonLoose[1] corresponds to mu-
				mz = [mz_e, mz_mu]
				PTsum = [PTsum_e, PTsum_mu]
				PT_z = [PT_Ze, PT_Zmu]
				Eta_z = [Eta_Ze, Eta_Zmu]
				Phi_z = [Phi_Ze, Phi_Zmu]
				z_pole_distance = metric(mz, Mz)
				if z_pole_distance[0] < z_pole_distance[1]:	# N.B. metric returns a negative value for practical reasons
					mz.reverse()
					PTsum.reverse()
					z_pole_distance.reverse()
				Lepton_PT = Electron_PT + MuonLoose_PT
				Lepton_PT.sort(reverse=True)
				mll = []
				mll.append((p_Electron[0] + p_MuonLoose[1]).M())
				mll.append((p_Electron[1] + p_MuonLoose[0]).M())
				mzz = (p_Electron[0] + p_Electron[1] + p_MuonLoose[0] + p_MuonLoose[1]).M()
				PT_ZZ = (p_Electron[0] + p_Electron[1] + p_MuonLoose[0] + p_MuonLoose[1]).Pt()
				p_Lepton = p_Electron + p_MuonLoose
				DeltaR = []
				for p in p_Jet:
					for q in p_Lepton:
						DeltaR.append(p.DeltaR(q))
				for (l, col) in enumerate(DeltaR_columns):
					df_Candidates[col][(i, n_Candidates)] = DeltaR[l]
				df_Candidates['Z2.ScalarPTsum'][(i,n_Candidates)] = PTsum[1]
				df_Candidates['Lepton.PT1'][(i,n_Candidates)] = Lepton_PT[0]
				df_Candidates['Lepton.PT2'][(i,n_Candidates)] = Lepton_PT[1]
				df_Candidates['Lepton.PT3'][(i,n_Candidates)] = Lepton_PT[2]
				df_Candidates['Lepton.PT4'][(i,n_Candidates)] = Lepton_PT[3]
				df_Candidates['Z1.Mass'][(i,n_Candidates)] = mz[0]
				df_Candidates['Z2.Mass'][(i,n_Candidates)] = mz[1]
				df_Candidates['e+mu-.Mass'][(i,n_Candidates)] = mll[0]
				df_Candidates['e-mu+.Mass'][(i,n_Candidates)] = mll[1]
				df_Candidates['Z1.PoleDistance'][(i,n_Candidates)] = np.abs(z_pole_distance[0])
				df_Candidates['mzz'][(i, n_Candidates)] = mzz
				df_Candidates['m4l'][(i, n_Candidates)] = mzz
				df_Candidates['PT_j'][(i, n_Candidates)] = Jet_PT
				df_Candidates['Eta_j'][(i, n_Candidates)] = [p_Jet[0].Eta(), p_Jet[1].Eta()]
				df_Candidates['Phi_j'][(i, n_Candidates)] = [p_Jet[0].Phi(), p_Jet[1].Phi()]
				df_Candidates['mjj'][(i, n_Candidates)] = mjj
				df_Candidates['DeltaEta_jj'][(i, n_Candidates)] = DeltaEta_jj
				df_Candidates['PT_l'][(i, n_Candidates)] = Lepton_PT
				df_Candidates['Eta_l'][(i, n_Candidates)] = [p_Electron[0].Eta(), p_Electron[1].Eta(), p_MuonLoose[0].Eta(), p_MuonLoose[1].Eta()]
				df_Candidates['mll'][(i, n_Candidates)] = mz
				df_Candidates['mz'][(i, n_Candidates)] = mz
				df_Candidates['PT_z'][(i, n_Candidates)] = PT_z
				df_Candidates['PT_Ze'][(i, n_Candidates)] = PT_Ze
				df_Candidates['PT_Zmu'][(i, n_Candidates)] = PT_Zmu
				df_Candidates['Eta_z'][(i, n_Candidates)] = Eta_z
				df_Candidates['Eta_Ze'][(i, n_Candidates)] = Eta_Ze
				df_Candidates['Eta_Zmu'][(i, n_Candidates)] = Eta_Zmu
				df_Candidates['Phi_z'][(i, n_Candidates)] = Phi_z
				df_Candidates['PT_ZZ'][(i, n_Candidates)] = PT_ZZ
				df_Candidates['Theta_e'][(i, n_Candidates)] = Theta_e
				df_Candidates['Theta_mu'][(i, n_Candidates)] = Theta_mu
				df_Candidates['ScalePDF'][(i, n_Candidates)] = df_Event['Event.ScalePDF'].loc[i]
				n_Candidates += 1
	df_Candidates.dropna(inplace=True)
	# Selections on PT1, PT2, PT3, PT4, mZ, mll
	print("candidates = %d" % df_Candidates.shape[0])
	df_Candidates = df_Candidates[df_Candidates['Lepton.PT1'] > 20]
	print("PT1 selection...")
	print("candidates = %d" % df_Candidates.shape[0])
	df_Candidates = df_Candidates[df_Candidates['Lepton.PT2'] > 10]
	print("PT2 selection...")
	print("candidates = %d" % df_Candidates.shape[0])
	df_Candidates = df_Candidates[df_Candidates['Lepton.PT3'] > 10]
	print("PT3 selection...")
	print("candidates = %d" % df_Candidates.shape[0])
	df_Candidates = df_Candidates[df_Candidates['Lepton.PT4'] > 10]
	print("PT4 selection...")
	print("candidates = %d" % df_Candidates.shape[0])
	df_Candidates = df_Candidates[df_Candidates['Z1.Mass'] > 60]
	df_Candidates = df_Candidates[df_Candidates['Z1.Mass'] < 120]
	df_Candidates = df_Candidates[df_Candidates['Z2.Mass'] > 60]
	df_Candidates = df_Candidates[df_Candidates['Z2.Mass'] < 120]
	print("mz selection...")
	print("candidates = %d" % df_Candidates.shape[0])
	df_Candidates = df_Candidates[df_Candidates['e+mu-.Mass'] > 4]
	df_Candidates = df_Candidates[df_Candidates['e-mu+.Mass'] > 4]
	print("mll selection...")
	print("candidates = %d" % df_Candidates.shape[0])
	df_Candidates = df_Candidates[df_Candidates['mzz'] > 180]
	print("mzz selection...")
	print("candidates = %d" % df_Candidates.shape[0])
	for (l, col) in enumerate(DeltaR_columns):
		df_Candidates = df_Candidates[df_Candidates[col] > 0.4]
	print("DeltaR selection...")
	print("candidates = %d" % df_Candidates.shape[0])
	df_Candidates = df_Candidates[df_Candidates['mjj'] > 100]
	print("mjj selection...")
	print("candidates = %d" % df_Candidates.shape[0])

	# Choose the candidate with the largest scalar sum of the pT of the Z2 leptons
	selected_entries = df_Candidates.index.get_level_values('entry').tolist()
	selected_entries = list( dict.fromkeys(selected_entries) )		# Remove duplicates from entries list
	indices_to_drop = []
	for i in selected_entries:
		if len(df_Candidates.loc[i]) > 1:
			row = df_Candidates.loc[i].sort_values(by='Z2.ScalarPTsum', ascending=False)
			index = row.index.tolist()
			for j in index[1:]:
				indices_to_drop.append((i,j))
	df_Candidates = df_Candidates.drop(indices_to_drop)
	df_Candidates.index = df_Candidates.index.droplevel("subentry")	# Drop meaningless subentry index

	# Select entries of Electron, MuonLoose and Jet dataframes
	selected_entries = df_Candidates.index.get_level_values('entry').tolist()
	selected_entries = list( dict.fromkeys(selected_entries) )		# Remove duplicates from entries list
	df_list_output = []
	for (i, df_part) in enumerate(df_list):
		df_list_output.append(df_part.loc[selected_entries])

	print(df_Candidates[:10])
	if save == True:
		if entrystop != None:
			path = path + "chunks/"
			out_file = path + 'df_delphes_' + label + '_' + str(entrystop) + '.csv'
		else:
			out_file = path + 'df_delphes_' + label + '.csv'
		if not os.path.exists(path):
			os.makedirs(path)
		print("Saving dataframe into " + out_file)
		df_Candidates.to_csv(out_file)

	#return df_Candidates, df_list_output
	return df_Candidates

# Functions for reading dataframe and cuts selections
def create_lists(df, varset, entrystart=None, entrystop=None):
	chunk = None
	if entrystop != None:
		print("Creating lists for dataframe chunk [%d, %d]" % (entrystart, entrystop))
		chunk = df.loc[entrystart:entrystop-1].copy()
	else:
		print("Creating lists for dataframe")
		chunk = df.copy()
	for index, row in chunk.iterrows():
		for var in varset:
			chunk[var][index] = list(map(float, row[var][1:-1].split(',')))
	return chunk

def threshold(column, val, cfr='>', abs=False):
	length = len(column.values[0])
	val_list = [length*[val]]
	val_list = column.shape[0]*val_list
	val_list = np.array(val_list)
	col_list = column.values.tolist()
	boolean_list = []
	if abs == True:
		x = np.abs(np.array(col_list))
	else:
		x = np.array(col_list)
		
	for (i, item) in enumerate(col_list):
		if cfr == '>':
			boolean_list.append((x[i] > np.array(val_list)[i]).all())
		else:
			if cfr == '<':
				boolean_list.append((x[i] < np.array(val_list)[i]).all())
	return boolean_list

def column_leading(column, hierarchy):
	leading_list = np.array([])
	if type(column) in [list, np.ndarray]:
		if hierarchy == 1:
			leading_list = np.array(column)[:,0]
		if hierarchy == 2:
			leading_list = np.array(column)[:,1]
		if hierarchy == 3:
			leading_list = np.array(column)[:,2]
		if hierarchy == 4:
			leading_list = np.array(column)[:,3]
	else:
		if hierarchy == 1:
			leading_list = np.array(column.values.tolist())[:,0]
		if hierarchy == 2:
			leading_list = np.array(column.values.tolist())[:,1]
		if hierarchy == 3:
			leading_list = np.array(column.values.tolist())[:,2]
		if hierarchy == 4:
			leading_list = np.array(column.values.tolist())[:,3]

	return leading_list

def threshold_leading(column, hierarchy, val, cfr='>', abs=False):
	leading_list = column_leading(column, hierarchy)
			
	val_list = val*np.ones_like(leading_list)
	
	boolean_list = []
	if abs == True:
		x = np.abs(np.array(leading_list))
	else:
		x = np.array(leading_list)
		
	for (i, item) in enumerate(leading_list):
		if cfr == '>':
			boolean_list = (x > np.array(val_list))
		else:
			if cfr == '<':
				boolean_list = (x < np.array(val_list))
	return boolean_list

def inclusive(column):		# Returns a flat list out of a column with list of lists
	return np.array(column.values.tolist()).flatten().tolist()

def sort_column(column, reverse=False, key=None):	# Returns the column with sorted lists
	sorted_column = []
	for item in column:
		row = item.copy()
		row.sort(reverse=reverse, key=key)
		sorted_column.append(row)
	return sorted_column

