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

def build_col(df, varset):   # Leading and sub-leading particles are defined with respect to
	varname = varset[0]      # the first variable of the list (varset[0]) (normally PT)
	x = []					 # The function generates lists of the type PT=[40.5, 36, 20, 10]
	list_var = []
	for (j, var) in enumerate(varset):
		list_var.append([])
	length = len(df.loc[0][varname].values)
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
	print("Computing 'vbs_tag', 'mjj', 'DeltaEtajj'")
	for index, row in df_j.iterrows():
		pid = row['PID']
		mjj = []
		DeltaEta_jj = []
		pairs = []
		p = []
		for i in range(len(pid)):
			p.append(TLorentzVector())
			p[i].SetPtEtaPhiM(row['PT'][i], row['Eta'][i], row['Phi'][i], 0.)
		for i in range(len(pid)-1):
			for j in range(i+1, len(pid)):
				pairs.append((i,j))
				mjj.append((p[i]+p[j]).M())
				DeltaEta_jj.append(np.abs(p[i].Eta()-p[j].Eta()))

		w_pole_distance = metric(mjj, Mw)
		i_w = argmaxs(w_pole_distance)[0]
		pair_w = pairs[i_w]
		pair_vbs = (-1,-1)
		mjj_ordered = []
		DeltaEta_jj_ordered = []
		for u,v in pairs:
			if not (u in pair_w) and not (v in pair_w):
				pair_vbs = (u,v)
		i_vbs = pairs.index(pair_vbs)
		mjj_ordered = [mjj[i_vbs], mjj[i_w]]
		DeltaEta_jj_ordered = [DeltaEta_jj[i_vbs], DeltaEta_jj[i_w]]
		pop_list = []
		for u,v in pairs:
			if ((u,v) == pair_w) or ((u,v) == pair_vbs):
				pop_list.append(pairs.index((u,v)))

		pop_list.reverse()
		for i in pop_list:
			mjj.pop(i)
			DeltaEta_jj.pop(i)

		i_max = argmaxs(mjj)	# 'other' jets ordered by mjj
		for i in i_max:
			mjj_ordered.append(mjj[i])
			DeltaEta_jj_ordered.append(DeltaEta_jj[i])

		df_j['vbs_tag'][index] = [ i in pair_vbs for i in range(4) ]
		df_j['mjj'][index] = mjj_ordered					# Order: VBS jet pair, W jet pair, others (ordered by mjj)
		df_j['DeltaEta_jj'][index] = DeltaEta_jj_ordered	# Order: VBS jet pair, W jet pair, others (ordered by mjj)

	print("Computing 'mll', 'PT_miss'")
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
				df_l[var][index].reverse()		# Revert lists in order to have [charged lepton, neutrino]
		df_l['mll'][index] = (p[0]+p[1]).M()
		df_l['PT_miss'][index] = PT_miss

	df_j.rename(columns={varset[i] : varset[i] + "_j" for i in range(len(varset))}, inplace=True)
	df_l.rename(columns={varset[i] : varset[i] + "_l" for i in range(len(varset))}, inplace=True)
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
	df['Theta_e'] = [None]*df.shape[0]
	df['Theta_mu'] = [None]*df.shape[0]

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

	return df	


def save_tree(filename, varset, part_type, save=False, label='', path='dataframes/'):
	print("Opening %s" % filename)
	file = uproot.open(filename)

	#events = file[b'Delphes;1/Event']
	tree = file[b'Delphes;1/Particle']
	#print(str(events.name) + " contains " + str(len(events)) + " entries")
	print(str(tree.name) + " contains " + str(len(tree)) + " entries")

	tree_delphes = file[b'Delphes;1']
	#df_part = tree_delphes.pandas.df([b'Particle.PID', b'Particle.PT', b'Particle.E', b'Particle.Eta', b'Particle.Phi', b'Particle.Mass', b'Event.ScalePDF'])
	#df_part = tree_delphes.pandas.df([b'Particle.PID', b'Particle.PT', b'Particle.E', b'Particle.Eta', b'Particle.Phi', b'Particle.Mass'])
	df_part = tree_delphes.pandas.df([b'Particle.fUniqueID', b'Particle.PID', b'Particle.M2', b'Particle.PT', b'Particle.E', b'Particle.Eta', b'Particle.Phi', b'Particle.Mass'])
	df_leading = None
	if label in ["zz_mumuee", "zz0_mumuee", "zzT_mumuee"]:
		if part_type == 'l':
			df_p = df_particle(df_part, part_type)
			df_leading = build_dataframe(df_p, varset)
			df_leading['mll'] = [[]]*df_leading.shape[0]
			df_leading['m4l'] = [None]*df_leading.shape[0]
			df_leading['2e2mu'] = [None]*df_leading.shape[0]
			df_leading = eemumu(df_leading)
		if part_type == 'j':
			df_p = df_particle(df_part, part_type)
			df_leading = build_dataframe(df_j, varset)
			df_leading['mjj'] = [None]*df_leading.shape[0]
			df_leading['DeltaEta_jj'] = [None]*df_leading.shape[0]
			df_leading = jj(df_j_leading)
		if part_type == "final":
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
			df_leading = zz(df_j_leading, df_l_leading, df_z_leading)

	if label == "wpwm0_jjmuvm":
		if part_type == "final":
			df_j = df_particle(df_part, 'j')
			df_l = df_particle(df_part, 'l')
			df_j_leading = build_dataframe(df_j, varset)
			df_l_leading = build_dataframe(df_l, varset)
			#varset = df_j_leading.columns
			df_j_leading['mjj'] = [[]]*df_j_leading.shape[0]
			df_l_leading['mll'] = [None]*df_l_leading.shape[0]
			df_j_leading['DeltaEta_jj'] = [[]]*df_j_leading.shape[0]
			df_l_leading['PT_miss'] = [None]*df_l_leading.shape[0]
			df_j_leading['vbs_tag'] = [[]]*df_j_leading.shape[0]
			df_leading = jjmuvm(df_j_leading, df_l_leading)

	if not label in ["zz_mumuee", "zz0_mumuee", "zzT_mumuee", "wpwm0_jjmuvm"]:
		print("Label name not available.")
		return

	print(df_leading[:10])
	if save == True:
		if not os.path.exists(path):
			os.makedirs(path)
		out_file = path + 'df_' + part_type + '_' + label + '.csv'
		print("Saving dataframe into " + out_file)
		df_leading.to_csv(out_file)
	return df_leading

# Functions for reading dataframe and cuts selections

def create_lists(df, varset):
	for index, row in df.iterrows():
		list_of_lists = []
		for var in varset:
			df[var][index] = list(map(float, row[var][1:-1].split(',')))

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

def inclusive(column):		# Return a flat list out of a column with list of lists
	return np.array(column.values.tolist()).flatten().tolist()
