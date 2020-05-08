import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

d = {"pt" : "$P_{T}$ spectrum", "eta" : "$\eta$ distribution", "mll" : "$M_{ll}$ spectrum", "mjj" : "$M_{jj}$ spectrum", "deltaetajj" : "$\Delta\eta_{jj}$ distribution",
	 "mww" : "$M_{WW}$ spectrum", "m4l" : "$M_{4l}$ spectrum", "mzz" : "$M_{ZZ}$ spectrum", "pdfscale" : "Renormalization scale spectrum",
	 "costheta" : "cosθ distribution", "deltaphijj" : "$\Delta\phi_{jj}$ distribution"}
d_label = {"pt" : "$P_{T}$ [GeV]", "eta" : "$\eta$", "mll" : "$M_{ll}$ [GeV]", "mjj" : "$M_{jj}$ [GeV]", "deltaetajj" : "$\Delta\eta_{jj}$",
		   "mww" : "$M_{WW}$ [GeV]", "m4l" : "$M_{4l}$ [GeV]", "mzz" : "$M_{ZZ}$ [GeV]", "pdfscale" : "$M_{ZZ}/\sqrt{2}$ [GeV]",
		   "costheta" : "cosθ", "deltaphijj" : "$\Delta\phi_{jj}$"}

class plot(object):
  
	def __init__(self, name, var, bin, label='', cut=None, color='blue', type='hist', plot_dir='plots/', filelabel='', save=False):
		self.name = name
		self.var = var
		self.bin = bin
		self.label = label
		self.cut = cut
		self.color = color
		self.type = type
		self.plot_dir = plot_dir
		self.filelabel = filelabel
		self.save = save
		self.varname = ""
		self.part = ""
		self.title = ""
		if len(name.split("_")) > 1:
			self.varname = self.name.split("_")[0]
			self.part = self.name.split("_")[1]
			self.title = (self.part + " " + d[self.varname]).capitalize()
		else:
			self.varname = self.name.split("_")[0]
			self.title = d[self.varname]
		self.xlabel = d_label[self.varname]
		self.ylabel = ""
		if type == "hist":
			self.ylabel = "# entries"
		if type == "scatter" and ("pdfscale_vs_" in self.name):
			if "wp" in self.plot_dir:
				self.xlabel = "$M_{WW}/\sqrt{2}$ [GeV]"
			if "z" in self.plot_dir:
				self.xlabel = "$M_{ZZ}/\sqrt{2}$ [GeV]"
			self.ylabel = "ScalePDF [GeV]"
		if type == "scatter" and self.name == "mzz_vs_m4l":
			if "z" in self.plot_dir:
				self.xlabel = "$M_{4l}$ [GeV]"
				self.ylabel = "$M_{ZZ}$ [GeV]"
		if type == "scatter" and self.name == "m4l_vs_mzz":
			if "z" in self.plot_dir:
				self.xlabel = "$M_{ZZ}$ [GeV]"
				self.ylabel = "$M_{4l}$ [GeV]"
		#print("Initializing " + self.name + " object.")

	def __del__(self):
		plt.close()
		#print("Deleting " + self.name + " object.")

	def draw(self):
		if self.type == "hist":
			plt.figure(figsize=[12, 9])
			n_list = []
			bins_list = []
			label_list = []
			if type(self.var[0]) is list:
				for (i, x) in enumerate(self.var):
					n, bins, patches = plt.hist(x, bins=self.bin, color=self.color[i], alpha=0.7, label=self.label[i], histtype="stepfilled", ec="black")
					n_list.append(n)
					bins_list.append([x + 0.5*(bins[1]-bins[0]) for x in bins[:-1]])
					label_list.append(self.label[i])
			else:
				n, bins, patches = plt.hist(self.var, bins=self.bin, color=self.color, alpha=0.7, label=self.label, histtype="stepfilled", ec="black")
				n_list.append(n)
				bins_list.append([x + 0.5*(bins[1]-bins[0]) for x in bins[:-1]])
				label_list.append(self.label)
			plt.title(self.title)
			if self.cut != None:
				cut_label = d_label[self.varname].split('[')[0].strip()
				if type(self.cut) is list:
					for (i, x) in enumerate(self.cut):
						if type(self.color) is list:
							plt.vlines(x, 0, max(max(n) for n in n_list), linestyles='dashed', label=cut_label + ' ' + self.label[i] + ' cut', color=self.color[i])
							if self.varname == "eta":
								plt.vlines(-x, 0, max(max(n) for n in n_list), linestyles='dashed', color=self.color[i])
						else:
							plt.vlines(x, 0, max(max(n) for n in n_list), linestyles='dashed', label=cut_label + ' ' + self.label[i] + ' cut', color='grey')
							if self.varname == "eta":
								plt.vlines(-x, 0, max(max(n) for n in n_list), linestyles='dashed', color='grey')
				else:
					plt.vlines(self.cut, 0, max(max(n) for n in n_list), linestyles='dashed', label=cut_label + ' ' + self.label + ' cut', color='grey')
					if self.varname == "eta":
						plt.vlines(-self.cut, 0, max(max(n) for n in n_list), linestyles='dashed', color='grey')

			plt.legend(loc="best")
			plt.xlabel(self.xlabel)
			plt.ylabel(self.ylabel)
			plot_file = self.plot_dir + self.name + ".png"
			if self.filelabel != "":
				plot_file = plot_file.replace(".png", "_" + self.filelabel + ".png")
			plt.savefig(plot_file, format="png")
			print("Saving " + plot_file)
			if self.save == True:
				if type(self.label) is not list:
					self.label = [self.label]
				bins_dict = {'bins' + str(i+1) : bins_list[i] for i in range(len(bins_list))}
				n_dict = {'n' + str(i+1) : n_list[i] for i in range(len(n_list))}
				label_dict = {'label' + str(i+1) : len(bins_list[0])*[self.label[i]] for i in range(len(self.label))}
				d = dict({})
				for dictionnary in [bins_dict, n_dict, label_dict]:
					d.update(dictionnary)
				df = pd.DataFrame(data=d)
				dataframes_dir = self.plot_dir.replace("plots", "dataframes")
				dataframes_dir = dataframes_dir + "histograms/"
				if not os.path.exists(dataframes_dir):
					os.makedirs(dataframes_dir)
				out_file = dataframes_dir + 'df_' + self.name + '_hist.csv'
				if self.filelabel != "":
					out_file = out_file.replace('_hist.csv', '_' + self.filelabel + '_hist.csv')
				print("Saving " + out_file)
				df.to_csv(out_file)
			plt.show()
			plt.close()
		else:
			if self.type == "scatter":
				plt.figure(figsize=[6, 6])
				plt.scatter(self.var[0], self.var[1], s=1, color="blue")
				plt.xlim(self.bin[0], self.bin[1])
				plt.ylim(self.bin[0], self.bin[1])
				plt.title(self.ylabel.split("[")[0].strip() + " vs " + self.xlabel.split("[")[0].strip())
				plt.xlabel(self.xlabel)
				plt.ylabel(self.ylabel)
				plot_file = self.plot_dir + self.name + ".png"
				if self.filelabel != "":
					plot_file = plot_file.replace(".png", "_" + self.filelabel + ".png")
				plt.savefig(plot_file, format="png")
				print("Saving " + plot_file + ".png")
				plt.show()
				plt.close()
		return


