import h5py
import numpy as np

with h5py.File('../data/output/combined.hdf5', 'r') as hf:
	mode = np.array(hf['mode'][()])
	ipnu = np.array(hf['ipnu'][()])
	oscw = np.array(hf['weightOsc_SKpaper'][()])
	weight = np.array(hf['weightReco'][()])
	itype = np.array(hf['itype'][()])

wght = np.multiply(oscw,weight)

ccnue = np.zeros(16)
ccnueb = np.zeros(16)
ccmu = np.zeros(16)
cctau = np.zeros(16)
nc = np.zeros(16)
totosc = np.zeros(16)
totnoosc = np.zeros(16)

j = -1
tot = 0.0
for mod, nu, w, typ in zip(mode, ipnu, wght, itype):
	tot = tot + w
	if (typ>=0 and typ<16):
		j = int(typ)
		totosc[j] = totosc[j] + w
		totnoosc[j] = totnoosc[j] + 1
		if (nu==12 and abs(mod)<30):
			ccnue[j] = ccnue[j] + w
		if (nu==-12 and abs(mod)<30):
			ccnueb[j] = ccnueb[j] + w
		if (abs(nu)==14 and abs(mod)<30):
			ccmu[j] = ccmu[j] + w
		if (abs(nu)==16 and abs(mod)<30):
			cctau[j] = cctau[j] + w
		if (abs(mod)>30):
			nc[j] = nc[j] + w



with open("TableOfContents.txt", "w") as table:
	table.write("\\begin{table*}[ht]\n")
# 	print(f'{value:.6f}')
	table.write("\\centering\n")
	table.write("\\begin{tabular}{p{0.2\\linewidth}p{0.09\\linewidth}p{0.09\\linewidth}p{0.12\\linewidth}p{0.09\\linewidth}p{0.09\\linewidth}p{0.09\\linewidth}}\n")
	table.write("Sample & CC $\\nu_e$ & CC $\\overline{\\nu}_e$ & CC $\\nu_\\mu+\\overline{\\nu}_\\mu$ & CC $\\nu_\\tau$ & NC& Fraction \\\\\n")
	table.write("\\hline\n")
	table.write("\\hline\n")
	table.write("\\vspace{0.1cm}\n")
	table.write("\\textbf{Sub-GeV}\\\\\n")
	table.write("\\text{  }Single-Ring, e-like\\\\\n")
	table.write("\\text{  }\\text{  }0 decay-e & \n")
	table.write(str("%.3f" % ( ccnue[0]/totosc[0])) + " & " + str("%.3f" % ( ccnueb[0]/totosc[0])) + " & " + str("%.3f" % ( ccmu[0]/totosc[0])) + " & " + str("%.3f" % ( cctau[0]/totosc[0])) + " & " + str("%.3f" % ( nc[0]/totosc[0])) + " & " + str("%.3f" % ( totosc[0]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }$\\ge$1 decay-e & \n")
	table.write(str("%.3f" % ( ccnue[1]/totosc[1])) + " & " + str("%.3f" % ( ccnueb[1]/totosc[1])) + " & " + str("%.3f" % ( ccmu[1]/totosc[1])) + " & " + str("%.3f" % ( cctau[1]/totosc[1])) + " & " + str("%.3f" % ( nc[1]/totosc[1])) +  " & " + str("%.3f" % ( totosc[1]/tot)) + " \\\\\n")
	table.write("\\text{  }Single-Ring, $\\mu$-like\\\\\n")
	table.write("\\text{  }\\text{  }0 decay-e & \n")
	table.write(str("%.3f" % ( ccnue[3]/totosc[3])) + " & " + str("%.3f" % ( ccnueb[3]/totosc[3])) + " & " + str("%.3f" % ( ccmu[3]/totosc[3])) + " & " + str("%.3f" % ( cctau[3]/totosc[3])) + " & " + str("%.3f" % ( nc[3]/totosc[3])) +  " & " + str("%.3f" % ( totosc[3]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }1 decay-e & \n")
	table.write(str("%.3f" % ( ccnue[4]/totosc[4])) + " & " + str("%.3f" % ( ccnueb[4]/totosc[4])) + " & " + str("%.3f" % ( ccmu[4]/totosc[4])) + " & " + str("%.3f" % ( cctau[4]/totosc[4])) + " & " + str("%.3f" % ( nc[4]/totosc[4])) +  " & " + str("%.3f" % ( totosc[4]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }$\\ge$2 decay-e & \n")
	table.write(str("%.3f" % ( ccnue[5]/totosc[5])) + " & " + str("%.3f" % ( ccnueb[5]/totosc[5])) + " & " + str("%.3f" % ( ccmu[5]/totosc[5])) + " & " + str("%.3f" % ( cctau[5]/totosc[5])) + " & " + str("%.3f" % ( nc[5]/totosc[5])) +  " & " + str("%.3f" % ( totosc[5]/tot)) + " \\\\\n")
	table.write("\\text{  }$\\pi^0$-like\\\\\n")
	table.write("\\text{  }\\text{  }Single-Ring & \n")
	table.write(str("%.3f" % ( ccnue[2]/totosc[2])) + " & " + str("%.3f" % ( ccnueb[2]/totosc[2])) + " & " + str("%.3f" % ( ccmu[2]/totosc[2])) + " & " + str("%.3f" % ( cctau[2]/totosc[2])) + " & " + str("%.3f" % ( nc[2]/totosc[2])) +  " & " + str("%.3f" % ( totosc[2]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }Multi-Ring  & \n")
	table.write(str("%.3f" % ( ccnue[6]/totosc[6])) + " & " + str("%.3f" % ( ccnueb[6]/totosc[6])) + " & " + str("%.3f" % ( ccmu[6]/totosc[6])) + " & " + str("%.3f" % ( cctau[6]/totosc[6])) + " & " + str("%.3f" % ( nc[6]/totosc[6])) +  " & " + str("%.3f" % ( totosc[6]/tot)) + " \\\\\n")
	table.write("\\vspace{0.1cm}\n")
	table.write("\\textbf{Multi-GeV}\\\\\n")
	table.write("\\text{  }Single-Ring\\\\\n")
	table.write("\\text{  }\\text{  }$\\nu_e$-like & \n")
	table.write(str("%.3f" % ( ccnue[7]/totosc[7])) + " & " + str("%.3f" % ( ccnueb[7]/totosc[7])) + " & " + str("%.3f" % ( ccmu[7]/totosc[7])) + " & " + str("%.3f" % ( cctau[7]/totosc[7])) + " & " + str("%.3f" % ( nc[7]/totosc[7])) +  " & " + str("%.3f" % ( totosc[7]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }$\\overline{\\nu}_e$-like & \n")
	table.write(str("%.3f" % ( ccnue[8]/totosc[8])) + " & " + str("%.3f" % ( ccnueb[8]/totosc[8])) + " & " + str("%.3f" % ( ccmu[8]/totosc[8])) + " & " + str("%.3f" % ( cctau[8]/totosc[8])) + " & " + str("%.3f" % ( nc[8]/totosc[8])) +  " & " + str("%.3f" % ( totosc[8]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }$\\mu$-like & \n")
	table.write(str("%.3f" % ( ccnue[9]/totosc[9])) + " & " + str("%.3f" % ( ccnueb[9]/totosc[9])) + " & " + str("%.3f" % ( ccmu[9]/totosc[9])) + " & " + str("%.3f" % ( cctau[9]/totosc[9])) + " & " + str("%.3f" % ( nc[9]/totosc[9])) +  " & " + str("%.3f" % ( totosc[9]/tot)) + " \\\\\n")
	table.write("\\text{  }Multi-Ring\\\\\n")
	table.write("\\text{  }\\text{  }$\\nu_e$-like & \n")
	table.write(str("%.3f" % ( ccnue[10]/totosc[10])) + " & " + str("%.3f" % ( ccnueb[10]/totosc[10])) + " & " + str("%.3f" % ( ccmu[10]/totosc[10])) + " & " + str("%.3f" % ( cctau[10]/totosc[10])) + " & " + str("%.3f" % ( nc[10]/totosc[10])) +  " & " + str("%.3f" % ( totosc[10]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }$\\overline{\\nu}_e$-like & \n")
	table.write(str("%.3f" % ( ccnue[11]/totosc[11])) + " & " + str("%.3f" % ( ccnueb[11]/totosc[11])) + " & " + str("%.3f" % ( ccmu[11]/totosc[11])) + " & " + str("%.3f" % ( cctau[11]/totosc[11])) + " & " + str("%.3f" % ( nc[11]/totosc[11])) +  " & " + str("%.3f" % ( totosc[11]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }$\\mu$-like & \n")
	table.write(str("%.3f" % ( ccnue[12]/totosc[12])) + " & " + str("%.3f" % ( ccnueb[12]/totosc[12])) + " & " + str("%.3f" % ( ccmu[12]/totosc[12])) + " & " + str("%.3f" % ( cctau[12]/totosc[12])) + " & " + str("%.3f" % ( nc[12]/totosc[12])) +  " & " + str("%.3f" % ( totosc[12]/tot)) + " \\\\\n")
	table.write("\\text{  }\\text{  }Other & \n")
	table.write(str("%.3f" % ( ccnue[13]/totosc[13])) + " & " + str("%.3f" % ( ccnueb[13]/totosc[13])) + " & " + str("%.3f" % ( ccmu[13]/totosc[13])) + " & " + str("%.3f" % ( cctau[13]/totosc[13])) + " & " + str("%.3f" % ( nc[13]/totosc[13])) +  " & " + str("%.3f" % ( totosc[13]/tot)) + " \\\\\n")
	table.write("\\vspace{0.1cm}\n")
	table.write("\\textbf{Partially Contained}\\\\\n")
	table.write("\\text{  }PC Stopping\\\\\n")
	table.write(str("%.3f" % ( ccnue[14]/totosc[14])) + " & " + str("%.3f" % ( ccnueb[14]/totosc[14])) + " & " + str("%.3f" % ( ccmu[14]/totosc[14])) + " & " + str("%.3f" % ( cctau[14]/totosc[14])) + " & " + str("%.3f" % ( nc[14]/totosc[14])) +  " & " + str("%.3f" % ( totosc[14]/tot)) + " \\\\\n")
	table.write("\\text{  }PC Through-going\\\\\n")
	table.write(str("%.3f" % ( ccnue[15]/totosc[15])) + " & " + str("%.3f" % ( ccnueb[15]/totosc[15])) + " & " + str("%.3f" % ( ccmu[15]/totosc[15])) + " & " + str("%.3f" % ( cctau[15]/totosc[15])) + " & " + str("%.3f" % ( nc[15]/totosc[15])) +  " & " + str("%.3f" % ( totosc[15]/tot)) + " \\\\\n")
	table.write("\\hline\n")
	table.write("\\end{tabular}\n")
	table.write("%\\caption{Sample  purity  broken  down  by  neutrino  flavor  assuming  neutrino  oscillations. These correspond to the fully-contained (FC) samples and from the implemented Super-Kamiokande reconstruction emulation using atmospheric neutrinos generated by GENIE, \\cite{}.}\n")
	table.write("\\label{table:sk_emulated}\n")
	table.write("\\end{table*}\n")
