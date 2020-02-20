#Basic stats, sliding window, SLIm ms output
#adding divergence to this
#python statistics_bigwindow_mystats_dpgp3.py -cutoff numbp50 -side 5p -state derived -noncodingLen 4000 -consCutoff 0.8
from __future__ import print_function
from libsequence.polytable import SimData
from libsequence.summstats import PolySIM
from libsequence.windows import Windows
from libsequence.summstats import ld
import sys
import pandas
import math
import argparse
import statistics

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
#parser.add_argument('-gen_post_burnin', dest = 'gen_post_burnin', action='store', nargs = 1, default = 100, type = int, choices = range(1,10000), help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-cutoff', dest = 'cutoff', action='store', nargs = 1, default = "numbp50", type = str, help = 'numbp50')#250 bp for small, 5000 bp for big
parser.add_argument('-side', dest = 'side', action='store', nargs = 1, type = str, help = '5 or 3 prime')
parser.add_argument('-state', dest = 'state', action='store', nargs = 1, type = str, help = 'derived or major')
parser.add_argument('-consCutoff', dest = 'consCutoff', action='store', nargs = 1, type = str, help = '0.0 or 0.8 or 0.5')
parser.add_argument('-noncodingLen', dest = 'noncodingLen', action='store', nargs = 1, type = int, choices = range(1,10001), help = 'noncoding length in bp')
args = parser.parse_args()
cutoff =  args.cutoff[0]
noncoding_size =  args.noncodingLen[0]
side = args.side[0]
state = args.state[0]
cons_cutoff = args.consCutoff[0]

#defined constants:
tot_genes = 94 #total number of replicates
#get length of coding from file:
f_sing = open("/Users/paruljohri/BgsDfeDemo/Drosophila/annotation/sing_exon_500_2000_4000_prop.txt", 'r')
d_codingsize = {}
l_genes = []
for line in f_sing:
	line1 = line.strip('\n')
	line2 = line1.split('\t')
	if line2[0] != "gene":
		d_codingsize[line2[0]] = line2[3]
		l_genes.append(line2[0])
f_sing.close()

#get length of noncoding region for each gene (this is different because of masking using phastcons
f_inter_mis =  open("/Users/paruljohri/BgsDfeDemo/Drosophila/dpgp3_intergenic/InsectCons_" + str(cons_cutoff) + "/noncoding_missing_dpgp3.txt", 'r')
d_inter_missing = {}
for line in f_inter_mis:
	line1 = line.strip('\n')
	line2 = line1.split('\t')
	if line2[0] != "gene":
		if line2[1] == side:
			d_inter_missing[line2[0]] = line2[2]
f_inter_mis.close()

#read ms file:
def read_subset_ms(f_ms, start, end):
	l_Pos = [] #list of positions of SNPs
	d_tmp = {}
	l_int_posn = []
	l_missing_rows = []
	row_num = 1
	for line in f_ms:#record relevant positions and missing rows 
		line1 = line.strip('\n')
		if "positions" in line1:
			line2 = line1.split()
			i = 0
			for x in line2:
				if "position" not in x:
					if (float(x) > float(start)) and (float(x) <= float(end)):
						l_Pos.append(float(x))
						d_tmp[str(i)] = ""
						l_int_posn.append(i)
					i = i + 1
		elif "//" not in line and "segsites" not in line:#check for missing rows
			s_gt = ""
			for j in l_int_posn:#moving across columns
				s_gt = s_gt + line1[j]#s_gt represent a row of genotypes
			non_missing_data = len(s_gt) - s_gt.count("N")
			#print (non_missing_data)
			if non_missing_data < 5:
				l_missing_rows.append(row_num)
		row_num = row_num + 1
	f_ms.seek(0)
	row_num = 1
	for line in f_ms:
		line1 = line.strip('\n')
		if "//" not in line and "segsites" not in line and "positions" not in line:
			if row_num not in l_missing_rows:
				for j in l_int_posn:#moving across columns
					d_tmp[str(j)] = d_tmp[str(j)] + line1[j]
		row_num = row_num + 1
	return (l_Pos, d_tmp, l_int_posn)
	
		
def read_fixed_mutations(f_fixed, ChrLen):
    d_subs = {}
    for line in f_fixed:
        line1 = line.strip('\n')
        posn = float(line1)/float(ChrLen)
        d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs #return a dictionary with base positions as keys and number of fixed substitutions as values

def read_subset_fixed_mutations(f_fixed, start, end, ChrLen):
    d_subs = {}
    for line in f_fixed:
        line1 = line.strip('\n')
        posn = float(line1)/float(ChrLen)
        if posn > float(start) and posn < float(end):
        	d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs #return a dictionary with base positions as keys and number of fixed substitutions as values


def avg_divergence_win(d_subs, start, end):
	s_sum = 0
	for posn in d_subs.keys():
		if float(posn) < float(end) and float(posn) > float(start):
			s_sum = s_sum + 1
	return s_sum

def get_S(l_int_Pos, d_tmp):
	S = 0
	for posn in l_int_Pos:
		s_col = d_tmp[str(posn)]
		count_0 = s_col.count("0")
		count_1 = s_col.count("1")
		if count_0 > 0 and count_1 > 0:
			S = S + 1
	return S

def get_tajimasd(l_int_Pos, l_Pos, d_tmp):
	S = get_S(l_int_Pos, l_Pos, d_tmp)
	pi = get_thetapi(l_int_Pos, l_Pos, d_tmp)
	l_td = []
	for posn in l_int_Pos:
		s_col = d_tmp[str(posn)]
		n = s_col.count("0") + s_col.count("1")
		if n >= 2:
			a1, a2 = 0.0, 0.0
			i = 1
			while i <= n-1:
				a1 = float(a1) + 1.0/float(i)
				a2 = float(a2) + 1.0/(float(i)*float(i))
				i = i + 1
			b1 = float(n+1)/float(3.0*float(n-1))
			b2 = float(2.0*(float(n*n)+n+3))/float(9.0*n*float(n-1))
			c1 = float(b1) - (1.0/float(a1))
			c2 = b2 - (float(n+2)/(a1*float(n))) + (a2/float(a1*a1))
			e1 = c1/float(a1)
			e2 = float(c2)/float((a1*a1) + a2)
	return (l_td)

def get_D(l_int_Pos, d_tmp):
	l_D = []
	l_rsq = []
	l_Dprime = []
	i = 0
	while i < len(l_int_Pos):
		j = i + 1
		while j < len(l_int_Pos):
			site1 = l_int_Pos[i]
			site2 = l_int_Pos[j]
			col1 = d_tmp[str(site1)]
			col2 = d_tmp[str(site2)]
			#convert columns to major and minor allele:
			col1_mod, col2_mod = "", ""
			if col1.count("1") > col1.count("0"):
				for x in col1:
					if x == "0":
						col1_mod = col1_mod + "1"
					elif x == "1":
						col1_mod = col1_mod + "0"
					else:
						col1_mod = col1_mod + "N"
			else:
				col1_mod = col1
			if col2.count("1") > col2.count("0"):
				for x in col2:
					if x == "0":
						col2_mod = col2_mod + "1"
					elif x == "1":
						col2_mod = col2_mod + "0"
					else:
						col2_mod = col2_mod + "N"
			else:
				col2_mod = col2
			#subsample to throw missing data:
			c_00, c_01, c_10, c_11 = 0, 0, 0, 0
			c_0_site1, c_1_site1, c_0_site2, c_1_site2 = 0, 0, 0, 0
			s_tot = 0
			m = 0
			while m < len(col1_mod):
				if col1_mod[m] != "N" and col2_mod[m] != "N":
					s_tot = s_tot + 1
					c_0_site1 = c_0_site1 + col1_mod[m].count("0")
					#c_1_site1 = c_1_site1 + col1[m].count("1")
					c_0_site2 = c_0_site2 + col2_mod[m].count("0")
					#c_1_site2 = c_1_site2 + col2[m].count("1")
					s_hap = col1_mod[m] + col2_mod[m]
					c_00 = c_00 + s_hap.count("00")
					#c_01 = c_01 + s_hap.count("01")
					#c_10 = c_10 + s_hap.count("10")
					#c_11 = c_11 + s_hap.count("11")
				m = m + 1
			if s_tot > 0:
				D_00 = (c_00/float(s_tot)) - (c_0_site1/float(s_tot))*(c_0_site2/float(s_tot))
				#D_01 = (c_01/float(s_tot)) - (c_0_site1/float(s_tot))*(c_1_site2/float(s_tot))
				#D_10 = (c_10/float(s_tot)) - (c_1_site1/float(s_tot))*(c_0_site2/float(s_tot))
				#D_11 = (c_11/float(s_tot)) - (c_1_site1/float(s_tot))*(c_1_site2/float(s_tot))
				l_D.append(D_00)
				pA = float(c_0_site1)/float(s_tot)
				pB = float(c_0_site2)/float(s_tot)
				if pA != 0.0 and pA != 1.0 and pB != 0.0 and pB!=1.0:
					rsq = (D_00*D_00)/(pA*(1.0-pA)*pB*(1.0-pB))
					l_rsq.append(rsq)
				if D_00 < 0:
					Dmax = max([(-1.0*pA*pB), (-1.0*(1-pA)*(1-pB))])
				elif D_00 > 0:
					Dmax = min([(pA*(1.0-pB)),((1.0-pA)*pB)])
				else:
					Dmax = 0.0
				if Dmax != 0.0:
					Dprime = D_00/float(Dmax)
					l_Dprime.append(Dprime)
			j = j + 1
		i = i + 1
	return (statistics.mean(l_D), statistics.mean(l_rsq), statistics.mean(l_Dprime))

#get numbpLinked:
f_rec = open("/Users/paruljohri/BgsDfeDemo/Drosophila/dpgp3_statistics/InsectCons_" + cons_cutoff + "/SingExon_200_" + side + "_" + state + ".recoverypi", 'r')
for line in f_rec:
	line1 = line.strip('\n')
	line2 = line1.split('\t')
	if line2[0] == "slope":
		i = 0
		for x in line2:
			if x == cutoff:
				col = i
			i = i + 1
	else:
		numbp = line2[col]
if float(numbp) <= noncoding_size:
    startLinked_size = noncoding_size - float(numbp)
else:
	startLinked_size = 0
f_rec.close()

result = open("/Users/paruljohri/BgsDfeDemo/Drosophila/dpgp3_statistics/InsectCons_" + cons_cutoff + "/SingExon_" + side + "_" + state + "_bigwindow_mystats.stats", 'w+')
result.write("simID" + '\t' + "WinType" + '\t' + "WinSize" + '\t' + "S" + '\t' +  "rsq" + '\t' + "D" + '\t' + "Dprime" + '\n')


#go through all simulation replicates and read data into pylibseq format
#addin the option of ignoring some files if they don't exist

s_absent = 0
for gene in l_genes:    
	coding_size = int(d_codingsize[gene])
	chr_len = int(noncoding_size) + coding_size
	start_neu = 2000/float(chr_len)    
	startLinked = float(startLinked_size)/float(chr_len)
	start_func = float(noncoding_size)/float(chr_len)
	end_func = 1.0
	#calculate effective denominator for the intergenic (masked) regions:
	neu_length, link_length = 0, 0
	i = 2000
	while i < int(startLinked_size):
		if d_inter_missing[gene][i] == "0":
			neu_length = neu_length + 1
		i = i + 1
	i = int(startLinked_size)+1
	while i < int(noncoding_size):
		if d_inter_missing[gene][i] == "0":
			link_length = link_length + 1
		i = i + 1
	print ("Reading file:" + gene)
	print ("bins are:" + '\t' + str(start_neu) + '\t' + str(startLinked) + '\t' + str(start_func) + '\t' + str(end_func))
	#try:
	if len("numsim") > 0:
		f_subs = open("/Users/paruljohri/BgsDfeDemo/Drosophila/dpgp3_intergenic/InsectCons_" + cons_cutoff + "/" + gene + "_" + side + ".fixed", 'r')
		d_subs_func = read_subset_fixed_mutations(f_subs, start_func, end_func, chr_len)
		f_subs.seek(0)
		d_subs_link = read_subset_fixed_mutations(f_subs, startLinked, start_func, chr_len)
		f_subs.seek(0)
		d_subs_neu = read_subset_fixed_mutations(f_subs, start_neu, startLinked, chr_len)
		f_subs.close()
		#reading ms file:
		f_ms = open("/Users/paruljohri/BgsDfeDemo/Drosophila/dpgp3_intergenic/InsectCons_" + cons_cutoff + "/" + gene + "_" + side + "_" + state + ".ms", 'r')
		t_ms_func = read_subset_ms(f_ms, start_func, end_func)
		f_ms.seek(0)
		t_ms_link = read_subset_ms(f_ms, startLinked, start_func)
		f_ms.seek(0)
		t_ms_neu = read_subset_ms(f_ms, start_neu, startLinked)
		f_ms.close()
		
		l_Pos_func, l_Pos_link, l_Pos_neu = t_ms_func[0], t_ms_link[0], t_ms_neu[0]
		d_tmp_func, d_tmp_link, d_tmp_neu = t_ms_func[1], t_ms_link[1], t_ms_neu[1]
		l_intPos_func, l_intPos_link, l_intPos_neu = t_ms_func[2], t_ms_link[2], t_ms_neu[2]
		if gene == "ato":
			print (l_Pos_func)
			print (d_tmp_func)
			print (l_intPos_func)
		#print (d_tmp_link)
		#print (d_tmp_link[str(l_intPos_link[0])])
		#print (len(d_tmp_link))
		#print (len(d_tmp_neu))
		#neu:
		l_data = []
		l_Genos_neu = []
		j = 0
		for i in l_intPos_neu:
			if d_tmp_neu[str(i)] != "":
				l_Genos_neu.append(d_tmp_neu[str(i)])
				t_tmp = (l_Pos_neu[j], d_tmp_neu[str(i)])
				l_data.append(t_tmp)
			j = j + 1
		sd_neu = SimData(l_data)
        #link:
		l_data = []
		l_Genos_link = []
		j = 0
		for i in l_intPos_link:
			if d_tmp_link[str(i)] != "":
				l_Genos_link.append(d_tmp_link[str(i)])
				t_tmp = (l_Pos_link[j], d_tmp_link[str(i)])
				l_data.append(t_tmp)
			j = j + 1
		sd_link = SimData(l_data)
		#func:
		l_data = []
		l_Genos_func = []
		j = 0
		for i in l_intPos_func:
			if d_tmp_func[str(i)] != "":
				l_Genos_func.append(d_tmp_func[str(i)])
				t_tmp = (l_Pos_func[j], d_tmp_func[str(i)])
				l_data.append(t_tmp)
			j = j + 1
		sd_func = SimData(l_data)		
		
		#calculating stats:
		div_func, div_link, div_neu = avg_divergence_win(d_subs_func, start_func, end_func), avg_divergence_win(d_subs_link, startLinked, start_func), avg_divergence_win(d_subs_neu, start_neu, startLinked)
		ps_func, ps_link, ps_neu = PolySIM(sd_func), PolySIM(sd_link), PolySIM(sd_neu)
		if sd_func.size() > 0: #To check if there is some data. We can change this cutoff size
			if len(sd_func.pos()) >= 5 and len(d_tmp_func[str(l_intPos_func[0])]) >= 10:
				LD_func = ld(sd_func)
				LDstats_func = pandas.DataFrame(LD_func)
				meanrsq_func = sum(LDstats_func['rsq'])/len(LDstats_func['rsq'])
				meanD_func = sum(LDstats_func['D'])/len(LDstats_func['D'])
				meanDprime_func = sum(LDstats_func['Dprime'])/len(LDstats_func['Dprime'])
				myLD_func = get_D(l_intPos_func, d_tmp_func)
			else:
				LD_func, meanrsq_func, meanD_func, meanDprime_func = "NA", "NA", "NA", "NA"
			result.write(str(gene) + '\t' + str("functional") + '\t' + str((end_func-start_func)*float(chr_len)) + '\t' + str("{:.5f}".format(get_S(l_intPos_func, d_tmp_func))) + '\t' +  str(myLD_func[1]) + '\t' + str(myLD_func[0]) + '\t' + str(myLD_func[2]) + '\n')
		else: #Most rows were empty, so there is no data
			result.write(str(gene) + '\t' + str("functional") + '\t' + str((end_func-start_func)*float(chr_len)) + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\n')
		
		if sd_link.size() > 0: #To check if there is some data. We can change this cutoff size
			if len(sd_link.pos()) >= 5 and len(d_tmp_link[str(l_intPos_link[0])]) >= 10:
				LD_link = ld(sd_link)
				LDstats_link = pandas.DataFrame(LD_link)
				meanrsq_link = sum(LDstats_link['rsq'])/len(LDstats_link['rsq'])
				meanD_link = sum(LDstats_link['D'])/len(LDstats_link['D'])
				meanDprime_link = sum(LDstats_link['Dprime'])/len(LDstats_link['Dprime'])
				myLD_link = get_D(l_intPos_link, d_tmp_link)
			else:
				LD_link, meanrsq_link, meanD_link, meanDprime_link = "NA", "NA", "NA", "NA"
			result.write(str(gene) + '\t' + str("linked") + '\t' + str(link_length) + '\t' + str("{:.5f}".format(get_S(l_intPos_link, d_tmp_link))) + '\t' +  str(myLD_link[1]) + '\t' + str(myLD_link[0]) + '\t' + str(myLD_link[2]) + '\n')
		else:
			result.write(str(gene) + '\t' + str("linked") + '\t' + str(link_length) + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\n')

		if sd_neu.size() > 0: #To check if there is some data. We can change this cutoff size
			if len(sd_neu.pos()) >= 5 and len(d_tmp_neu[str(l_intPos_neu[0])]) >= 10:
				LD_neu = ld(sd_neu)
				LDstats_neu = pandas.DataFrame(LD_neu)
				meanrsq_neu = sum(LDstats_neu['rsq'])/len(LDstats_neu['rsq'])
				meanD_neu = sum(LDstats_neu['D'])/len(LDstats_neu['D'])
				meanDprime_neu = sum(LDstats_neu['Dprime'])/len(LDstats_neu['Dprime'])
				myLD_neu = get_D(l_intPos_neu, d_tmp_neu)
			else:
				LD_neu, meanrsq_neu, meanD_neu, meanDprime_neu = "NA", "NA", "NA", "NA"
			result.write(str(gene) + '\t' + str("neutral") + '\t' + str(neu_length) + '\t' + str("{:.5f}".format(get_S(l_intPos_neu, d_tmp_neu))) + '\t' +  str(myLD_neu[1]) + '\t' + str(myLD_neu[0]) + '\t' + str(myLD_neu[2]) + '\n')
		else:
			result.write(str(gene) + '\t' + str("neutral") + '\t' + str(neu_length) + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\n')

        
		#write results:
        #result.write(str(gene) + '\t' + str("functional") + '\t' + str((end_func-start_func)*float(chr_len)) + '\t' + str("{:.5f}".format(ps_func.thetapi())) + '\t' + str("{:.5f}".format(ps_func.thetaw())) + '\t' + str("{:.5f}".format(ps_func.thetah())) + '\t' + str("{:.5f}".format(ps_func.hprime())) + '\t' + str("{:.5f}".format(ps_func.tajimasd())) + '\t' + str(ps_func.numexternalmutations()) + '\t' + str("{:.5f}".format(ps_func.hapdiv())) + '\t' + str(meanrsq_func) + '\t' + str(meanD_func) + '\t' + str(meanDprime_func) + '\t' + str("{:.5f}".format(div_func)) + '\n')
        #result.write(str(gene) + '\t' + str("linked") + '\t' + str(link_length) + '\t' + str("{:.5f}".format(ps_link.thetapi())) + '\t' + str("{:.5f}".format(ps_link.thetaw())) + '\t' + str("{:.5f}".format(ps_link.thetah())) + '\t' + str("{:.5f}".format(ps_link.hprime())) + '\t' + str("{:.5f}".format(ps_link.tajimasd())) + '\t' + str("{:.5f}".format(ps_link.numexternalmutations())) + '\t' + str("{:.5f}".format(ps_link.hapdiv())) + '\t' + str(meanrsq_link) + '\t' + str(meanD_link) + '\t' + str(meanDprime_link) + '\t' + str("{:.5f}".format(div_link)) + '\n')
        #result.write(str(gene) + '\t' + str("neutral") + '\t' + str(neu_length) + '\t' + str("{:.5f}".format(ps_neu.thetapi())) + '\t' + str("{:.5f}".format(ps_neu.thetaw())) + '\t' + str("{:.5f}".format(ps_neu.thetah())) + '\t' + str("{:.5f}".format(ps_neu.hprime())) + '\t' + str("{:.5f}".format(ps_neu.tajimasd())) + '\t' + str(ps_neu.numexternalmutations()) + '\t' + str("{:.5f}".format(ps_neu.hapdiv())) + '\t' + str(meanrsq_neu) + '\t' + str(meanD_neu) + '\t' + str(meanDprime_neu) + '\t' + str("{:.5f}".format(div_neu)) + '\n')

	#except:
	else:
		s_absent = s_absent + 1
		print ("This file does not exist or cannot be read or is empty")

result.close()
print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






