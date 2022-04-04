# Tidal Evolution
#
# Author: Joseph A'Hearn
# Created 02/28/2019
#
# This program makes tidal evolution plots 

# Zhang 2009 assumes a Q of 18000

import numpy as np 
import matplotlib.pyplot as plt 
import plot_assistant as pa 
import useful as uf 
from fractions import Fraction 
from matplotlib.ticker import MultipleLocator

def load_data(filename):
	with open(filename) as f:
		content = f.readlines()
	content = [x.strip() for x in content]
	central_body = content[0][14:]

	N = uf.rawcount(filename) - 3
	labels = []
	ms = np.zeros(N)
	af = np.zeros(N)
	for i in range(N):
		line_with_body = content[4+i]
		name = line_with_body[0:9].strip().split("\t")[0] 
		stuff = line_with_body.strip().split("\t")
		labels.append(stuff[0])
		ms[i] = float(stuff[1])
		af[i] = float(stuff[2])
	return central_body, N, labels, ms, af  

def load_central_body(central_body):
	if(central_body == 'Saturn'):
		k2 = 0.390 	# Lainey et al 2017
		Q = 2000	# Lainey et al 2017
		C = 6.0330E+07 # planetary radius in m from Cooper et al 2015
		M = 5.6836E+26 # mass in kg from Cooper et al 2015
		#Omega = 1.631E-04 # rad/s 10.7 hr
		Omega = 1.65269E-04 # rad/s 10.56056 hr from Mankovich et al 2019
	elif(central_body == 'Uranus'):
		k2 = 0.104 # Lainey 2016 cites Gavrilov & Zharkov (1977)
		Q = 500 # higher limit according to Lainey 2016	# the range estimate by Tittemore & Wisdom 1990 was 11,000 < Q < 39,000
		C = 2.5559E+07 
		M = 8.6810E+25
		Omega = 1.015E-04 # 17.2 hr
	elif(central_body == 'Neptune'):
		k2 = 0.127 # Lainey 2016 cites Gavrilov & Zharkov (1977); 0.41 (Bursa 1992); 0.39 (Banfield & Murray 1992)
		Q = 12000 # the range estimate of Showalter et al 2019 was 12,000 < Q < 330,000	(This is from Banfield & Murray 1992). Zhang & Hamilton 2008 suggest 9,000 < Q < 36,000
		C = 2.4764E+07 # Archinal et al. 2018
		M = 1.024126E+26 # Jacobson 2009
		Omega = 1.084E-04 # 16.1 hr; Lecacheux et al. 1993
	elif(central_body == 'Earth'):
		k2 = 0.302 # Lerch et al 1992, Lemoine et al 1998
		Q = 280 # Ray et al 2001  
		C = 6.3781E+06
		M = 5.97237E+24
		Omega = 7.2921159E-05 # 23.93446945 hr, i.e 23 hr, 56 m, 4.09 s
	elif(central_body == 'Sun'):
		k2 = 1 # unknown
		Q = 1 # unknown
		C = 6.957E+08 
		M = 1.98955E+30
		Omega = 2.972E-06 # 24.27 days = 587.28 hr

	return k2, Q, C, M, Omega

def compute_a(t, a_sr, a0, m, name, k2=0.390, Q=2000, C=6.0268E+07, M=5.6846E+26, G=6.673E-11, custom_Q=True, future=False): # SSD 4.213 with a little bit of algebra, default values are for Saturn and its satellites in the Mimas-Enceladus regime
	# C is planetary radius
	if(custom_Q == True):
		if(name == 'Titan'):
			Q = 100 # Lainey et al 2020
		elif(name == 'Rhea'):
			Q = 300 # Lainey et al 2020
		elif(name == 'Tethys'):
			Q = 3000 # Lainey et al 2020

	# This code works but I'm experimenting below with try and except for cases when a satellite gets to a_sr
	#if(a0 > a_sr):
	#	a = a0 * ((1 - ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t * 8.64E+04 * 365.25))**(2/13))
	#else:
	#	a = a0 * ((1 + ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t * 8.64E+04 * 365.25))**(2/13))

	if(a0 > a_sr):
		if(future == True):
			if(name == 'Triton'): # retrograde orbit, looking into the future
				try:
					a = a0 * ((1 - ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t * 8.64E+04 * 365.25))**(2/13))
				except RuntimeWarning:
					a = np.zeros(len(t))
					for j in range(len(t)):
						try: 
							a[j] = a0 * ((1 - ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t[j] * 8.64E+04 * 365.25))**(2/13))
						except:
							a[j] = a_sr
			else:
				try:
					a = a0 * ((1 + ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t * 8.64E+04 * 365.25))**(2/13))
				except RuntimeWarning:
					a = np.zeros(len(t))
					for j in range(len(t)):
						try: 
							a[j] = a0 * ((1 + ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t[j] * 8.64E+04 * 365.25))**(2/13))
						except:
							a[j] = a_sr
		elif(name == 'Triton'): # retrograde orbit, looking into the past (no try and except because at present it's beyond a_sr)
			print('Triton here')
			a = a0 * ((1 + ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t * 8.64E+04 * 365.25))**(2/13))
		else:
			try:
				a = a0 * ((1 - ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t * 8.64E+04 * 365.25))**(2/13))
			except RuntimeWarning:
				a = np.zeros(len(t))
				for j in range(len(t)):
					try: 
						a[j] = a0 * ((1 - ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t[j] * 8.64E+04 * 365.25))**(2/13))
					except:
						a[j] = a_sr
	else:
		if(future == True):
			a = a0 * ((1 - ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t * 8.64E+04 * 365.25))**(2/13))
		else:
			a = a0 * ((1 + ((39 * k2 / (2 * Q * (a0**(13/2)))) * np.sqrt(G/M) * (C**5) * m * t * 8.64E+04 * 365.25))**(2/13))
	return a 

def find_resonances(t, a, N, Nt, mu=6.673E-11*5.6846E+26, maxrange=12, maxorder=2):
	res_o = []
	res_k = []
	res_l = []
	res_i = []
	res_j = []
	n = np.zeros((N, Nt))
	for i in range(N):
		n[i] = 1 / (np.sqrt(a[i])**3 / mu) # units of rad/s

	for i in range(N):
		for j in range(N):
			if(i < j):
				for k in range(1, maxrange):
					for l in range(1, maxrange):
						if(0 < np.absolute(k - l) < (maxorder + 1)):
							phi_dot = k * n[i] - l * n[j]
							for o in range(1, Nt):
								if(np.sign(phi_dot[o]) != np.sign(phi_dot[o-1])):
									x = Fraction(k, l)
									if(x.numerator != k):
										continue
									else:
										res_o.append(o)
										res_k.append(k)
										res_l.append(l)
										res_i.append(i)
										res_j.append(j)

										#print(str(k) + ':' + str(l) + ' resonance')
										#for z in range(N):
										#	print((a[z][o] + a[z][o-1]) / 2)
										#print((a[j][o] + a[j][o-1]) / 2)
										#print((a[i][o] + a[i][o-1]) / 2)

	return res_o, res_k, res_l, res_i, res_j

def print_important_data(t, a):
	print('Final: ' + str(t[0] / 1.0E+06) + ' Ma, ' + str(a[0][0] / 1.0E+03) + ' km')
	print('Initial: ' + str(t[-1] / 1.0E+06) + ' Ma, ' + str(a[0][-1] / 1.0E+03) + ' km')
	for i in range(1, len(t)):
		if((t[i-1] < 0.8E+06) and (t[i] > 0.8E+06)):
			print('MPT: ' + str(t[i] / 1.0E+06) + ' Ma, ' + str(a[0][i] / 1.0E+03) + ' km')
		if((t[i-1] < 0.9E+06) and (t[i] > 0.9E+06)):
			print('MPT: ' + str(t[i] / 1.0E+06) + ' Ma, ' + str(a[0][i] / 1.0E+03) + ' km')

def make_plot(central_body, N, t, a_sr, a_rrl, a_frl, a_Lind, a, rho_s, labels, res_o, res_k, res_l, res_i, res_j, filename, black=False, future=False, show_Roche_limits=True):
	fig, ax = pa.initialize_plot(black)

	# determine time unit: yr, ka, Ma, or Ga
	yr = 1
	ka = 1.0E+03
	Ma = 1.0E+06
	Ga = 1.0E+09

	time_unit = Ma

	if(black == True):
		c = 'w'
	else:
		c = 'k'

	if(central_body == 'Saturn'):
		#fig, ax = pa.title_and_axes(fig, ax, 'Tidal evolution', 't [Ma]', 'a [10' + r'$^5$' + ' km]', t[0], t[-1] / 1.0E+06, np.mean(a[0]) * 0.8 / 1.0E+08, np.mean(a[-1]) * 1.1 / 1.0E+08)
		#fig, ax = pa.title_and_axes(fig, ax, '', 't [Ma]', 'a [10' + r'$^5$' + ' km]', t[0], t[-1] / time_unit, 1.7, 2.4)
		if(filename == 'pallene_vicinity'):
			a_min = 1.7
			#a_min = 1.2
			#a_max = 2.4
			a_max = 4.0
		elif(filename == 'titan_vicinity'):
			a_min = 5
			a_max = 15
		fig, ax = pa.title_and_axes(fig, ax, '', 't [Ma]', 'a [10' + r'$^5$' + ' km]', t[0], t[-1] / time_unit, a_min, a_max)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		if(filename == 'pallene_vicinity'):
			ax.set_xticks([0, 50, 100, 150])
			ax.set_yticks([1.8, 2.0, 2.2, 2.4])
			#ax.set_xticks([0, 100, 200, 300, 400])
			#ax.set_yticks([1.2, 1.6, 2.0, 2.4])
			ax.tick_params(which='minor', length=4)
			ax.tick_params(which='major', length=8)
			ax.tick_params(which='both', top=False, right=False)
			ax.xaxis.set_minor_locator(MultipleLocator(10))
			ax.yaxis.set_minor_locator(MultipleLocator(0.05))
	elif(central_body == 'Neptune'):
		time_unit = Ga
		#time_unit = Ma
		#a_min = 0.24764 # Neptune equatorial radius
		a_min = 0.46
		#a_max = 0.75 # 0.8 when Larissa is the outermost moon
		a_max = 1.2 # 1.2 when Proteus is included and is the outermost moon
		#a_max = 4.0 # when Triton is included
		if(time_unit == Ga):
			fig, ax = pa.title_and_axes(fig, ax, '', 't [Ga]', 'a [10' + r'$^5$' + ' km]', t[0], t[-1] / time_unit, a_min, a_max)
		else:
			fig, ax = pa.title_and_axes(fig, ax, '', 't [Ma]', 'a [10' + r'$^5$' + ' km]', t[0], t[-1] / time_unit, a_min, a_max)
	elif(central_body == 'Earth'):
		time_unit = Ga
		a_min = 3.75
		a_max = 3.85
		#fig, ax = pa.title_and_axes(fig, ax, 'Tidal evolution', 't [Ma]', 'a [10' + r'$^5$' + ' km]', t[0], t[-1] / 1.0E+06, 3.8439, 3.84405)
		fig, ax = pa.title_and_axes(fig, ax, '', 't [Ga]', 'a [10' + r'$^5$' + ' km]', t[0], t[-1] / time_unit, a_min, a_max)
		#ax.plot([0.9,  0.9], [0, 5], c='r', linestyle='dashed')
		#ax.plot([0.8,  0.8], [0, 5], c='r', linestyle='dashed')
		#ax.text(0.83, 3.84394, 'mid-Pleistocene transition', horizontalalignment='right', verticalalignment='center', rotation='vertical', color=c, fontsize=8)
		#ax.text(1.26, 3.84396, r'$a_0=$' + '{0:.6f}'.format(a[0][-1] / 1.0E+08) + r'$\times 10^5$' + ' km', horizontalalignment='center', verticalalignment='center', color=c, fontsize=9)
		#ax.text(0.25, 3.843995, r'$a_f=$' + '{0:.6f}'.format(a[0][0]  / 1.0E+08) + r'$\times 10^5$' + ' km', horizontalalignment='center', verticalalignment='center', color=c, fontsize=9)
		ax.text(3.26, 3.77, r'$a_0=$' + '{0:.6f}'.format(a[0][-1] / 1.0E+08) + r'$\times 10^5$' + ' km', horizontalalignment='center', verticalalignment='center', color=c, fontsize=9)
		ax.text(0.25, 3.843995, r'$a_f=$' + '{0:.6f}'.format(a[0][0]  / 1.0E+08) + r'$\times 10^5$' + ' km', horizontalalignment='center', verticalalignment='center', color=c, fontsize=9)
		
	for j in range(N):
		if(filename == 'pallene_vicinity'):
			c = ['b', 'g', 'r', 'y', 'm']
			ax.plot(t / time_unit, a[j] / 1.0E+08, linewidth=5.0, color=c[j])
			#ax.text(0.5 * t[-1] / time_unit, (a[j][0] + a[j][-1]) / 2.01E+08, labels[j], verticalalignment='top', weight='bold', fontsize=20)
			if(j == 0):
				#ax.text(0.5 * t[-1] / time_unit, (a[j][0] + a[j][-1]) / 1.88E+08, labels[j], verticalalignment='top', weight='bold', fontsize=20)
				ax.text(0.5 * t[-1] / time_unit, (a[j][0] + a[j][-1]) / 2.00E+08, labels[j], verticalalignment='top', weight='bold', fontsize=20)
			elif(j == 1):
				#ax.text(0.9 * t[-1] / time_unit, (a[j][0] + a[j][-1]) / 2.03E+08, labels[j], verticalalignment='top', weight='bold', fontsize=20)
				ax.text(0.5 * t[-1] / time_unit, (a[j][0] + a[j][-1]) / 2.01E+08, labels[j], verticalalignment='top', weight='bold', fontsize=20)
			elif(j == 2):
				#ax.text(0.98 * t[-1] / time_unit, (a[j][0] + a[j][-1]) / 1.95E+08, labels[j], verticalalignment='top', weight='bold', fontsize=20)
				ax.text(0.9 * t[-1] / time_unit, (a[j][0] + a[j][-1]) / 1.96E+08, labels[j], verticalalignment='top', weight='bold', fontsize=20)
		else:
			ax.plot(t / time_unit, a[j] / 1.0E+08, linewidth=5.0)
			ax.text(0.5 * t[-1] / time_unit, (a[j][0] + a[j][-1]) / 2.01E+08, labels[j], horizontalalignment='center' , verticalalignment='top', weight='bold')
		#ax.plot(t / time_unit, a[j] / 1.0E+08, label=labels[j])
	if(a_min < (a_sr / 1.0E+08) < a_max):
		ax.plot([t[0], t[-1] / time_unit], [a_sr / 1.0E+08, a_sr / 1.0E+08], linestyle='dashed', color=c)
		ax.text(0.5 * t[-1] / time_unit, a_sr / 1.0E+08, 'synchronous orbit', verticalalignment='top', horizontalalignment='center')
		
	for i in range(len(res_o)):
		ax.plot([t[res_o[i]] / time_unit,  t[res_o[i]] / time_unit], [a[res_i[i]][res_o[i]] / 1.0E+08, a[res_j[i]][res_o[i]] / 1.0E+08], c=c)
		ax.text(t[res_o[i]] / time_unit, (a[res_i[i]][res_o[i]] + a[res_j[i]][res_o[i]]) / 2.0E+08, str(res_k[i]) + ':' + str(res_l[i]), 
			horizontalalignment='right', verticalalignment='center', rotation='vertical', color=c, fontsize=12)

	if(show_Roche_limits == True):
		if(a_min < (a_rrl / 1.0E+08) < a_max):
			ax.plot([t[0], t[-1] / time_unit], [a_rrl / 1.0E+08, a_rrl / 1.0E+08], linestyle='dashed', color=c)
			ax.text(0.5 * t[-1] / time_unit, a_rrl / 1.0E+08, 'rigid Roche limit for ' + r'$\rho=$' + str(np.round(rho_s / 1000, decimals=1)) + ' g/cm' + r'$^3$', verticalalignment='top', horizontalalignment='center')
		if(a_min < (a_frl / 1.0E+08) < a_max):
			ax.plot([t[0], t[-1] / time_unit], [a_frl / 1.0E+08, a_frl / 1.0E+08], linestyle='dashed', color=c)
			ax.text(0.5 * t[-1] / time_unit, a_frl / 1.0E+08, 'fluid Roche limit', verticalalignment='top', horizontalalignment='center')
		if(a_min < (a_Lind / 1.0E+08) < a_max):
			ax.plot([t[0], t[-1] / time_unit], [a_Lind / 1.0E+08, a_Lind / 1.0E+08], linestyle='dashed', color=c)
			ax.text(0.5 * t[-1] / time_unit, a_Lind / 1.0E+08, 'maximum orbit reachable via Lindblad torques', verticalalignment='top', horizontalalignment='center')
		
	if(future == False):
		ax.invert_xaxis()
	#plt.tight_layout()

	#ax.legend(loc=(0.01, 0.94), fontsize=7)
	for_paper = False
	if(for_paper == True):
		ax.set_xlabel(r'$t$' + ' [Ma]', fontsize=20)
		ax.set_ylabel(r'$a$' + ' [10' + r'$^5$' + ' km]', fontsize=20)
		#plt.rc('xtick', labelsize=26)    # fontsize of the tick labels
		#plt.rc('ytick', labelsize=26)    # fontsize of the tick labels
		ax.tick_params(axis='both', labelsize=18)
		ax.tick_params(which='both', axis='both', direction='in')
		ax.figure.savefig(filename + "_orbital_evolution.pdf")
		plt.clf()
	else:
		pa.save_and_clear_plot(fig, ax, filename=filename + '_orbital_evolution', black=black, dpi=300)

def set_to_empty():
	return [], [], [], [], []


def SSD4_150(Omega, G, M, ms, af, a0, I):
	return Omega + ((np.sqrt(G) / I) * (((M**(3/2)) * ms) / (M + ms)) * (np.sqrt(af) - np.sqrt(a0))) 

def main(filename, yrs=1.0E+06, locate_resonances=False, Nt=10**5, G=6.673E-11, black=False, future=False):
	central_body, N, labels, ms, af = load_data(filename + '.txt')
	k2, Q, C, M, Omega = load_central_body(central_body)
	a_sr = (G * M / (Omega**2))**(1/3) # synchronous rotation
	rho_s = 1200
	a_rrl = ((9 * M) / (4 * np.pi * rho_s))**(1/3) # solid Roche limit 
	rho = (3 * M) / (4 * np.pi * (C**3))
	a_frl = 2.456 * C * ((rho / rho_s)**(1/3)) # fluid Roche limit
	a_Lind = (4**(1/3)) * a_frl # maximum orbit a satellite may migrate to via Lindblad torques
	a = np.zeros((N, Nt))
	t = np.linspace(0, yrs, Nt) # units of years
	for i in range(N):
		a[i] = compute_a(t, a_sr, af[i], ms[i], labels[i], k2, Q, C, M, G, False, future)
		#print(a[i][0], a[i][-1])
		#print(100 * (a[i][0] - a[i][-1]) / (a[i][0]))
		print(a[i][-1])
	if(locate_resonances == True):
		res_o, res_k, res_l, res_i, res_j = find_resonances(t, a, N, Nt, G * M)
	else:
		res_o, res_k, res_l, res_i, res_j = set_to_empty()
	make_plot(central_body, N, t, a_sr, a_rrl, a_frl, a_Lind, a, rho_s, labels, res_o, res_k, res_l, res_i, res_j, filename, black, future)
	#print_important_data(t, a)
	#Omega_initial = SSD4_150(Omega, G, M, ms[0], a[0][0], a[0][-1], 8.04E+37)
	#nomega_factor = np.sqrt(G * M) * np.absolute((Omega / ((a[0][0])**(3/2))) - (Omega_initial / ((a[0][-1])**(3/2))))
	#print(2 * np.pi / (Omega_initial * 3600))
	#print(nomega_factor)
	#print(39 * k2 * np.sqrt(G / M) * (C**5) * 4.2E+18 / (2 * Q * (7.3548E+07**(13/2))))

#main('pallene_vicinity', 0.5E+08, True)
#main('pallene_vicinity', 1.75E+08, True)
#main('pallene_vicinity', 5.0E+06, True)
#main('pallene_vicinity', 4.2E+08, True)
#main('titan_vicinity', 2.00E+09, True)
#main('hippocamp_vicinity')
#main('galatea_vicinity', 0.5E+09, True)
# for the 5:4 resonance, 0.9308 was not far back enough, but 0.9309 was
# for the 11:9 resonance, 0.6568 was not far back enough, but 0.6569 was
# for the 6:5 resonance, 0.46145 was not far back enough, but 0.4615 was 
main('galatea_vicinity', 0.9309E+09, True, future=False) 
#main('earthmoo')
#main('solar_system')
#main('D68_vicinity', 1.0E+06)
# Larissa	4.2E+18	7.3548E+07
# Galatea	2.12E+18	6.1953E+07
