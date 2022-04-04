# Milankovitch Cycles
#
# Author: Joseph A'Hearn
# Created 04/23/2019
#
# This program attempts to recreate Earth's orbital parameters

import numpy as np 
import matplotlib.pyplot as plt 
import plot_assistant as pa 

def generate_data(t):
	e = 0.005 + np.absolute(0.055 * ((np.sin(1.5E-01 * np.pi * t))**2) * np.sin(2.5E-03 * np.pi * t))
	obliquity = 0.406 + (0.015 * np.sin(4.5E-05 * np.pi * t) * np.sin(5.0E-03 * np.pi * t))
	omega = np.sin(1.5E-04 * np.pi * t)
	return e, obliquity, omega

def make_plot(t, e, obliquity, omega, black=True):
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(8,10))
	else:
		fig = plt.figure(figsize=(8,10))
	rowspan = 5
	colspan = 15
	ax1 = plt.subplot2grid((15,15),(0, 0), rowspan=rowspan, colspan=colspan)
	ax1.set_title('Milankovitch cycles', fontsize=24)
	ax1.set_ylabel('eccentricity', fontsize=22)
	ax1.set_xlim([t[0],t[-1] / 1.0E+06])
	ax1.set_ylim([0,0.06])
	ax1.plot(t / 1.0E+06, e)

	plt.setp(ax1.get_xticklabels(), visible=False) # make these tick labels invisible
	ax2 = plt.subplot2grid((15,15),(5, 0), rowspan=rowspan, colspan=colspan)
	ax2.set_ylabel('obliquity [radians]', fontsize=22)
	ax2.set_xlim([t[0],t[-1] / 1.0E+06])
	ax2.set_ylim([0.38,0.43])
	plt.setp(ax2.get_xticklabels(), visible=False) # make these tick labels invisible
	ax2.plot(t / 1.0E+06, obliquity)

	ax3 = plt.subplot2grid((15,15),(10, 0), rowspan=rowspan, colspan=colspan)
	ax3.set_xlabel('t [Ma]', fontsize=22)
	ax3.set_ylabel(r'$e$' + 'sin' + r'$\omega$', fontsize=22)
	ax3.set_xlim([t[0],t[-1] / 1.0E+06])
	ax3.set_ylim([-0.06,0.06])
	ax3.plot(t / 1.0E+06, e * np.sin(omega))

	
	ax1.invert_xaxis()
	ax2.invert_xaxis()
	ax3.invert_xaxis()
	plt.tight_layout()

	if(black==True):
		ax1.spines['top'].set_visible(False)
		ax1.spines['right'].set_visible(False)
		ax1.spines['bottom'].set_linewidth(0.5)
		ax1.spines['left'].set_linewidth(0.5)
		ax1.spines['bottom'].set_color('white')
		ax1.spines['left'].set_color('white')
		ax1.title.set_color('white')
		ax1.yaxis.label.set_color('white')
		ax1.xaxis.label.set_color('white')
		ax1.tick_params(axis='x', colors='white')
		ax1.tick_params(axis='y', colors='white')
		ax1.tick_params(axis='both', direction='in')
		ax1.get_xaxis().tick_bottom()
		ax1.get_yaxis().tick_left()
		ax2.spines['top'].set_visible(False)
		ax2.spines['right'].set_visible(False)
		ax2.spines['bottom'].set_linewidth(0.5)
		ax2.spines['left'].set_linewidth(0.5)
		ax2.spines['bottom'].set_color('white')
		ax2.spines['left'].set_color('white')
		ax2.title.set_color('white')
		ax2.yaxis.label.set_color('white')
		ax2.xaxis.label.set_color('white')
		ax2.tick_params(axis='x', colors='white')
		ax2.tick_params(axis='y', colors='white')
		ax2.tick_params(axis='both', direction='in')
		ax2.get_xaxis().tick_bottom()
		ax2.get_yaxis().tick_left()
		ax3.spines['top'].set_visible(False)
		ax3.spines['right'].set_visible(False)
		ax3.spines['bottom'].set_linewidth(0.5)
		ax3.spines['left'].set_linewidth(0.5)
		ax3.spines['bottom'].set_color('white')
		ax3.spines['left'].set_color('white')
		ax3.title.set_color('white')
		ax3.yaxis.label.set_color('white')
		ax3.xaxis.label.set_color('white')
		ax3.tick_params(axis='x', colors='white')
		ax3.tick_params(axis='y', colors='white')
		ax3.tick_params(axis='both', direction='in')
		ax3.get_xaxis().tick_bottom()
		ax3.get_yaxis().tick_left()
		ax1.figure.savefig('earth_milankovitch_cycles_black.png', facecolor=fig.get_facecolor(), transparent=True)
	else:
		ax1.figure.savefig('earth_milankovitch_cycles_white.png')
	plt.clf()

def calculate_precession_period(TE, aM2RE):
	Omega = 2 * np.pi / (TE * 3600)
	aM = aM2RE * 3.84399E+08 / 60
	J2 = 0.00335
	cosPsi = 0.915
	mu_m = 0.0123
	G = 6.67384E-11
	ME = 5.974E+24
	MS = 1.9885E+30
	aE = 1.495978707E+11
	nm = np.sqrt(G * ME / (aM**3))
	ns = np.sqrt(G * MS / (aE**3))
	return 2 * Omega / (3 * J2 * cosPsi * ns * (1 + (mu_m * ((nm / ns)**2))))

def main(Nt=10000):
	t = np.linspace(0, 1.0E+06, Nt) # units of years
	e, obliquity, omega = generate_data(t)
	make_plot(t, e, obliquity, omega)

#main()
#TE = 24
#aM2RE = 60
#TE = 10
#aM2RE = 30
#TE = 23.99926944
#aM2RE = 60 * 0.9999937565
TE = 22
aM2RE = 58

print(calculate_precession_period(TE, aM2RE))