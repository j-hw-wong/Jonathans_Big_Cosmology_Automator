import os
import csv
import matplotlib
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/3x2pt_pipeline_new/angular_binning/')
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/3x2pt_pipeline_new/gaussian_cl_likelihood/')
#sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/3x2pt_sim')
sys.path.insert(1, '/raid/scratch/wongj/mywork/3x2pt/3x2pt_pipeline_new/')

import gaussian_cl_likelihood
from gaussian_cl_likelihood.python import cosmosis_utils, simulation, posteriors
from angular_binning import loop_likelihood_nbin, posterior, mask, param_grids, covariance, error_vs_nbin, like_bp_gauss_mix

np.seterr(divide='ignore')

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['mathtext.fontset'] = 'cm'

#nbins = int(os.environ['NBIN'])
#nside = int(os.environ['NSIDE'])

nbins = 3
nside = 256

pi = np.pi

bins = np.arange(1,nbins+1,1)

n_bandpowers = 8
lmin = 50
lmax = 500
bandpower_spacing = 'log'

pbl = gaussian_cl_likelihood.python.simulation.get_binning_matrix(
    n_bandpowers=n_bandpowers,
    output_lmin=lmin,
    output_lmax=lmax,
    bp_spacing=bandpower_spacing)

#print(pbl)

sz = 1.0/(nbins+2)

fig = matplotlib.pyplot.figure(figsize=(10,10))

def open_dat(fname):
    dat_arr = []
    with open(fname) as f:
        for line in f:
            column = line.split()
            if not line.startswith('#'):
                dat_i = float(column[0])
                dat_arr.append(dat_i)
    dat_arr = np.asarray(dat_arr)
    return dat_arr

save_dir = '/raid/scratch/wongj/mywork/3x2pt/2209_256_DZ01_MEASUREMENT/'

colors = ['#0077BB','#33BBEE','#009988','#EE7733','#CC3311','#EE3377','#BBBBBB']

for i in bins:
    for j in bins:
        if i>=j:

            ell = open_dat(save_dir + 'theory_cls/shear_cl/ell_measured.txt')
            bp = open_dat(save_dir + 'theory_cls/shear_cl/PCl_Bandpowers_EE_bin_{}_{}.txt'.format(i,j))

            ell_measured = open_dat(save_dir + 'measured_3x2pt_bps/shear_bp/Cl_EE/ell_measured.txt')
            bp_measured = open_dat(save_dir + 'measured_3x2pt_bps/shear_bp/Cl_EE/bin_{}_{}.txt'.format(i, j))

            rect = ((i * sz) + (0.08 * i) - 0.15, (j * sz) + (0.04 * j) - 0.0875, 1.25*sz, sz)
            #rect = (i*sz,j*sz,sz,sz)
            ax =fig.add_axes(rect)

            plt.plot(ell, bp, label='Theoretical Spectra', color='black',zorder=1)
            plt.plot(ell,bp_measured,color=colors[3],label='Measured Spectra',linestyle='--',zorder=10)
            plt.xscale('log')

            
            if i == 1 and j == 1:
                plt.xlabel("$\\ell$",fontsize=15)
                labelstr = str("$\\ell (\\ell+1) C_\\ell^{\delta_{g}\delta_{g}} / 2 \\pi$")
                plt.ylabel(labelstr,fontsize=15)
                ax.legend(bbox_to_anchor=(0.9, 1.6),fontsize=13.5)

            if j == 1:
                plt.xlabel("$\\ell$",fontsize=15)

            if i == j:
                labelstr = str("$\\ell (\\ell+1) C_\\ell^{\delta_{g}\delta_{g}} / 2 \\pi$")
                plt.ylabel(labelstr,fontsize=15)
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            '''
            if j!=1:
                plt.gca().xaxis.set_ticklabels([])

            if i!=j:
                plt.gca().yaxis.set_ticklabels([])

            if j ==1:
                ax.set_ylim([-2e-5,4e-4])

            if j == 2:
                ax.set_ylim([-2e-5,2e-4])

            if j == 3:
                ax.set_ylim([-2e-5,1.75e-4])
                #ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='y')
            '''           

            ax.minorticks_on()

            ax.tick_params(which='both',axis='both',right=True, top=True, labelright=False, labeltop=False,left=True,bottom=True,labelleft=True,labelbottom=True,direction='in')
            ax.tick_params(length=2.5, which='minor')
            ax.tick_params(length=5.5, which='major')
            ax.tick_params(labelsize=12.5)

            plt.text(0.125,0.75, "("r'$z_{%d}$' ", "r'$z_{%d}$'")" % (i, j), fontsize=15, color='black',transform=ax.transAxes)
            #plt.title('NGal = 5e6, 10 Realisations')

plt.show()

fig = matplotlib.pyplot.figure(figsize=(10,10))

for i in bins:
    for j in bins:
        if i>=j:

            ell = open_dat(save_dir + 'theory_cls/shear_cl/ell_measured.txt')
            bp = open_dat(save_dir + 'theory_cls/shear_cl/PCl_Bandpowers_EE_bin_{}_{}.txt'.format(i,j))

            ell_measured = open_dat(save_dir + 'measured_3x2pt_bps/shear_bp/Cl_EE/ell_measured.txt')
            bp_measured = open_dat(save_dir + 'measured_3x2pt_bps/shear_bp/Cl_EE/bin_{}_{}.txt'.format(i, j))

            rect = ((i * sz) + (0.08 * i) - 0.15, (j * sz) + (0.04 * j) - 0.0875, 1.25*sz, sz)
            #rect = (i*sz,j*sz,sz,sz)
            ax =fig.add_axes(rect)

            plt.plot(ell, ((bp_measured)/(bp))-1, label='Fractional Difference', color=colors[3],zorder=1,linestyle='--')
            plt.axhline(y=0,color='black')
            plt.xscale('log')

            
            if i == 1 and j == 1:
                plt.xlabel("$\\ell$",fontsize=15)
                labelstr = str("$\\ell (\\ell+1) C_\\ell^{\delta_{g}\delta_{g}} / 2 \\pi$")
                plt.ylabel(labelstr,fontsize=15)
                ax.legend(bbox_to_anchor=(0.85, 1.6),fontsize=13.5)

            if j == 1:
                plt.xlabel("$\\ell$",fontsize=15)

            if i == j:
                labelstr = str("$\\ell (\\ell+1) C_\\ell^{\delta_{g}\delta_{g}} / 2 \\pi$")
                plt.ylabel(labelstr,fontsize=15)
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

            if j!=1:
                plt.gca().xaxis.set_ticklabels([])

            if i!=j:
                plt.gca().yaxis.set_ticklabels([])
            ax.set_ylim([-0.2,0.2])

            '''
            if j ==1:
                ax.set_ylim([-2e-5,4e-4])

            if j == 2:
                ax.set_ylim([-2e-5,2e-4])

            if j == 3:
                ax.set_ylim([-2e-5,1.75e-4])
                #ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='y')
            '''           

            ax.minorticks_on()

            ax.tick_params(which='both',axis='both',right=True, top=True, labelright=False, labeltop=False,left=True,bottom=True,labelleft=True,labelbottom=True,direction='in')
            ax.tick_params(length=2.5, which='minor')
            ax.tick_params(length=5.5, which='major')
            ax.tick_params(labelsize=12.5)

            plt.text(0.125,0.75, "("r'$z_{%d}$' ", "r'$z_{%d}$'")" % (i, j), fontsize=15, color='black',transform=ax.transAxes)
            #plt.title('NGal = 5e6, 10 Realisations')

plt.show()


