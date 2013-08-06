#!/usr/bin/env python

from scipy.interpolate import griddata
import numpy as np
import numpy.ma as ma
import re

from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.animation as animation

import argparse
import subprocess

import graybody as gray
import tlusty_runs as tr


import warnings

import shlex

from types import *



##Defining physical constants
c=3*10**10
h=6.67*10**-27
kb=1.38*10**-16


##Run a command from the bash shell
def bash_command(cmd):
     process=subprocess.Popen(['/bin/bash', '-c', cmd],  stdout=subprocess.PIPE)
     process.wait()
     return process

##Function to extract frequency from wavelength in angstroms
def get_freq(w):
    return c*10**8/w


# Function to parse file name in order to get list of parameters
def parse_file(f):
    f2=re.sub('.14\Z', '', f)
    refloat =r'[+-]?\d+\.?\d*'
    t=re.search("t"+refloat,f2)
    t=re.search(refloat,t.group(0))
    
    m=re.search("m"+refloat,f2)
    m=re.search(refloat,m.group(0))
    
    q=re.search("q"+refloat,f2)
    q=re.search(refloat,q.group(0))
    
    
    params=np.array([float(t.group(0)), float(m.group(0)), float(q.group(0))])
    params=params/10.
    return params

# Extract spectrum from unit 14 type output file, if mu  is set to -1 then get the flux
def get_spec(file, nfreq=300, mu=-1, nmu=10):
    #Read file assumed to be in the format of tlusty uni 14 file. File has lines with wavelengths and h followed by blocks with 
    #intensities and polariztions
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec=np.genfromtxt(file, skip_header=1, invalid_raise=False)
        wlh=np.genfromtxt(file, invalid_raise=False)
    
    spec=spec.flatten()
    wlh=wlh.flatten()
    #Skip over fluxes and polarizations
    spec=spec[::2]
    wl=wlh[::2]

    #mu is and index for the polar angle to use. -1 corresponds to returning emergent flux from the disk
    if mu == -1:
        h=wlh[1::2]
        spec=4*np.pi*h
        spec=np.vstack((wl, spec))        
        return spec
        
    spec=spec.reshape(( nfreq, nmu))   
    spec=spec[:, mu]

    print np.shape(wl), np.shape(spec), f
    
    spec=np.vstack((wl,spec)) 
    return spec


# Regrid spectrum in wavelength space. Want all spectra to be on the same grid in wavelength space when we interpolate; In frequency the default limits for the 
# the wavelength correspond to 2.4e14 and 5e19 Hz.
def regrid(spec, wlo=0.06, whi=12400, nws=300):
    wgrid=np.log(wlo)+np.log(whi/wlo)*np.arange(0, nws)/(nws-1)
    spec[0]=np.log(spec[0])
    
    newspec=griddata(spec[0], spec[1], wgrid, fill_value=1.e-36)
    
    wgrid=np.exp(wgrid)
    freq=map(get_freq, wgrid)
    
    freq=freq[::-1]
    newspec=newspec[::-1]
    return np.vstack((freq, newspec))



# Function calculates brightness temperature for a given frequncy and intensity
def Tb(freq, intens, fcol=2):
    return np.log10(h*freq/kb/fcol/np.log(1+(2*h*freq**3/fcol**4/intens/c/c)))


# Invert brightness temperature to obtain the intensity
def Tb_inv(freq, Tb, fcol=2):
    if np.isinf(np.exp(h*freq/kb/fcol/10**Tb)):
        return 0
    try:
        return (2./fcol**4)*(h*(freq**3)/c/c)/(np.exp(h*freq/kb/fcol/10**Tb)-1)
    except:
        return 0

    
# Get brightness temperature from a spectrum, with default color correction factor of 2 as described in ApJS 164 530D
# returns the log of the brightness temperature
def bright_spec(spec, fcol=2):
    freq=spec[0]
    intens=spec[1]
    farr=np.empty(len(freq))
    farr.fill(fcol)
    
    bright=map(Tb, freq, intens, farr)
    return np.vstack((freq, bright))


# Invert brightness temperature to obtain the spectrum
def bright_inv(spec_inv, fcol=2):
    freq=spec_inv[0]
    temp=spec_inv[1]
    
    farr=np.empty(len(freq))
    farr.fill(fcol)
    
    intens=map(Tb_inv, freq, temp, farr)
    
    spec=np.vstack((freq, intens))
    return spec


# Constructing table of spectra from unit 14 files
def construct_table(models, logi=False):
    models=np.genfromtxt(models, dtype='string')
    params=map(parse_file,models)
    params=np.array(params)
    params=params[:, ::2]
    
    spec=map(get_spec, models)
    spec=np.array(spec)
    spec=map(regrid, spec)
    spec=np.array(spec)

    if logi:
        spec[:, 1]=np.log10(spec[:, 1])
        spec=np.array(spec)
    else:
        spec=map(bright_spec, spec)
        spec=np.array(spec)

    return (params, spec)

    
##Get parameters corresponding to our disk
def get_params(file):
    return np.genfromtxt(file)


##Use table construct spectra from list of parameters
def params_to_spec(params, table, method='', logi=False):
    # method='cubic'
    # if linear:
    #     method='linear'
    if method:
        print method
        grid2=griddata(table[0], table[1], params, method=method)
    else:
        grid2=griddata(table[0], table[1], params)
    #grid2=ma.array(grid2)
    good=np.empty(len(grid2), dtype=bool)
    for i in range(len(grid2)):
        #print np.any(np.isnan(grid2[i]))
        if np.any(np.isnan(grid2[i])):
            continue
            #grid2[i]=ma.masked
        elif logi:
            grid2[i,1]=10.**grid2[i,1]
            #good[i]=True
        else:
            #print params[i]
            grid2[i]=bright_inv(grid2[i])
            #good[i]=True
    
    #Return the interpolated spectra for parameters within our grid
    return grid2
    #return grid2[good]

# Computes composite spectrum given array of radii, spectra. Also computes corresponding blackbody spectrum using list of Teff that are passed to the function
def sum_spec(r, specs, Teff, Qg):   
    dr=r[:, 0]
    r=r[:, 1]
    #print specs.shape
    #Implicitly assuming that the first entry is not outside our grid -- not ideal!
    nu=specs[0,0]
    #print nu
    f=specs[:, 1]

    L=np.zeros(len(nu))
    L2=np.zeros(len(nu))
    bb=np.zeros(len(nu))
    gb=np.zeros(len(nu))
    for i in range(len(r)):
        #Check for any invalid entries -- these correspond to the parameters outside the edge of our table
        valid=not np.any(np.isnan(specs[i]))
        rad=2*np.pi*r[i]*dr[i]
        for j in range(len(nu)):
            #bb[j]+=2*np.pi*r[i]*dr[i]*(np.pi*Tb_inv(nu[j], Teff[i], fcol=1))
            #gb[j]+=2*np.pi*r[i]*dr[i]*(gray.gb(nu[j], 10.**Teff[i], 10.**Qg[i]))
            if valid:
                bb[j]+=rad*(np.pi*Tb_inv(nu[j], Teff[i], fcol=1))
                gb[j]+=rad*(gray.gb(nu[j], 10.**Teff[i], 10.**Qg[i]))
                L[j] +=rad*f[i,j]
                L2[j]+=rad*f[i,j]
            else:
                bb[j]+=rad*(np.pi*Tb_inv(nu[j], Teff[i], fcol=1))
                gb[j]+=rad*(gray.gb(nu[j], 10.**Teff[i], 10.**Qg[i]))
                L[j] +=rad*(np.pi*Tb_inv(nu[j], Teff[i], fcol=1))
            
            #bb[j]*=2*np.pi*r[i]*dr[i]
        
    return (np.vstack([nu, gb]),np.vstack([nu, bb]), np.vstack([nu, L]),np.vstack([nu, L2])) 
    #return np.vstack([nu, (nu*2*np.pi*r*dr*f).sum(axis=0)])


##Compares tlusty spectrum to one that is interpolated from a table 
def test_spec(f, table=[], tablef='tmpd', method='', logi=False):
    if table==[]:
        table=construct_table(tablef, logi=logi)
    #Read tlusty spectrum from file. 
    testspec=get_spec(f)
    testspec=regrid(testspec)
    
    #Construct interpolated spectrum for the params corresponding to file.
    params=parse_file(f)
    params_s='t'+str(10*params[0])+'m'+str(10*params[1])+'q'+str(10*params[2])
    print params
    params=params[::2]
    testspec_interp=params_to_spec([params], table, method=method, logi=logi)
    testspec_interp=testspec_interp[0]

    #Construct figure, comparing the interpolated spectrum to that directly computed from TLUSTY. 
    fig,ax=plt.subplots()
    plt.loglog()
    plt.xlabel(r"$\nu$ [hz]")
    plt.ylabel(r"$\nu F_{\nu}$ [ergs s$^{-1}$ cm$^{-2}$ ]")
    plt.axis([10.**14, 2*10.**18, 10.**6, 10.**16])
    plt.title(r'%s'%params_s)
    #plt.plot(diff[0], diff[1])
    ax.plot(testspec[0], testspec[0]*testspec[1])
    ax.plot(testspec_interp[0], testspec_interp[0]*testspec_interp[1])


    #plt.plot(testspec[0], 1.1*testspec[0]*testspec[1])
    #plt.plot(testspec_interp[0], testspec_interp[0]*testspec_interp[1])
    if np.any(np.isnan(testspec_interp)):
        print "Warning -- unable to interpolate spectrum for specified parameters."
        return [testspec_interp[0], testspec_interp[0]*testspec_interp[1], -1]
    #Finding max fractional deviation of the interpolated spectrum from the tlusty spectrum.
    peak=np.max(testspec[1])
    cut=testspec[1]>0.01*peak

    #Comparing the tlusty and interpolated spectra
    compare=[testspec[0,cut],(testspec_interp[1,cut])/(testspec[1,cut])]
    args=np.argsort(np.abs(compare[1]-1))
    max_deviation=compare[1][args[-1]]
    return [testspec_interp[0], testspec_interp[0]*testspec_interp[1], max_deviation]


# Calculates a composite disk spectrum given an file containing input radial parameters.
def disk_spec(f, table=[], tablef='tmpd', method='', logi=False):
    #Construct table
    if table==[]:
        table=construct_table(tablef, logi=logi)
    bin_params=np.fromfile(f, dtype=float, sep=' ',count=3)
    disk_params=np.genfromtxt(f, skip_header=1)


    specs=params_to_spec(disk_params[:, 2::2], table, method=method, logi=logi)

    r=disk_params[:, 0:2]
    Teff=disk_params[:, 2]
    Qg=disk_params[:, 4]

    #Finding the total flux
    totf=sum_spec(r, specs, Teff, Qg)
    totfg=totf[0]
    totfb=totf[1]
    totft=totf[2]
    totft2=totf[3]
    #totfg2=np.genfromtxt('gray_test')
    

    fig=plt.figure()
    plt.loglog()
    #plt.figsize(20, 8)
    plt.xlabel(r"$\nu$ [hz]")
    plt.ylabel(r"$\nu L_{\nu}$ [ergs s$^{-1}$]")
    plt.axis([10.**14, 2*10.**18, 10.**38, 10.**45]) 
    plt.title(str(bin_params[0])+" "+str(bin_params[1])+" "+str(bin_params[2]))

    plt.plot(totfg[0], totfg[0]*totfg[1])
    plt.plot(totfb[0], totfb[0]*totfb[1])
    plt.plot(totft[0], totft[0]*totft[1])
    plt.plot(totft2[0], totft2[0]*totft2[1])
    #plt.plot(totfg2[:,0], totfg2[:,0]*totfg2[:,1])
    #plt.show()
    return fig

def main():
    parser=argparse.ArgumentParser(
        description='Either takes list of disk parameters and computes composite disk spectrum from table or compares spectra interpolated'+
          ' from table to those in a test directory')
    parser.add_argument('-f', '--file',
        help='file with list of files containing radial disk profiles ',
        default='')
    parser.add_argument('-d', '--dir',
        help='directory containing the location of models to test',
        default='')
    parser.add_argument('-sf', '--skip',
        help='in case d argument is used this can be used to specify file name containing the names of models to be skipped',
        default='')
    parser.add_argument('-tf', '--tablefile',
        help='name of file containing list of table models',
        default='tmpd')
    parser.add_argument('-m', '--method',
        help='Specify order of interpolation should be linear for table',
        choices=['linear', 'cubic', 'nearest'])
    parser.add_argument('-li', '--logi',
        help='Specifies that the interpolation should be done in terms of log intensity instead of brightness temp.',
        action='store_true')
    parser.add_argument('-a', '--animate',
        help='For the case of test spectra, specifies that a movie should be made rather than a static pdf.',
        action='store_true')

    args=parser.parse_args()
    f=args.file
    d=args.dir
    method=args.method
    tablef=args.tablefile
    logi=args.logi
    skip=args.skip
    animate=args.animate
  
    if f:
        table=construct_table(tablef, logi=logi)
        pdf_pages = PdfPages('composite.pdf')
        param_files=np.genfromtxt(f, dtype=str)
        for pf in param_files:
            #print param_files
            fig=disk_spec(pf, table=table, tablef=tablef, method=method)
            pdf_pages.savefig(fig)
        pdf_pages.close()
    elif d:
        table=construct_table(tablef, logi=logi)
        #pdf_pages = PdfPages('interp_test.pdf')
        process=bash_command('echo '+d+'/*14')

        bash_command('rm interp_log')
        logfile=open('interp_log', 'a')

        deviation_list=np.empty(0)
        models=process.stdout.readlines()[0]
        models=shlex.split(models)
        if skip:
            skip_models=np.genfromtxt(d+'/'+skip, dtype=str)
            models=np.setdiff1d(models, skip_models)

        params=np.array(map(tr.parse_file, models))
        teffs=params[:,0]
        order=np.argsort(teffs)
        models=np.array(models, dtype=str)
        models[order]

        fig, ax=plt.subplots()
        plt.loglog()
        plt.xlabel(r"$\nu$ [hz]")
        plt.ylabel(r"$\nu F_{\nu}$ [ergs s$^{-1}$ cm$^{-2}$ ]")
        plt.axis([10.**14, 2*10.**18, 10.**6, 10.**16])

        tmpspec=(test_spec(models[0], table=table, tablef=tablef, method=method, logi=logi)[1])
        nu=test_spec(models[0], table=table, tablef=tablef, method=method, logi=logi)[0]
        spec_plot,=ax.plot(nu, tmpspec)

        spec=[]
        print len(models)
        for m in models:
            m=m.rstrip()
            spec.append(test_spec(m, table=table, tablef=tablef, method=method, logi=logi)[1])
            nu=test_spec(m, table=table, tablef=tablef, method=method, logi=logi)[1]
            #ims.append((plt.plot(x, y)),)
            
            # deviation=spec[-1]
            # logfile.write(m+" "+str(deviation)+"\n")
            # if deviation!=-1:
            #     deviation_list=np.append(deviation_list,deviation) 
        spec=np.array(spec)
        print spec[0]
        print spec.shape
        #pdf_pages.close()
        def update_img(n):
            spec_plot.set_ydata(spec[n])
            #spec_plot.set_xdata(nu)
            return spec_plot


        #Deviation histogram
        # fig2=plt.figure()
        # plt.hist(deviation_list, bins=50, normed=1)
        # fig2.savefig('dev_hist.pdf')

        ani = animation.FuncAnimation(fig,update_img,len(spec),interval=1)
        writer = animation.writers['ffmpeg'](fps=10) 
        #im_ani = animation.ArtistAnimation(fig3, ims, interval=1, repeat_delay=3000,

        #plt.show()
        ani.save('test.mp4',writer=writer,dpi=100)
    else:
        parser.print_help()

    #pdf_pages.close()



if __name__ == '__main__':
    main()
    
    

    

# def Knu(es, f):
#     return  np.sqrt(-4.*(-1. + 2.*es - 1.*es**2)**6 + \
#     (2. - 12.*es + 30.*es**2 - 40.*es**3 + 30.*es**4 - 12.*es**5 + \
#     2.*es**6 - 27.*es**2*f + 108.*es**3*f - 162.*es**4*f + 108.*es**5*f - \
#     27.*es**6*f)**2)
    

# <codecell>

# table=construct_table('tmpd')
# disk_params=get_params('dparams_pert_total')


# specs=params_to_spec(disk_params, table)

# loglog()
# figsize(20, 8)

# r=disk_params[:, 0:2]
# Teff=disk_params[:, 2]


# totf=sum_spec(r, specs, Teff)
# totf1=totf[1]
# totf2=totf[0]

# grayf=np.genfromtxt("graydisk")


# peak=10.**44
# #peak=np.max([max1, max2])
# #totf2=sum_spec(r[2:], specs[2:])
# #totf3=sum_spec(r[5:], specs[5:])

# loglog()
# figsize(20, 8)

# xlabel("nu [hz]")
# ylabel("nu L_nu [ergs/s]")

# axis([10.**14, 2*10.**18, 10**-6*peak, 2*peak]) 

# plot(totf1[0], totf1[0]*totf1[1])
# plot(totf2[0], totf2[0]*totf2[1])
# plot(grayf[:, 0], grayf[:, 0]*grayf[:, 1])
# #plot(totf2[0], totf2[0]*totf2[1])
# #plot(totf3[0], totf3[0]*totf3[1])

# savefig("composite_pert.pdf")


