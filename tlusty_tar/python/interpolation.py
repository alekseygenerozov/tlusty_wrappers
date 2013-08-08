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



##Defining physical constants
G=6.67*10**-8
c=3*10**10
h=6.67*10**-27
kb=1.38*10**-16
M_sun=2.*10**33


##Run a command from the bash shell
def bash_command(cmd):
     process=subprocess.Popen(['/bin/bash', '-c', cmd],  stdout=subprocess.PIPE)
     process.wait()
     return process

##Function to extract frequency from wavelength in angstroms
def get_freq(w):
    return c*10**8/w


# Function to parse file name in order to get list of parameters.
def parse_file(f, string=False):
    f2=re.sub('.14\Z', '', f)
    refloat =r'[+-]?\d+\.?\d*'
    t=re.search("t"+refloat,f2)
    t=re.search(refloat,t.group(0))
    
    m=re.search("m"+refloat,f2)
    m=re.search(refloat,m.group(0))
    
    q=re.search("q"+refloat,f2)
    q=re.search(refloat,q.group(0))
    
    if string:
        params=(t.group(0), m.group(0), q.group(0))
    else:
        params=(float(t.group(0))/10., float(m.group(0))/10., float(q.group(0))/10.)
    
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
    #Second moment of radiative intensity
    h=wlh[1::2]
    flux=4*np.pi*h
    flux.shape=(300,1)

    #Add flux to end of intensity arrays for every frequency
    spec.shape=(nfreq, nmu)
    spec=np.append(spec, flux, axis=1)

    return (wl,spec) 



# Regrid spectrum in wavelength space. Want all spectra to be on the same grid in wavelength space when we interpolate; In frequency the default limits for the 
# the wavelength correspond to 2.4e14 and 5e19 Hz.
def regrid(spec, wlo=0.06, whi=12400, nws=300):
    wgrid=np.log(wlo)+np.log(whi/wlo)*np.arange(0, nws)/(nws-1)
    newspec=griddata(np.log(spec[0]), spec[1], wgrid, fill_value=1.e-36)
    newspec=newspec[::-1]

    wgrid=np.exp(wgrid)
    freq=map(get_freq, wgrid)
    freq=freq[::-1]
    
    freq2=np.empty_like(newspec)
    for i in range(len(freq2)):
        freq2[i].fill(freq[i])
    spec=np.vstack([freq2, newspec])
    spec.shape=(2,300,11)

    return spec



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

    
# # Get brightness temperature from a spectrum, with default color correction factor of 2 as described in ApJS 164 530D
# # returns the log of the brightness temperature
# def bright_spec(spec, fcol=2):
#     freq=spec[0]
#     intens=spec[1]
#     farr=np.empty(len(freq))
#     farr.fill(fcol)
    
#     bright=map(Tb, freq, intens, farr)
#     return np.vstack((freq, bright))


# # Invert brightness temperature to obtain the spectrum
# def bright_inv(spec_inv, fcol=2):
#     freq=spec_inv[0]
#     temp=spec_inv[1]
    
#     farr=np.empty(len(freq))
#     farr.fill(fcol)
    
#     intens=map(Tb_inv, freq, temp, farr)
#     spec=np.vstack((freq, intens))
#     return spec


##Constructing table of spectra from unit 14 files
def construct_table(models, logi=False):
    models=np.genfromtxt(models, dtype='string')
    params=map(parse_file,models)
    params=np.array(params)
    #params=params[:, ::2]
    
    spec=map(get_spec, models)
    spec=map(regrid, spec)
    spec=np.array(spec)

    if logi:
        spec[:, 1]=np.log10(spec[:, 1])
        spec=np.array(spec)
    else:
        spec[:, 1]=map(Tb, spec[:,0],spec[:,1])
        spec=np.array(spec)

    return (params, spec)

    
##Get parameters corresponding to our disk
def get_params(file):
    return np.genfromtxt(file)


##Use table construct spectra from list of parameters
def params_to_spec(params, table, method='', logi=False, mu=-1):
    intens=table[1]
    #Interpolating based on table
    if method:
        print method
        grid2=griddata(table[0], intens[:,:,:,mu], params, method=method)
    else:
        grid2=griddata(table[0], intens[:,:,:,mu], params)

    good=np.empty(len(grid2), dtype=bool)
    for i in range(len(grid2)):
        #print np.any(np.isnan(grid2[i]))
        if np.any(np.isnan(grid2[i])):
            continue
        elif logi:
            print logi
            grid2[i,1]=10.**grid2[i,1]
        else:
            grid2[i,1]=np.array(map(Tb_inv,grid2[i,0],grid2[i,1]))

    
    #Return the interpolated spectra for parameters within our grid
    return grid2


##Computes composite spectrum given array of radii, spectra. Also computes corresponding blackbody spectrum using list of Teff that are passed to the function
def sum_spec(r, specs, Teff, Qg, M=10.**9): 
    r=r*(G*M*M_sun/c**2)  
    lr=np.log10(r)
    #For now assuming a logarithmically evenly spaced grid--for simplicity
    dlr=np.diff(lr)[0]
    print dlr

    dr={}
    for i in range(len(lr)):
        dr[i]=10.**(lr[i]+(dlr/2))-10.**(lr[i]-(dlr/2))
    #Implicitly assuming that the first entry is not outside our grid -- not ideal!
    nu=specs[0,0]
    f=specs[:, 1]

    L=np.zeros(len(nu))
    L2=np.zeros(len(nu))
    bb=np.zeros(len(nu))
    gb=np.zeros(len(nu))
    for i in range(len(r)):
        #Check for any invalid entries -- these correspond to the parameters outside the edge of our table
        valid=not np.any(np.isnan(specs[i]))
        rad=2*np.pi*r[i]*dr[i]
        #For every frequncy under consideration
        for j in range(len(nu)):
            #If we have a valid interpolated spectrum...
            if valid:
                bb[j]+=rad*(np.pi*Tb_inv(nu[j], Teff[i], fcol=1))
                gb[j]+=rad*(gray.gb(nu[j], 10.**Teff[i], 10.**Qg[i]))
                L[j] +=rad*f[i,j]
                L2[j]+=rad*f[i,j]
            #Otherwise...
            else:
                bb[j]+=rad*(np.pi*Tb_inv(nu[j], Teff[i], fcol=1))
                gb[j]+=rad*(gray.gb(nu[j], 10.**Teff[i], 10.**Qg[i]))
                L[j] +=rad*(np.pi*Tb_inv(nu[j], Teff[i], fcol=1))

    return (nu, gb, bb, L, L2) 

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

    if np.any(np.isnan(testspec_interp)):
        print "Warning -- unable to interpolate spectrum for specified parameters."
        return [testspec[0], testspec[0]*testspec[1], testspec_interp[0]*testspec_interp[1]]
    return [testspec[0], testspec[0]*testspec[1], testspec_interp[0]*testspec_interp[1]]


#takes list of models creates animation of computed spectrum plus interpolated spectrum
def animate_test_spec(models, table=[], tablef='tmpd', method='', logi=False):
    def update_img(n):
        spec_plot.set_ydata(spec[n])
        spec_plot_interp.set_ydata(spec_interp[n])
        label.set_text('log(teff)= '+str(teffs[n]))
        return [spec_plot, spec_plot_interp]

    fig, ax=plt.subplots()
    plt.loglog()
    plt.xlabel(r"$\nu$ [hz]")
    plt.ylabel(r"$\nu F_{\nu}$ [ergs s$^{-1}$ cm$^{-2}$ ]")
    plt.axis([10.**14, 2*10.**18, 10.**6, 10.**16])
    label=ax.text(0.02, 0.95, '', transform=ax.transAxes)

    (spec, spec_interp, teffs)=([],[],[])

    nu=test_spec(models[0], table=table, tablef=tablef, method=method, logi=logi)[0]
    tmpspec=(test_spec(models[0], table=table, tablef=tablef, method=method, logi=logi)[1])
    tmpspec_interp=(test_spec(models[0], table=table, tablef=tablef, method=method, logi=logi)[1])
    spec_plot,=ax.plot(nu, tmpspec)
    spec_plot_interp,=ax.plot(nu, tmpspec_interp)


    (ts,ms,qs)=parse_file(models[0], True)
    writer = animation.writers['ffmpeg'](fps=10) 
    for m in models:
        m=m.rstrip()
        (ts2,ms2,qs2)=parse_file(m, True)
        teffs.append(float(ts2)/10)
        if not(ms2==ms) or not(qs2==qs):
            print len(spec)
            ani = animation.FuncAnimation(fig,update_img,len(spec),interval=50)
            ani.save('interp_m'+ms+'_q'+qs+'.mp4',writer=writer,dpi=100)
            (spec, spec_interp, teffs)=([],[],[])
            (ms,qs)=(ms2,qs2)

        #spec=test_spec(m, table=table, tablef=tablef, method=method, logi=logi) 
        spec.append(test_spec(m, table=table, tablef=tablef, method=method, logi=logi)[1])
        spec_interp.append(test_spec(m, table=table, tablef=tablef, method=method, logi=logi)[2])


##Calculates a composite disk spectrum given an file containing input radial parameters.
def disk_spec(f, table=[], tablef='tmpd', method='', logi=False, mu=-1):
    #Construct table if necessary
    if table==[]:
        table=construct_table(tablef, logi=logi)
    disk_params=np.genfromtxt(f)#, skip_header=1)
    specs=params_to_spec(disk_params[:, 1:4], table, method=method, logi=logi, mu=mu)

    r=disk_params[:, 0]
    Teff=disk_params[:, 1]
    Qg=disk_params[:, 3]

    #Finding the total flux, for the given parameters
    totf=sum_spec(r, specs, Teff, Qg)
    nu=totf[0]
    totfg=totf[1]
    totfb=totf[2]
    totft=totf[3]
    totft2=totf[4]

    fig,ax=plt.subplots(nrows=2, ncols=1, figsize=(6,8),sharex=True, subplot_kw=dict(adjustable='datalim'))
    #plt.title(str(bin_params[0])+" "+str(bin_params[1])+" "+str(bin_params[2]))
    plt.xlabel(r"$\nu$ [hz]")

    #Plotting the composite disk spectrum
    ax[0].set_ylabel(r"$\nu L_{\nu}$ [ergs s$^{-1}$]")
    ax[0].set_xlim(10.**14, 10.**17)
    ax[0].set_ylim(10.**40, 10.**47)
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    
    ax[0].plot(nu, nu*totfg)
    ax[0].plot(nu, nu*totfb)
    ax[0].plot(nu, nu*totft)
    ax[0].plot(nu, nu*totft2)

    #Plotting the contributions of individual annuli
    ax[1].set_ylim(10.**-5, 1.)
    ax[1].set_xlim(10.**14, 10.**17)
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    for i in range(len(specs)):
        ax[1].plot(nu, specs[i, 1])

    return fig

def main():
    parser=argparse.ArgumentParser(
        description='Either takes list of disk parameters and computes composite disk spectrum from table or compares spectra interpolated'+
          ' from table to those in a test directory')
    parser.add_argument('-d', '--disk',
        help='file with list of files containing radial disk profiles ',
        default='')
    parser.add_argument('-mu', '--mu',
        help='cosine of inslincation angle for the case of a disk ',
        type=int,
        default=-1)
    parser.add_argument('-t', '--test',
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
    # parser.add_argument('-a', '--animate',
    #     help='For the case of test spectra, specifies that a movie should be made rather than a static pdf.'
    #     a)

    args=parser.parse_args()
    t=args.test
    d=args.disk
    method=args.method
    tablef=args.tablefile
    logi=args.logi
    skip=args.skip
    mu=args.mu

    print logi
  
    if d:
        table=construct_table(tablef, logi=logi)
        pdf_pages = PdfPages('composite.pdf')
        param_files=np.genfromtxt(d, dtype=str)
        for pf in param_files:
            #print param_files
            fig=disk_spec(pf, table=table, tablef=tablef, method=method,logi=logi, mu=mu)
            pdf_pages.savefig(fig)
        pdf_pages.close()
    elif t:
        table=construct_table(tablef, logi=logi)
        # pdf_pages = PdfPages('interp_test.pdf')
        process=bash_command('echo '+t+'/*14')

        # bash_command('rm interp_log')
        # logfile=open('interp_log', 'a')

        # deviation_list=np.empty(0)
        models=process.stdout.readlines()[0]
        models=shlex.split(models)
        if skip:
            skip_models=np.genfromtxt(t+'/'+skip, dtype=str)
            models=np.setdiff1d(models, skip_models)
        models=np.array(models, dtype=str)
        #Sort models used by teff, then m, and finally q
        params=map(parse_file, models)
        dtype=[('teff',float), ('m', float), ('q', float)]
        params=np.array(params, dtype=dtype)
        order=np.argsort(params, order=['q', 'm', 'teff'])
        models=models[order]

        #Create animations comparing the interpolated spectra to tlusty spectra found in the specified directory
        animate_test_spec(models, table=table, tablef=tablef, method=method, logi=logi)

        
        # for o in order:
        #     params2=params[o]
        #     print params2
        #     diff=np.diff(params2[:, 0])
        #     br=np.where((diff>10^-6)

        #     print br

            # print np.split(params2, br)
        
        #animate_test_spec(models, table=table, tablef=tablef, method=method, logi=logi)


        # params=np.array(map(tr.parse_file, models))
        # params=map(tr.parse_file, models)
        # params=np.array(params)
        # teffs=params[:,0]
        # order=np.argsort(teffs)

        # models[order]
        # for m in models:
        #     m=m.rstrip()
        #     spec=test_spec(m, table=table, tablef=tablef, method=method, logi=logi)
        #     fig=spec[0]
        #     deviation=spec[1]
        #     logfile.write(m+" "+str(deviation)+"\n")
        #     if deviation!=-1:
        #         deviation_list=np.append(deviation_list,deviation)
        #     #Save file 
        #     pdf_pages.savefig(fig)

        # print deviation_list
        # pdf_pages.close()
        # fig2=plt.figure()
        # plt.hist(deviation_list, bins=50, normed=1)
        # fig2.savefig('dev_hist.pdf')
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


