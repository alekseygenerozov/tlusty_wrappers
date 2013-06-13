from scipy.interpolate import griddata
import numpy as np
import re
import matplotlib.pyplot as plt
import argparse



#Defining physical constants
c=3*10**10
h=6.67*10**-27
kb=1.38*10**-16


# Function to extract frequency from wavelength in angstroms
def get_freq(w):
    return c*10**8/w


# Function to parse file name in order to get list of parameters
def parse_file(file):
    #Extracting parameters from file names
    t=re.search("t[-+]?\d+",file)
    t=re.search("[-+]?\d+",t.group(0))
    
    m=re.search("m[-+]?\d+",file)
    m=re.search("[-+]?\d+",m.group(0))
    
    q=re.search("q[-+]?\d+",file)
    q=re.search("[-+]?\d+",q.group(0))
    
    params=np.array([float(t.group(0)), float(m.group(0)), float(q.group(0))])
    params=params/10.
    
    return params


# Extract spectrum from unit 14 type output file, if mu  is set to -1 then get the flux
def get_spec(file, nfreq=300, mu=-1, nmu=10):
    #Read file assumed to be in the format of tlusty uni 14 file. File has lines with wavelengths and h followed by blocks with 
    #intensities and polariztions
    spec=np.genfromtxt(file, skip_header=1, invalid_raise=False)
    wlh=np.genfromtxt(file, invalid_raise=False)
    

    spec=spec.flatten()
    wlh=wlh.flatten()
    
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
    return (2./fcol**4)*(h*(freq**3)/c/c)/(np.exp(h*freq/kb/fcol/10**Tb)-1)
   

    
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
def construct_table(models):
    models=np.genfromtxt(models, dtype='string')
    params=map(parse_file,models)
    params=np.array(params)
    
    spec=map(get_spec, models)
    spec=np.array(spec)
    spec=map(regrid, spec)

    spec_inv=map(bright_spec, spec)
    spec_inv=np.array(spec_inv)
    
    return (params, spec_inv)

    
# Get parameters corresponding to our disk
def get_params(file):
    return np.genfromtxt(file)


# Use table construct spectra from list of parameters
def params_to_spec(params, table):
    params2=params[:, 2:]
    grid2=griddata(table[0], table[1], params2)#, fill_value=1.e-36)
    grid2=map(bright_inv,grid2)
    
    return np.array(grid2)

# Computes composite spectrum given array of radii, spectra. Also computes corresponding blackbody spectrum using list of Teff that are passed to the function
def sum_spec(r, specs, Teff):   
    dr=r[:, 0]
    r=r[:, 1]
    nu=specs[0,0]
    f=specs[:, 1]
        
    L=np.zeros(len(nu))
    bb=np.zeros(len(nu))
    for i in range(len(r)):
        for j in range(len(nu)):
            L[j]+=2*np.pi*r[i]*dr[i]*f[i,j]
            bb[j]+=2*np.pi*r[i]*dr[i]*(np.pi*Tb_inv(nu[j], Teff[i], fcol=1))
            #bb[j]*=2*np.pi*r[i]*dr[i]
        
    return (np.vstack([nu, bb]), np.vstack([nu, L])) 
    #return np.vstack([nu, (nu*2*np.pi*r*dr*f).sum(axis=0)])


# Calculates a composite disk spectrum given an file containing input radial parameters.
def disk_spec(file):
    #Construct table
    table=construct_table('tmpd')
    disk_params=get_params(file)
    specs=params_to_spec(disk_params, table)

    r=disk_params[:, 0:2]
    Teff=disk_params[:, 2]

    #Finding the total flux
    totf=sum_spec(r, specs, Teff)
    totf1=totf[1]
    totf2=totf[0]

    #max1=np.max(totf1[1])
    #max2=np.max(totf2[1])
    #peak=np.max(max1, max2)

    plt.loglog()
    #plt.figsize(20, 8)

    plt.xlabel("nu [hz]")
    plt.ylabel("nu L_nu [ergs/s]")

    plt.axis([10.**14, 2*10.**18, 10.**38, 10.**44]) 

    plt.plot(totf1[0], totf1[0]*totf1[1])
    plt.plot(totf2[0], totf2[0]*totf2[1])
    plt.show()

def main():
    parser=argparse.ArgumentParser(
        description='Computes composite disk spectrum based on TLUSTY disk spectra')
    parser.add_argument('f',
        help='file name containing radial disk profile ',
        type=str)

    args=parser.parse_args()
    f=args.f

    disk_spec(f)

if __name__ == '__main__':
    main()
    
    
# def Knu(es, f):
#     return -0.3333333333333333 + (0.41997368329829105*(-1. + 2.*es - \
#     1.*es**2))/(2. - 12.*es + 30.*es**2 - 40.*es**3 + 30.*es**4 - \
#     12.*es**5 + 2.*es**6 - 27.*es**2*f + 108.*es**3*f - 162.*es**4*f + \
#     108.*es**5*f - 27.*es**6*f + np.sqrt(-4.*(-1. + 2.*es - 1.*es**2)**6 + \
#     (2. - 12.*es + 30.*es**2 - 40.*es**3 + 30.*es**4 - 12.*es**5 + \
#     2.*es**6 - 27.*es**2*f + 108.*es**3*f - 162.*es**4*f + 108.*es**5*f - \
#     27.*es**6*f)**2))**0.3333333333333333 + (0.26456684199469993*(2. - \
#     12.*es + 30.*es**2 - 40.*es**3 + 30.*es**4 - 12.*es**5 + 2.*es**6 - \
#     27.*es**2*f + 108.*es**3*f - 162.*es**4*f + 108.*es**5*f - \
#     27.*es**6*f + np.sqrt(-4.*(-1. + 2.*es - 1.*es**2)**6 + (2. - 12.*es + \
#     30.*es**2 - 40.*es**3 + 30.*es**4 - 12.*es**5 + 2.*es**6 - \
#     27.*es**2*f + 108.*es**3*f - 162.*es**4*f + 108.*es**5*f - \
#     27.*es**6*f)**2))**0.3333333333333333)/(-1. + 2.*es - 1.*es**2)
    

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


