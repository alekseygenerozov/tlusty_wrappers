# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

from scipy.interpolate import griddata
import numpy as np
import re

# <codecell>

#Defining physical constants
c=3*10**10
h=6.67*10**-27
kb=1.38*10**-16

# Function to extract frequency from wavelength in angstroms
def get_freq(w):
    return c*10**8/w

# Function to parse file name in order to get list of parameters
def parse_file(file):
    t=re.search("t[\-0-9]*",file)
    t=re.search("[-+]?\d+",t.group(0))
    
    m=re.search("m[\-0-9]*",file)
    m=re.search("[-+]?\d+",m.group(0))
    
    q=re.search("q[\-0-9]*",file)
    q=re.search("[-+]?\d+",q.group(0))
    
    params=np.array([float(t.group(0)), float(m.group(0)), float(q.group(0))])
    params=params/10.
    
    return params


# Extract spectrum from unit 14 type output file, if mu  is set to -1 then get the flux
def get_spec(file, nfreq=300, mu=-1, nmu=10):
    spec=np.genfromtxt(file, skip_header=1, invalid_raise=False)
    wlh=np.genfromtxt(file, invalid_raise=False)
    
    spec=spec.flatten()
    wlh=wlh.flatten()
    
    spec=spec[::2]
    wl=wlh[::2]
    
    if mu == -1:
        h=wlh[1::2]
        spec=4*pi*h
        spec=np.vstack((wl, spec))
        
        return spec
        
    
    spec=spec.reshape(( nfreq, nmu))   
    spec=spec[:, mu]
    
    spec=np.vstack((wl,spec)) 
    return spec


# Regrid spectrum in wavelength space. Want all spectra to be on the same grid in wavelength space when we interpolate
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
    grid2=griddata(table[0], table[1], params2)
    grid2=map(bright_inv,grid2)
    
    return np.array(grid2)


# Takes composite interpolated 
def sum_spec(r, specs):
    dr=r[:, 0]
    r=r[:, 1]
    
    nu=specs[0,0]
    f=specs[:, 1]
    
    return np.vstack([nu, (2*pi*r*dr*f).sum(axis=0)])
    
    
    

# <codecell>

disk_params=get_params('disk_params')

t=construct_table('tmpd')
specs=params_to_spec(disk_params, t)

# <codecell>

tmp=grid2.reshape((len(disk_params)*300*2, 1))

for i in range(len(tmp)):
    if (np.isnan(tmp[i])):
        print "test"
        

# <codecell>

#x=zeros(30)
#x=x.reshape(10, 3)

#for i in range(0, 10):
 #   x[i]=np.array([5, 3+0.1*i, -1])

#grid1=griddata(params, spec, [[5.05, 4, -1]])
#grid2=griddata(params, spec_inv, [[5.05, 4, -1]])

#grid1=griddata(params, spec, [[5.05, 4, -1]])
#grid2=griddata(params, spec_inv, [[5.05, 4, -1]])


#
#grid2=np.array(grid2)
#spec2[1, 150],grid1[:, 1, 150]
#grid2=map(bright_inv,grid2)
#grid2=np.array(grid2)
#spec2[1, 150],grid1[:, 1, 150]

# <codecell>

spec1=get_spec('t40m40q-100.14')
spec1=regrid(spec1)
spec2=get_spec('t41m40q-100.14')
spec2=regrid(spec2)
spec3=get_spec('test.14')
spec3=regrid(spec3)

loglog()
figsize(20, 8)
axis([10**14, 10**18, 10**10, 10**15])


#plot(spec2[0], spec2[0]*spec2[1])
#plot(spec1[0], spec1[0]*spec1[1])
#plot(spec1[0], spec1[0]*spec1[1])
plot(specs[0, 0], specs[0,0]*specs[0,1])
#plot(spec3[0], spec3[0]*spec3[1])

#spec1[ 150],grid1[:, ])
#spec1[ 150],grid1[:, 150]
#savefig('interpolated1.pdf')

# <codecell>

totf=sum_spec(disk_params[0:1], specs)

loglog()
figsize(20, 8)

xlabel("Freq [hz]")
ylabel("Flux [ergs/s]")

axis([10.**14, 2*10.**16, 10.**35, 10.**42]) 

plot(totf[0], totf[0]*totf[1])

