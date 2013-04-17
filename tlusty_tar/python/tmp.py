from astropy.io  import ascii

dat=ascii.read('params.in')



if dat[0][3]:
    print 'hello world'