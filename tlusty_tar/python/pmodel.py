#!/usr/bin/env python
import argparse
import numpy as np
import re
import tlusty_runs as tr
import matplotlib.pyplot as plt
import shlex

from matplotlib.backends.backend_pdf import PdfPages



def sanitize(f, end='/fort.7'):
    to_replace=re.compile('fort.*|/*')
    f=re.sub(to_replace,'',f)
    if f:
        return f+end
    else:
        return 'fort.7'


def pmodel(f, col=0, z=False):
	f=sanitize(f)

	to_read=open(f, 'r')
	dat=to_read.read()
	dat=np.array(shlex.split(dat), dtype=float)
	nd=dat[0]
	blocks=dat[1]

	m=dat[2:nd+2]
	atm=dat[nd+2:]
	atm.shape=[nd,blocks]

	plt.loglog()
	if z:
		z=atm[:,-1]
		plt.plot(z,atm[:,col])
	else:
		plt.plot(m,atm[:,col])




	
def main():
	parser=argparse.ArgumentParser(
        description='Wrapper for plotting TLUSTY output')
	parser.add_argument('-f', '--files',
		help='name of dir containing tlusty output to plot',
		nargs='*')
	parser.add_argument('-c', '--col',
		help='column number to plot',
		default=1,
		type=int)
	parser.add_argument('-z', '--z',
		help='flag to switch ind var from m to z in plot',
		default=False,
		action=store_true
		)

	parser.parse_args()
	args=parser.parse_args()
	files=args.files
	col=args.col
	z=args.z


	if files:
		for f in files:
			pmodel(f,col=col, z=z)
		plt.show()
	else:
		parser.print_help()




if __name__ == '__main__':
    main()
