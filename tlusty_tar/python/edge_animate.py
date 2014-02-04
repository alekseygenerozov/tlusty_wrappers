#!/usr/bin/env python

import interpolation as interp
import tlusty_runs as tr

from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.animation as animation

import shlex
import numpy as np

import sys

#Importing list of files. Idea of this function is to show
def spec_animate(pat, dpi=100, ly_left=3.26E15, ly_right=3.35E15):
	files=interp.bash_command('echo '+pat+'.14')
	files=shlex.split(files)
	print files

	params=map(tr.parse_file, files)
	params=np.array(params)
	teffs=params[:,0]
	order=np.argsort(teffs)

	m='{0:.5g}'.format(params[0,1]*10)
	q='{0:.5g}'.format(params[0,2]*10)


	#Importing and reordering spectra by teff
	spec=map(interp.get_spec, files)
	for i in range(len(spec)):
		spec[i]=interp.regrid(spec[i], keep=True)
	spec=np.array(spec)
	spec=spec[order]
	teffs=teffs[order]

	
	spec=spec[:,:,:,-1]
	print spec[0,0]

	fig, ax=plt.subplots()
	label=ax.text(0.02, 0.95, '', transform=ax.transAxes)
	plt.loglog()
	plt.xlabel(r"$\nu$ [hz]")
	plt.ylabel(r"$\nu F_{\nu}$ [ergs s$^{-1}$ cm$^{-2}$ ]")
	plt.axis([2.*10.**15, 4.*10.**15, 10.**6, 10.**16])


	spec_plot,=ax.plot(spec[0,0], spec[0,0]*spec[0,1], 'bs-')
	spec_left=spec[0, 1, spec[0,0]<ly_left]
	spec_right=spec[0, 1, spec[0,0]>ly_right]
	edge_pts,=ax.plot([ly_left], [spec_left[-1]*ly_left], 'rs')
	edge_pts2,=ax.plot([ly_right], [spec_right[0]*ly_right], 'rs')
	ratio=spec_left[-1]/spec_right[0]
	label.set_text(str(teffs[0])+" "+str(ratio)+" "+str(spec_left[-1]))

	def update_img(n):
		#print spec[n,0]<ly_left
		spec_left=spec[n, 1, spec[n,0]<ly_left]
		spec_right=spec[n, 1, spec[n,0]>ly_right]
		ratio=spec_left[-1]/spec_right[0]

		spec_plot.set_ydata(spec[n,0]*spec[n,1])
		edge_pts.set_ydata(ly_left*spec_left[-1])
		edge_pts2.set_ydata(ly_right*spec_right[0])

		label.set_text(str(teffs[n])+" "+str(ratio)+" "+str(spec_left[-1]))
		return spec_plot


	if len(spec)==1:
		plt.savefig('m'+m+'q'+q+'.png')
	else:
		ani = animation.FuncAnimation(fig,update_img,len(spec),interval=30)
		writer = animation.writers['ffmpeg'](fps=1)
		ani.save('m'+m+'q'+q+'.mp4',writer=writer,dpi=dpi)

#return ani


# for i in range(len(spec)):
# 	spec_plot.set_ydata(spec[i,0]*spec[i,1])
# 	plt.draw()

def main():
	pat=sys.argv[1]
	spec_animate(pat)

if __name__ == '__main__':
    main()




