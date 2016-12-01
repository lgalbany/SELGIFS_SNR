import numpy as np
import pyfits
import matplotlib.pyplot as plt
import pystarlight.io.starlighttable
import pystarlight
import atpy
import matplotlib.mlab as mlab
from astropy.io import fits
from scipy.misc import imread # Cargo imread de scipy.misc
from matplotlib import cm
from astropy.utils.data import download_file
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
from scipy.stats import norm

from mpl_toolkits.mplot3d import Axes3D         # Cargo Axes3D de mpl_toolkits.mplot3d
import matplotlib as mpl
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

import fitting4gauss

#SN identification. Colaborando con Lluis Galbany
#El codigo esta escrito en sucio porque tampoco contaba compartirlo hasta dar una version final.

#Parte de datos de CALIFA (Comun para Javi y Lluis)
#Cargamos datos.
cube = 'NGC2906.V500.rscube.fits'
califa_cube = pyfits.open(cube)
califa_cube.info()
flux = califa_cube[0].data #Nos quedamos con el PRIMARY

print "-------------------------------------------------------------"

print "Galaxy: ",cube
print " "
print "Possible SNs: "
print " "
#Esto era porque reciclo codigo y lo comento.
#Defino posicion espacial mxn
#m = 4; #Range from 0 to 77
#n = 4; #Range from 0 to 71

vx0=[] #defino los vectores a rellenar antes de los bucles
#zrot=[]

for m in range(78):
	for n in range(72):
		SNdef = fitting4gauss.fitting4gauss(flux,n,m)
		vx0.append([SNdef[10], SNdef[1], SNdef[2]]) #Posiciones de la linea H alpha (x0,m,n)
		if SNdef[0]==1: #Si es una supernova, saca todo para los plots
			mm = SNdef[1] #Posicion m
			nn = SNdef[2] #Posicion n
			lambv4 = SNdef[3] #Vector de longitudes de onda
			spect1v4 = SNdef[4] #Vector de flujos
			ajA4 = SNdef[5] #Primera gaussiana - H alpha
			ajB4 = SNdef[6] #Segunda gaussiana - [NII] de la izq
			ajC4 = SNdef[7] #Tercera gaussiana - [NII] de la der
			ajD4 = SNdef[8] #Cuarta gaussiana - La del exceso por SN
			C4 = SNdef[9] #El offset
			EW1 = SNdef[11] #EW de H alpha
			EW2 = SNdef[12] #EW de [NII] de la izq
			EW3 = SNdef[13] #EW de [NII] de la der
			EW4 = SNdef[14] #EW del exceso
			print "There is a SN in this Galaxy"
			print "Position (m,n): ",mm,nn
			#print "Equivalent widths of the lines"
			#print "H alpha"
			#print EW1
			#print "[NII] left"
			#print EW2
			#print "[NII] right"
			#print EW3
			#print "SN excess"
			#print EW4
			#print "  "

print "---------------------------------------------------------------------------"
Posit = np.asarray(vx0) #Esto es para pasar de lista a array
#np.set_printoptions(threshold=np.nan) #Esto es para el print siguiente, para ver el vector entero
#print Posit[:,0] #Para comprobar lo que salia aqui

#Las 4 gaussianas. Comparten el C4.
gauss1 = ajA4 + C4
gauss2 = ajB4 + C4
gauss3 = ajC4 + C4
gauss4 = ajD4 + C4
ajuste4= ajA4 + ajB4 + ajC4 + ajD4 + C4 #La suma de las 4 mas el C4

#Primer plot, el del fit.
fig1=plt.figure(num=1,figsize=(8,5))
plt.plot(lambv4, spect1v4, label='Spectra')
plt.plot(lambv4, ajuste4, label='4 Gaussians')
plt.plot(lambv4, gauss1, label='Gaussian Halpha')
plt.plot(lambv4, gauss2, label='Gaussian [NII]')
plt.plot(lambv4, gauss3, label='Gaussian [NII]')
plt.plot(lambv4, gauss4, label='Gaussian SN')
plt.xlabel('Wavelength')
plt.ylabel('Flux')
plt.title('NGC 2906 (Ha-N fitting 4 Gaussians)')
plt.legend(loc=2)
plt.savefig('NGC2906spectraHaNGaussiana4prueba')
fig1.show() #Pongo fig1 en lugar de plt porque asi se me ven todas juntas al final
plt.close()

#Calculamos el EW de la linea de Halpha.
#EW = (1-C/(A+C))*np.sqrt(2*np.pi*B)
#print "EW Halpha"
#print EW

#Checking
#print len(Posit[:,1])
#print Posit[:,1]
#print len(Posit[:,2])
#print Posit[:,2]
#print len(Posit[:,0])
#print Posit[:,0]

#print "Max"
#print np.max(Posit[:,0])
#print "Min"
#print np.min(Posit[:,0])

#Segundo plot, el de campo de posiciones de la linea de H alpha. Esta cutremente sacado.
#Si alguno sabe como fittearlo mejor se agradeceria.
fig2=plt.figure(num = 2, figsize = (10, 8), facecolor="white")
ax2 = fig2.gca(projection='3d')
pd = ax2.scatter(Posit[:,1],Posit[:,2],Posit[:,0], marker='o', s=10, cmap='seismic')
#plot2 = ax2.plot_surface(xx,yy,Posit[:,0], rstride = 2, cstride = 2, cmap = cm.coolwarm, linewidth = 0.5, antialiased = True)
m = cm.ScalarMappable(cmap='seismic')
m.set_array([6490,6515])
plt.colorbar(m)
ax2.view_init(elev=89.9, azim=-90.1)
ax2.axes.zaxis.set_ticklabels([])
plt.tight_layout()
plt.draw()
plt.savefig('Posiciones.pdf')
fig1.show()
plt.show()
plt.close()
