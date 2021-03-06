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

Result=[]
Resultb=[]
vx0=[] #defino los vectores a rellenar antes de los bucles
vx0b=[]
PosSN=[]
FittingSN=[]
FittingSNb=[]
Spectras=[]
SpectrasSN=[]
Spectrasb=[]
SpectrasSNb=[]

for m in range(78):
	for n in range(72):
		SNdef = fitting4gauss.fitting4gauss(flux,n,m)
		#SN, m, n, lambv4, spect1v4, A1, B1, A2, B2, A3, B3, A4, B4, C4, x0, A1b, B1b, A2b, B2b, A3b, B3b, A4b, B4b, C4b, x0b, lambv4b, spect1v4b
		Result.append([SNdef[0], SNdef[1], SNdef[2], SNdef[5], SNdef[6], SNdef[7], SNdef[8], SNdef[9], SNdef[10], SNdef[11], SNdef[12], SNdef[13], SNdef[14]])
		Resultb.append([SNdef[0], SNdef[1], SNdef[2], SNdef[15], SNdef[16], SNdef[17], SNdef[18], SNdef[19], SNdef[20], SNdef[21], SNdef[22], SNdef[23], SNdef[24]])
		Spectras.append([SNdef[3], SNdef[4]])
		Spectrasb.append([SNdef[25], SNdef[26]])
		vx0.append([SNdef[14], SNdef[1], SNdef[2]]) #Posiciones de la linea H alpha (x0,m,n)
		vx0b.append([SNdef[24], SNdef[1], SNdef[2]])
		if SNdef[0]==1: #Si es una supernova, saca todo para los plots
			PosSN.append([SNdef[1], SNdef[2]]) #Posiciones SN
			FittingSN.append([SNdef[0], SNdef[1], SNdef[2], SNdef[5], SNdef[6], SNdef[7], SNdef[8], SNdef[9], SNdef[10], SNdef[11], SNdef[12], SNdef[13], SNdef[14]])
			FittingSNb.append([SNdef[0], SNdef[1], SNdef[2], SNdef[15], SNdef[16], SNdef[17], SNdef[18], SNdef[19], SNdef[20], SNdef[21], SNdef[22], SNdef[23], SNdef[24]])
			SpectrasSN.append([SNdef[3], SNdef[4]])
			SpectrasSNb.append([SNdef[25], SNdef[26]])
			mm = SNdef[1] #Posicion m
			nn = SNdef[2] #Posicion n
			print "There is a SN in this Galaxy"
			print "Position (m,n): ",mm,nn

print "---------------------------------------------------------------------------"
Resultados = np.asarray(Result) #Esto es para pasar de lista a array
Resultadosb = np.asarray(Resultb) #Esto es para pasar de lista a array

Posit = np.asarray(vx0) #Esto es para pasar de lista a array
Positbeta = np.asarray(vx0b) #Esto es para pasar de lista a array

PosicionesSN = np.asarray(PosSN) #Esto es para pasar de lista a array
AjusteSN = np.asarray(FittingSN) #Esto es para pasar de lista a array
AjusteSNb = np.asarray(FittingSNb) #Esto es para pasar de lista a array

Ventana = np.asarray(SpectrasSN)
Ventanab = np.asarray(SpectrasSNb)

#print "vx0"
#print vx0

#print "FittingSN"
#print FittingSN

#print "Ventana"
#print Ventana

#np.set_printoptions(threshold=np.nan) #Esto es para el print siguiente, para ver el vector entero
#print Posit[:,0] #Para comprobar lo que salia aqui

ajA4=[]
ajB4=[]
ajC4=[]
ajD4=[]
ajA4b=[]
ajB4b=[]
ajC4b=[]
ajD4b=[]

for i in range(len(AjusteSN[:,0])):
	ajA4.append(AjusteSN[i,3] * np.exp(-(Ventana[i,0,:]-AjusteSN[i,12])**2/(2.0*AjusteSN[i,4]))) #Primera gaussiana, picada en x0
	ajB4.append(AjusteSN[i,5] * np.exp(-(Ventana[i,0,:]-(AjusteSN[i,12]-15))**2/(2.0*AjusteSN[i,6]))) #Segunda gaussiana, a la izquierda
	ajC4.append(AjusteSN[i,7] * np.exp(-(Ventana[i,0,:]-(AjusteSN[i,12]+21))**2/(2.0*AjusteSN[i,8]))) #Tercera gaussiana, a la derecha
	ajD4.append(AjusteSN[i,9] * np.exp(-(Ventana[i,0,:]-(AjusteSN[i,12]-7))**2/(2.0*AjusteSN[i,10]))) #Cuarta gaussiana, entre H alpha y [NII]

for i in range(len(AjusteSNb[:,0])):
	ajA4b.append(AjusteSNb[i,3] * np.exp(-(Ventanab[i,0,:]-AjusteSNb[i,12])**2/(2.0*AjusteSNb[i,4]))) #Primera gaussiana, picada en x0
	ajB4b.append(AjusteSNb[i,5] * np.exp(-(Ventanab[i,0,:]-(AjusteSNb[i,12]+100))**2/(2.0*AjusteSNb[i,6]))) #Segunda gaussiana, a la izquierda
	ajC4b.append(AjusteSNb[i,7] * np.exp(-(Ventanab[i,0,:]-(AjusteSNb[i,12]+148))**2/(2.0*AjusteSNb[i,8]))) #Tercera gaussiana, a la derecha
	ajD4b.append(AjusteSNb[i,9] * np.exp(-(Ventanab[i,0,:]-(AjusteSNb[i,12]-7))**2/(2.0*AjusteSNb[i,10]))) #Cuarta gaussiana, entre H alpha y [NII]
	

#Calculamos el EW de las lineas
#EW1 = (1-Resultados[:,11]/(Resultados[:,3]+Resultados[:,11]))*np.sqrt(2*np.pi*Resultados[:,4])
#EW2 = (1-Resultados[:,11]/(Resultados[:,5]+Resultados[:,11]))*np.sqrt(2*np.pi*Resultados[:,6])
#EW3 = (1-Resultados[:,11]/(Resultados[:,7]+Resultados[:,11]))*np.sqrt(2*np.pi*Resultados[:,8])
#EW4 = (1-Resultados[:,11]/(Resultados[:,9]+Resultados[:,11]))*np.sqrt(2*np.pi*Resultados[:,10])

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

#Las 4 gaussianas. Comparten el C4.
gauss1=[]
gauss2=[]
gauss3=[]
gauss4=[]
ajuste4=[]

gauss1b=[]
gauss2b=[]
gauss3b=[]
gauss4b=[]
ajuste4b=[]

for i in range(len(AjusteSN[:,0])):
	gauss1.append(ajA4[i] + AjusteSN[i,11])
	gauss2.append(ajB4[i] + AjusteSN[i,11])
	gauss3.append(ajC4[i] + AjusteSN[i,11])
	gauss4.append(ajD4[i] + AjusteSN[i,11])
	ajuste4.append(ajA4[i] + ajB4[i] + ajC4[i] + ajD4[i] + AjusteSN[i,11]) #La suma de las 4 mas el C4


for i in range(len(AjusteSNb[:,0])):
	gauss1b.append(ajA4b[i] + AjusteSNb[i,11])
	gauss2b.append(ajB4b[i] + AjusteSNb[i,11])
	gauss3b.append(ajC4b[i] + AjusteSNb[i,11])
	gauss4b.append(ajD4b[i] + AjusteSNb[i,11])
	ajuste4b.append(ajA4b[i] + ajB4b[i] + ajC4b[i] + ajD4b[i] + AjusteSNb[i,11]) #La suma de las 4 mas el C4
#print "AjusteSN"
#print AjusteSN

longs=Ventana[:,0,:] #Vector de vectores de longitud de onda
spects=Ventana[:,1,:] #Vector de spectros

longsb=Ventanab[:,0,:] #Vector de vectores de longitud de onda
spectsb=Ventanab[:,1,:] #Vector de spectros

lambv4 = longs[3,:] #Vector de longitudes de onda donde elijo uno particular (de momento)
spect1v4 = spects[3,:] #Vector de flujos

lambv4b = longsb[3,:] #Vector de longitudes de onda donde elijo uno particular (de momento)
spect1v4b = spectsb[3,:] #Vector de flujos

#Primer plot, el del fit.
fig1=plt.figure(num=1,figsize=(8,5))
plt.plot(lambv4, spect1v4, label='Spectra')
plt.plot(lambv4, ajuste4[3], label='4 Gaussians')
plt.plot(lambv4, gauss1[3], label='Gaussian Halpha')
plt.plot(lambv4, gauss2[3], label='Gaussian [NII]')
plt.plot(lambv4, gauss3[3], label='Gaussian [NII]')
plt.plot(lambv4, gauss4[3], label='Gaussian SN')
plt.xlabel('Wavelength')
plt.ylabel('Flux')
plt.title('NGC 2906 (Ha and 2[NII] fitting 4 Gaussians)')
plt.legend(loc=2, fontsize = 'small')
plt.savefig('NGC2906spectraHaNs4Gaussianas')
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
#print "Todo"
#np.set_printoptions(threshold=np.nan)
#print Posit[:,0]

#Segundo plot, el de campo de posiciones de la linea de H alpha.
fig2=plt.figure(num = 2, figsize = (10, 8), facecolor="white")
cm = plt.cm.get_cmap('seismic')
sc = plt.scatter(Posit[:,1],Posit[:,2], c=Posit[:,0], vmin=6580, vmax=6610, marker='s', s=20, cmap=cm, label='H alpha Position')
sc2 = plt.scatter(PosicionesSN[:,0],PosicionesSN[:,1], c='g', marker='s', s=20, cmap=cm, label='SN Position')
plt.legend(loc = 2)
plt.colorbar(sc)
plt.savefig('Posiciones.pdf')
fig1.show()
fig2.show()
plt.close()

#Tercer plot, la varianza.
fig3=plt.figure(num = 3, figsize = (10, 8), facecolor="white")
cm = plt.cm.get_cmap('seismic')
sc = plt.scatter(Resultados[:,1],Resultados[:,2], c=Resultados[:,10], vmin=np.min(Resultados[:,10]), vmax=np.max(Resultados[:,10]), marker='s', s=20, cmap=cm, label='Variance value 4G Ha')
plt.legend(loc = 2)
plt.colorbar(sc)
plt.savefig('FitsB4.pdf')
fig1.show()
fig2.show()
fig3.show()
plt.close()

fig4=plt.figure(num = 4, figsize = (10, 8), facecolor="white")
cm = plt.cm.get_cmap('seismic')
sc = plt.scatter(Resultados[:,1],Resultados[:,2], c=Resultados[:,4], vmin=np.min(Resultados[:,4]), vmax=np.max(Resultados[:,4]), marker='s', s=20, cmap=cm, label='Variance value Ha')
plt.legend(loc = 2)
plt.colorbar(sc)
plt.savefig('FitsVarHa.pdf')
fig1.show()
fig2.show()
fig3.show()
fig4.show()
plt.close()

fig5=plt.figure(num=5,figsize=(8,5))
plt.plot(lambv4b, spect1v4b, label='Spectra')
plt.plot(lambv4b, ajuste4b[3], label='4 Gaussians')
plt.plot(lambv4b, gauss1b[3], label='Gaussian Hbeta')
plt.plot(lambv4b, gauss2b[3], label='Gaussian [OIII]')
plt.plot(lambv4b, gauss3b[3], label='Gaussian [OIII]')
plt.plot(lambv4b, gauss4b[3], label='Gaussian SN (beta)')
plt.xlabel('Wavelength')
plt.ylabel('Flux')
plt.title('NGC 2906 (Hb and 2 [OIII] fitting 4 Gaussians)')
plt.legend(loc=2, prop={'size':10})
plt.savefig('NGC2906spectraHbOs4Gaussianas')
fig1.show()
fig2.show()
fig3.show()
fig4.show()
fig5.show()
plt.close()

fig6=plt.figure(num = 6, figsize = (10, 8), facecolor="white")
cm = plt.cm.get_cmap('seismic')
sc = plt.scatter(Resultadosb[:,1],Resultadosb[:,2], c=Resultadosb[:,4], vmin=np.min(Resultadosb[:,4]), vmax=np.max(Resultadosb[:,4]), marker='s', s=20, cmap=cm, label='Variance value Hb (Var)')
plt.legend(loc = 2)
plt.colorbar(sc)
plt.savefig('FitsVarHb.pdf')
fig1.show()
fig2.show()
fig3.show()
fig4.show()
fig5.show()
fig6.show()
plt.close()

fig7=plt.figure(num = 7, figsize = (10, 8), facecolor="white")
cm = plt.cm.get_cmap('seismic')
sc = plt.scatter(Resultadosb[:,1],Resultadosb[:,2], c=Resultadosb[:,10], vmin=np.min(Resultadosb[:,10]), vmax=np.max(Resultadosb[:,10]), marker='s', s=20, cmap=cm, label='Variance value 4Gauss Hb')
plt.legend(loc = 2)
plt.colorbar(sc)
plt.savefig('FitsB4Hb.pdf')
fig1.show()
fig2.show()
fig3.show()
fig4.show()
fig5.show()
fig6.show()
fig7.show()
plt.close()

fig8=plt.figure(num=8,figsize=(8,5))
plt.plot(lambv4b, spect1v4b, label='Spectra')
plt.plot(lambv4b, ajuste4b[6], label='4 Gaussians')
plt.plot(lambv4b, gauss1b[6], label='Gaussian Hbeta')
plt.plot(lambv4b, gauss2b[6], label='Gaussian [OIII]')
plt.plot(lambv4b, gauss3b[6], label='Gaussian [OIII]')
plt.plot(lambv4b, gauss4b[6], label='Gaussian SN (beta)')
plt.xlabel('Wavelength')
plt.ylabel('Flux')
plt.title('NGC 2906 (Hb and 2 [OIII] fitting 4 Gaussians)')
plt.legend(loc=2, prop={'size':10})
plt.savefig('NGC2906spectraHbOs4Gaussianas2')
fig1.show()
fig2.show()
fig3.show()
fig4.show()
fig5.show()
fig6.show()
fig7.show()
fig8.show()
plt.close()

fig9=plt.figure(num=9,figsize=(8,5))
plt.plot(lambv4b, spect1v4b, label='Spectra')
plt.plot(lambv4b, ajuste4b[9], label='4 Gaussians')
plt.plot(lambv4b, gauss1b[9], label='Gaussian Hbeta')
plt.plot(lambv4b, gauss2b[9], label='Gaussian [OIII]')
plt.plot(lambv4b, gauss3b[9], label='Gaussian [OIII]')
plt.plot(lambv4b, gauss4b[9], label='Gaussian SN (beta)')
plt.xlabel('Wavelength')
plt.ylabel('Flux')
plt.title('NGC 2906 (Hb and 2 [OIII] fitting 4 Gaussians)')
plt.legend(loc=2, prop={'size':10})
plt.savefig('NGC2906spectraHbOs4Gaussianas3')
fig1.show()
fig2.show()
fig3.show()
fig4.show()
fig5.show()
fig6.show()
fig7.show()
fig8.show()
plt.show()
plt.close()


