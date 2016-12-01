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
#SN identification. Col. Lluis Galbany
#Parte de CALIFA (Comun para Javi y Lluis)
#Cargamos datos.
#califa_cube = pyfits.open('NGC2906.V500.rscube.fits')
#califa_cube.info()
#flux = califa_cube[0].data #Nos quedamos con el PRIMARY
#Defino posicion espacial mxn
#m = 34; #Range from 0 to 77
#n = 47; #Range from 0 to 71

#Funcion de ajuste para la linea de H alfa
def gaussiana(x, A, B, C):
    """Modelo para nuestros datos."""
    return A * np.exp(-(x-x0)**2/(2.0*B)) + C


def gauss4(x, A1, B1, A2, B2, A3, B3, A4, B4, C4):
    """Modelo para nuestros datos."""
    return A1 * np.exp(-(x-x0)**2/(2.0*B1)) + A2 * np.exp(-(x-(x0-15))**2/(2.0*B2)) + A3 * np.exp(-(x-(x0+21))**2/(2.0*B3)) + A4 * np.exp(-(x-(x0-7))**2/(2.0*B4)) + C4

#Aqui viene lo gordo:
def fitting4gauss(flux,n,m):
	#Saca el espectro en ese spaxel
	spect1 = flux[:,n,m] #Remember (F,y,x)
	#Checkpoint1
	#print "Longitud del vector flujo (check point): "
	#print len(spect1)
	#Salto entre dato y dato en wavelength
	salto = (7501.0-3749.0)/len(spect1)
	#Limites en la longitud de onda (sacado de la web de CALIFA) Angstroms
	lamb0 = 3749.0
	lambf = 7501.0
	lamb = np.linspace(lamb0, lambf, len(spect1)) #Defino vector en ese tramo

	#Ventana de datos para ajustar linea - Necesito algo generico para cuidar la cinematica
	lamb00 = 3749 #Posicion 0 en el array tiene esta longitud de onda
	P = 1400
	Pini = P - 30
	Pfin = P + 30
	#Saco el flujo y la longitud de onda
	spect1max= spect1[Pini:Pfin] #Vector de flujo en ventana
	#print spect1max
	#Posicion de la linea de Halpha en este (mxn)
	PK = np.argmax(spect1max) #Posicion del maximo en la ventana, le sumo Pini
	P0 = PK + Pini +1 #Posicion del maximo en el espectro (el 1 es necesario)
	#print P0
	#Espero que sea siempre la linea de H alfa..........
	ruido = 0.1 #Esto es para que no me use maximos falsos o feos
	if spect1max[PK]<ruido: 
		P0 = Pini + 1

	#print "Posicion de la linea de H alpha"
	#print P0
	global x0 #Esto es importante para que se ejecute en las dos funciones de arriba
	x0 = P0*salto + lamb00 + 2 #El 2 es necesario porque soy un noob de python
	#print x0
	#Aqui comento la ventana para fitear solo la linea de H alpha
	#Ventana de datos para ajustar H alpha y cercanias.
	#Pini1 = P0 - 6
	#Pfin1 = P0 + 6
	#lambdaemit = 6629 #Lambda in the galactic center
	#zrot = (x0 - lambdaemit)/lambdaemit

	#Saco el flujo y la longitud de onda
	#spect1v1= spect1[Pini1:Pfin1] #Vector de flujo en ventana
	#lamb0v1 = Pini1*salto + lamb00
	#lambfv1 = Pfin1*salto + lamb00
	#difv1 = Pfin1-Pini1 #Ventana espectro
	#lambv1 = np.linspace(lamb0v1, lambfv1, num=difv1)

	#print "Ventana considerada para Halpha"
	#print Pini1*salto + lamb00
	#print Pfin1*salto + lamb00

	
#	(A, B, C), pcov = curve_fit(gaussiana, lambv1, spect1v1)

#	print "Parametros del fit"
#	print(A, B, C)
	#print pcov
	#perr = np.sqrt(np.diag(pcov)) # one standard deviation errors on the parameters
	#print perr

	#Ajusto las 4 gaussianas:
	#Ventana de datos para ajustar las 3 lineas
	Pini4 = P0 - 25
	Pfin4 = P0 + 30
	#print Pini4
	#print spect1
	spect1v4= spect1[Pini4:Pfin4]
	#print spect1v4

	lamb0v4 = Pini4*salto + lamb00
	lambfv4 = Pfin4*salto + lamb00
	difv4 = Pfin4-Pini4 #Ventana espectro
	lambv4 = np.linspace(lamb0v4, lambfv4, num=difv4)
	#print (x0 - lamb00)/salto
	#print "Ventana considerada para Halpha y los N"

	#print lamb0v4
	#print lambfv4

	#print spect1v4
	ajA4=[] #Todos estos vectores son ncesarios para rellenarlos con los resultados del bucle
	ajB4=[]
	ajC4=[]
	ajD4=[]
	C4=[]
	EW1=[]
	EW2=[]
	EW3=[]
	EW4=[]
	if spect1v4[np.argmax(spect1v4)]>=0.4: #Esto es, si el valor donde esta el maximo es considerable, tira para delante.
		param_bounds=([0,0,0,0,0,0,0,100.0,-np.inf],[5000,5000,5000,5000,5000,5000,5000,10000,np.inf])#Estos limites son para que no rompa por exceso de tiempo en el ajuste. El 100 en el limite inferior es importante para que sea una SN
		(A1, B1, A2, B2, A3, B3, A4, B4, C4), pcov = curve_fit(gauss4, lambv4, spect1v4, bounds=param_bounds)
		#print "Parametros del fit"
		#print (A1, B1, A2, B2, A3, B3, A4, B4, C4)

		ajA4=A1 * np.exp(-(lambv4-x0)**2/(2.0*B1)) #Primera gaussiana, picada en x0
		ajB4=A2 * np.exp(-(lambv4-(x0-15))**2/(2.0*B2)) #Segunda gaussiana, a la izquierda
		ajC4=A3 * np.exp(-(lambv4-(x0+21))**2/(2.0*B3)) #Tercera gaussiana, a la derecha
		ajD4=A4 * np.exp(-(lambv4-(x0-7))**2/(2.0*B4)) #Cuarta gaussiana, entre Halpha y [NII]
		ajuste4= ajA4 + ajB4 + ajC4 + ajD4 + C4
		#Calculamos el EW de las lineas
		EW1 = (1-C4/(A1+C4))*np.sqrt(2*np.pi*B1)
		EW2 = (1-C4/(A2+C4))*np.sqrt(2*np.pi*B2)
		EW3 = (1-C4/(A3+C4))*np.sqrt(2*np.pi*B3)
		EW4 = (1-C4/(A4+C4))*np.sqrt(2*np.pi*B4)
		#print "EW Halpha"
		#print EW
		Param1 = 2*A2
		Param2 = 2*ruido
		if B4>150 and A4>=Param1 and A4>Param2 and B4<300: #Esto es para identificar las SN
			#print "-----------------------"
			#print "There is SN emission."
			#print "Position: ( m , n )=(",m,",",n,")"
			SN = 1 #Identificador de supernova
			return (SN, m, n, lambv4, spect1v4, ajA4, ajB4, ajC4, ajD4, C4, x0, EW1, EW2, EW3, EW4)
		else:
			SN = 0 #Esto es que no hay SN
			return (SN, m, n, lambv4, spect1v4, ajA4, ajB4, ajC4, ajD4, C4, x0, EW1, EW2, EW3, EW4)
	else:
		SN = 0 #Esto es que no hay SN
		#print "Position: ( m , n )=(",m,",",n,")"
		#print "No enough signal to fit. Probably you are far from the galactic center."
		ajA4 = np.linspace(0,0,len(lambv4))
		ajB4 = np.linspace(0,0,len(lambv4))
		ajC4 = np.linspace(0,0,len(lambv4))
		ajD4 = np.linspace(0,0,len(lambv4))
		C4 = np.linspace(0,0,len(lambv4))
		return (SN, m, n, lambv4, spect1v4, ajA4, ajB4, ajC4, ajD4, C4, x0, EW1, EW2, EW3, EW4)


	

