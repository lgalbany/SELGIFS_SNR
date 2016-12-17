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

#Incluyo: Hbeta+OIII4961+OIII5007. Ha 6563, 4861, 4961, 5007
#Hb is around 4895
#Primero, xbeta = x0 - 1715
#Y ahora saco las demas: xOIIIL = x0 - 1649
#xOIIIR = x0 - 1603
#
#La cuarta gaussiana la voy a situar como la de Ha pero respecto a Hb: xbeta4G = x0 - 1722

#Funcion de ajuste para la linea de H alfa
def gaussiana(x, A, B, C):
    """Modelo para nuestros datos."""
    return A * np.exp(-(x-x0)**2/(2.0*B)) + C


def gauss4(x, A1, B1, A2, B2, A3, B3, A4, B4, C4):
    """Modelo para nuestros datos."""
    return A1 * np.exp(-(x-x0)**2/(2.0*B1)) + A2 * np.exp(-(x-(x0-15))**2/(2.0*B2)) + A3 * np.exp(-(x-(x0+21))**2/(2.0*B3)) + A4 * np.exp(-(x-(x0-7))**2/(2.0*B4)) + C4


def gauss4b(x, A1b, B1b, A2b, B2b, A3b, B3b, A4b, B4b, C4b):
    """Modelo para nuestros datos."""
    return A1b * np.exp(-(x-(x0-1714))**2/(2.0*B1b)) + A2b * np.exp(-(x-(x0-1614))**2/(2.0*B2b)) + A3b * np.exp(-(x-(x0-1566))**2/(2.0*B3b)) + A4b * np.exp(-(x-(x0-1721))**2/(2.0*B4b)) + C4b


#Aqui viene lo gordo:
def fitting4gauss(flux,n,m):
	#Saca el espectro en ese spaxel
	spect1 = flux[:,n,m] #Remember (F,y,x)
	#Checkpoint1
	#print "Longitud del vector flujo (check point): "
	#print len(spect1)
	#Salto entre dato y dato en wavelength
    
    ######
    ###### ESTO HAY QUE LEERLO DEL HEADER, PARA CUBOS QUE NO SON DE CALIFA NO FUNCIONARÍA
    ###### TAMBIEN, LA POSICIÓN DE Halpha y Hbeta HAY QUE LEERLA DEL CUBO
    ###### EL KEYWORD MED_VEL TE DA LA VELOCIDAD DE RECESIÓN QUE PUEDES
    ###### CONVERTIR A REDSHIFT, LAS LINEAS ESTAN EN LAMBDA_EM*(1+Z)
    ######
    
	salto = (7501.0-3749.0)/len(spect1)
	#Limites en la longitud de onda (sacado de la web de CALIFA) Angstroms
	lamb0 = 3749.0
	lambf = 7501.0
	lamb = np.linspace(lamb0, lambf, len(spect1)) #Defino vector en ese tramo
	#Ventana de datos para ajustar linea - Necesito algo generico para cuidar la cinematica
	lamb00 = 3749 #Posicion 0 en el array tiene esta longitud de onda

    ######## PORQUE 1400 ?????? 

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
		P0 = Pini + 45

	#print "Posicion de la linea de H alpha"
	#print P0
	global x0 #Esto es importante para que se ejecute en las dos funciones de arriba
	x0 = P0*salto + lamb00 + 2 #El 2 es necesario porque soy un noob de python

	x0b = x0 - 1714 #Posicion de la linea de Hb
    
    ######## PORQUE 1714 ??????
    
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

	#Ajusto las 4 gaussianas para Hb:
	#Ventana de datos para ajustar las 3 lineas
	Pini4b = P0 - 900
	Pfin4b = P0 - 750
	#print Pini4
	#print spect1
	spect1v4b= spect1[Pini4b:Pfin4b]
	#print spect1v4
	#print spect1v4b

	lamb0v4b = Pini4b*salto + lamb00
	lambfv4b = Pfin4b*salto + lamb00
	difv4b = Pfin4b-Pini4b #Ventana espectro
	lambv4b = np.linspace(lamb0v4b, lambfv4b, num=difv4b)
	#print (x0 - lamb00)/salto
	#print "Ventana considerada para Halpha y los N"

	#print spect1v4
	A1=[] #Todos estos vectores son ncesarios para rellenarlos con los resultados del bucle
	B1=[]	
	A2=[]
	B2=[]
	A3=[]
	B3=[]
	A4=[]
	B4=[]
	C4=[]
	A1b=[] #Todos estos vectores son ncesarios para rellenarlos con los resultados del bucle
	B1b=[]	
	A2b=[]
	B2b=[]
	A3b=[]
	B3b=[]
	A4b=[]
	B4b=[]
	C4b=[]
	param_bounds=([0,0,0,0,0,0,0,120.0,-np.inf],[100,1000,100,1000,100,1000,100,10000,np.inf])#Estos limites son para que no rompa por exceso de tiempo en el ajuste. El 100 en el limite inferior es importante para que sea una SN
	param_boundsb=([0,0,0.05,0,0.1,5,0,0,-1],[400,600,400,600,400,600,300,1000,10])#Estos limites son para que no rompa por exceso de tiempo en el ajuste. El 100 en el limite inferior es importante para que sea una SN	
	if spect1v4[np.argmax(spect1v4)]>=0.4: #Esto es, si el valor donde esta el maximo es considerable, tira para delante.
		(A1, B1, A2, B2, A3, B3, A4, B4, C4), pcov = curve_fit(gauss4, lambv4, spect1v4, bounds=param_bounds)
		try:
			(A1b, B1b, A2b, B2b, A3b, B3b, A4b, B4b, C4b), _ = curve_fit(gauss4b, lambv4b, spect1v4b, bounds=param_boundsb)
		except RuntimeError:
			print "Error - Fitting in H beta failed (RuntimeError). Position: ", m, n
			A1b = 0
			B1b = 0
			A2b = 0
			B2b = 0
			A3b = 0
			B3b = 0
			A4b = 0
			B4b = 1000
			C4b = 0

		#print "Parametros del fit"
		#print (A1, B1, A2, B2, A3, B3, A4, B4, C4)
		#print "EW Halpha"
		#print EW
		Param1 = 1.2*A2
		Param2 = 2*ruido
		if B4>=150 and A4>=Param1 and A4>=Param2 and B4<=300: #Esto es para identificar las SN
			#print "-----------------------"
			#print "There is SN emission."
			#print "Position: ( m , n )=(",m,",",n,")"
			SN = 1 #Identificador de supernova
			print "------"
			print "Parameters of the peaks in beta area: "
			print A1b, B1b
			print A2b, B2b
			print A3b, B3b
			print A4b, B4b
			print "------"
			#print lambv4b
			return (SN, m, n, lambv4, spect1v4, A1, B1, A2, B2, A3, B3, A4, B4, C4, x0, A1b, B1b, A2b, B2b, A3b, B3b, A4b, B4b, C4b, x0b, lambv4b, spect1v4b)
		else:
			SN = 0 #Esto es que no hay SN
			B4 = 10000
			B4b = 1000
			return (SN, m, n, lambv4, spect1v4, A1, B1, A2, B2, A3, B3, A4, B4, C4, x0, A1b, B1b, A2b, B2b, A3b, B3b, A4b, B4b, C4b, x0b, lambv4b, spect1v4b)
	else:
		SN = 0 #Esto es que no hay SN
		#print "Position: ( m , n )=(",m,",",n,")"
		#print "No enough signal to fit. Probably you are far from the galactic center."
		A1 = 0
		B1 = 0
		A2 = 0
		B2 = 0
		A3 = 0
		B3 = 0
		A4 = 0
		B4 = 10000
		C4 = 0
		A1b = 0
		B1b = 0
		A2b = 0
		B2b = 0
		A3b = 0
		B3b = 0
		A4b = 0
		B4b = 1000
		C4b = 0
		return (SN, m, n, lambv4, spect1v4, A1, B1, A2, B2, A3, B3, A4, B4, C4, x0, A1b, B1b, A2b, B2b, A3b, B3b, A4b, B4b, C4b, x0b, lambv4b, spect1v4b)


	

