Repositorio para el proyecto de SNR de la escuela de SELGIFS.

Código actualizado a fecha de 7 de diciembre.

Se debe ejecutar el script VPSNfinder.py estando en la misma carpeta que fitting4gauss.py y que el cubo de datos de CALIFA de la galaxia NGC 2906.

El programa localiza los spaxels con características de supernova en la galaxia considerada, y saca los ajustes del espectro en las regiones de la línea de Ha (con los respectivos [NII]) y de la línea Hb (junto con los [OIII]). La posicion de Ha la saca localizando el máximo de una región próxima a donde podemos encontrar esa línea. Luego el resto de posiciones vienen descritos en función de este valor. La suposición es que el máximo en esa región coincide con la línea de Ha, y de momento no ha tenido problemas. En las regiones de SN hacemos un plot de uno de los spaxels, y vemos que hay un exceso de emisión al lado de las líneas Ha-Hb, que ajustamos con una cuarta gaussiana.

Por terminal nos saca las localizaciones (m,n) de cada spaxel considerado SN, y hasta que lo modifique (a modo de ayuda visual) he incluido que me saque los parámetros del ajuste entorno a la línea de Hb.

El programa también saca un mapa de la posición de Ha, simulando un mapa de velocidades, y otro de sigma cuadrado (la varianza) de la cuarta gaussiana, con lo que vemos que solo se marca la zona que coincide con la supernova.

De momento eso es todo hasta la vuelta.

