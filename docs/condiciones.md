## Condiciones Iniciales
En la ecuacion de calor en 2D, las condiciones inciales definen como esta distribuida la temperatura en todo el dominio espacial en un tiempo inicial. Las condiciones inciales y de frontera son esenciales porque determinan completamente la evolucion temporal de la temperatura. 
Las condiciones iniciales que se utilizaron para la resolucion de la ecuacion de calor en 2D son: 

**1) Campana Gaussiana Centrada:**

    u[:, :] = np.exp(-100 * ((X - 0.5)**2 + (Y - 0.5)**2)) 

Esta condicion simula una fuente de calor muy localizada, la temperatura es mas alta en el centro y disminuye rapidamente hacia los bordes. EL numero 100 en el exponente controla la concentracion de calor. Por lo que la campana Gaussiana centrada simula un calor concentrado, el cual es ideal para estudiar como se difunde el calor desde un punto. 

**2) Paraboloide Centrado:**

    u[:, :] = 10 * ((X - 0.5)**2 + (Y - 0.5)**2)

En esta condicion inicial, el centro esta mas frio y la temperatura aumenta hacia los bordes, como si las paredes externas calentaran la placa. El 10 es un factor de escala que controla que tan caliente estan los bordes. EL Paraboloide centrado es bueno para estudiar el flujo del calor del exterior hacia el centro. Cabe resaltar que esta condicion inicial es simetrica radialmente. 

**3) Onda Senosoidal:**

    u[:, :] = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)

La onda Senoidal simula el patron oscilante de calor; en donde la temperatura del dominio subiera y bajara de forma periodica. Esta condicion representa una superposicion de modos termicos y es muy util para ver como se disipan modos oscilatorios de calor en el tiempo. 
    
----------------------------------------------------------------------------------------------------------------------

## Condiciones de Frontera 

Al resolver la ecuacion de calor en 2D, se debe de decirle al programa que es lo que pasa con los bordes de la placa, y esas son las condiciones de frontera. A continuacion, se explicaran las condiciones de frontera que se utilizaron para resolver dicha ecuacion: 

**1) Dirichlet:**

 Dirichlet simula una placa conectada a un material o liquido que absorbe todo el calor que llega, los bordes del dominio estan en contacto con reservorios termicos que los mantienenm a una temperatura fija. Esta condicion de frontera es muy estable y comun en simulaciones.

    u_proxima[0, :] = 0; u_proxima[-1, :] = 0; u_proxima[:, 0] = 0; u_proxima[:, -1] = 0


**2) Neumann:**

Neumann simula que los bordes estan aislados termicamente, el calor no puede salir ni entrar por los bordes y no hay flujo de calor a traves de los bordes. Esto implica que no se conoce la temperatura en la frontera, sino el flujo de calor a través de ella. En este caso, lo que se fija es la tasa de variación de la temperatura en dirección perpendicular al borde del dominio. 
 
    #u_proxima[0, :] = u_proxima[1, :]
    #u_proxima[-1, :] = u_proxima[-2, :]
    #u_proxima[:, 0] = u_proxima[:, 1]
    #u_proxima[:, -1] = u_proxima[:, -2]


**3) Robin:**

Robin es una condición mixta de Neumann y Dirichlet que modela un intercambio de calor con el ambiente exterior. El parametro beta es el coeficiente de conveccion, por lo que si beta es grande, se pierde calor mas rapido y si beta es pequeño, se pierde calor lentamente. Se establece una relación entre la temperatura en la frontera y su derivada normal.

    beta = 3.0 # Se puede modificar el valor  (Coeficiente de transferencia de calor en fronteras [W/m²K])
    u_proxima[0, :] = u_proxima[1, :] / (1 + beta * dx)
    u_proxima[-1, :] = u_proxima[-2, :] / (1 + beta * dx)
    u_proxima[:, 0] = u_proxima[:, 1] / (1 + beta * dy)
    u_proxima[:, -1] = u_proxima[:, -2] / (1 + beta * dy)



