# Soluión en Python 

En el presente apartado explicaremos a detalle cada linea de codigo de la resolucion de la Ecuacion de calor realizada en Python.

El codigo tiene 3 diferentes condiciones inciales y 3 diferentes condiciones de frontera, por lo que para obtener los 3 graficos esperados, debemos de escoger una opcion para cada condicion inicial y de frontera. 

Iniciamos importando las librerias, en donde numpy nos ayuda a realizar calculos numericos y el manejo de arreglos/matirces; matplotlib.pyplot, grafica y visualiza datos; 
scipy.sparse.diags, permite crear matrices dispersas que ahorran memoria y son mas eficientes; scipy.sparse.linalg.spsolve, resuelve sistemas lineales que usan matrices dispersas; FuncAnimation, permite crear animaciones y Image, sirve para mostrar imagenes en entornos como Jupyter Notebook. 

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    from matplotlib.animation import FuncAnimation
    from IPython.display import Image
    
Parámetros; alpha es la constante de difusión térmica, controla que tan rapido se difunde el calor. 

    alpha = 1.0 
    
Configuración del dominio; se define un dominio cuadrado [0,1] x [0,1], se divide ese dominio en 20 segmentos en cada direccion, dx y dy corresponde al tamano de cada celda, es decir, el espaciado entre nodos. 

    Longitud_x, Longitud_y = 1.0, 1.0  # Estos son los largos del dominio en x e y respectivamente 
    num_divisiones_x, num_divisiones_y = 20, 20  # Cantidad de divisiones en las direcciones x e y
    
    #Tamaño de cada celda (Resolución espacial)
    dx, dy = Longitud_x / num_divisiones_x, Longitud_y / num_divisiones_y  # Distancia entre nodos consecutivos en cada eje

Malla de coordenadas; el np.linspace(...) genera nodos equiespaciados entre 0 y 1, X, Y son las matrices que representan la malla completa del dominio, las cuales son utiles para graficar. 

    x = np.linspace(0, Longitud_x, num_divisiones_x+1)  # Puntos equiespaciados de 0 a Lx (incluyendo extremos)
    y = np.linspace(0, Longitud_y, num_divisiones_y+1)  # Puntos equiespaciados de 0 a Ly (incluyendo extremos)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
Configuración temporal

    T = 0.2      # Tiempo total de simulación
    dt = 0.001   # Paso temporal
    nt = int(T / dt)  # Número total de pasos de tiempo
    
Cálculo de parámetros auxiliares para Crank-Nicolson en x en y; se divide entre 2 para la estabilidad y precision (Crank-Nicolson), estos parametros se usan para construir las matrices tridiagonales. 

    parametro_x = alpha * dt / (2 * dx**2)  # Parámetro de difusión en x, dividido entre 2 por Crank-Nicolson
    parametro_y = alpha * dt / (2 * dy**2)  # Igual pero en y
    
Construcción de matrices tridiagonales para el método implícito en x

    Ax = diags([[-parametro_x]*(num_divisiones_x-1), [1+2*parametro_x]*(num_divisiones_x-1), [-parametro_x]*(num_divisiones_x-1)], [-1,0,1], shape=(num_divisiones_x-1, num_divisiones_x-1))  # Matriz del lado izquierdo (implícito) en x
    Bx = diags([[parametro_x]*(num_divisiones_x-1), [1-2*parametro_x]*(num_divisiones_x-1), [parametro_x]*(num_divisiones_x-1)], [-1,0,1], shape=(num_divisiones_x-1, num_divisiones_x-1))    # Matriz del lado derecho (explícito) en x
    
Construcción de matrices tridiagonales para el método implícito en y
    
    Ay = diags([[-parametro_y]*(num_divisiones_y-1), [1+2*parametro_y]*(num_divisiones_y-1), [-parametro_y]*(num_divisiones_y-1)], [-1,0,1], shape=(num_divisiones_y-1, num_divisiones_y-1))  # Implícito en y
    By = diags([[parametro_y]*(num_divisiones_y-1), [1-2*parametro_y]*(num_divisiones_y-1), [parametro_y]*(num_divisiones_y-1)], [-1,0,1], shape=(num_divisiones_y-1, num_divisiones_y-1))    # Explícito en y
    
Condición Inicial; se inicializa toda la matriz u con ceros y luego la llena con una onda senoidal doble.

    u = np.zeros((num_divisiones_x+1, num_divisiones_y+1))  # Inicialización de la matriz de temperatura u en todo el dominio

Ahora, vienen las 3 opciones de condiciones inciales, las cuales son Pulso Gaussiano centrado, Paraboloide centrado y Onda senosoidal.

Opción 1: Pulso gaussiano centrado en (0.5, 0.5). Se crea un pulso de calor concentrado en el centro del dominio, X y Y son las coordenadas de cada punto en la malla, el factor 100 controla que tan centrado esta el calor. 

    u[:, :] = np.exp(-100 * ((X - 0.5)**2 + (Y - 0.5)**2))
    
Opción 2: Paraboloide centrado (Alternativa físicamente consistente). Esta centrado en (0.5, 0.5) y tiene su valor minimo en el centro. 

    u[:, :] = 10 * ((X - 0.5)**2 + (Y - 0.5)**2)
    
Opción 3: Onda senosoidal suave (Para patrones periódicos), se usa una onda senoidal en ambas direcciones y y y, el 2 pi garantiza que la onda tiene un ciclo completo entre 0 y 1. 

    u[:, :] = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)
    
Matriz para almacenar la próxima temperatura

    u_proxima = np.zeros_like(u) # Matriz que contendrá la temperatura actualizada
    
Simulacion
    
Inicializacion, esta es la lista para guardar estados intermedios.

    historial_temperaturas = []
    intervalo_guardado = 10 # Guarda cada 10 pasos de tiempo
    
Bucle de tiempo, comienza el lazo de integración temporal; el for n in range(nt) itera desde t=0 hasta t=T

    for n in range(nt):
    
Paso intermedio: resolvemos en la dirección x, manteniendo y fijo. Es decir, fijamos una fila j y resolvemos en x, se forma el lado derecho con Bx y se resuelve el sistema lineal con Ax. 

        u_intermedia = np.zeros_like(u)  # Matriz temporal para almacenar resultados intermedios
        # Recorremos cada fila fija (j) y resolvemos en x (columnas)
        for j in range(1, num_divisiones_y):
            rhs = Bx.dot(u[1:num_divisiones_x, j])              # Lado derecho del sistema: combinación explícita
            u_intermedia[1:num_divisiones_x, j] = spsolve(Ax, rhs)    # Resolvemos el sistema lineal en x
    
Paso final: resolvemos en la dirección y, manteniendo x fijo. 

        for i in range(1, num_divisiones_x):
            rhs = By.dot(u_intermedia[i, 1:num_divisiones_y])  # Lado derecho para dirección y
            u_proxima[i, 1:num_divisiones_y] = spsolve(Ay, rhs)   # Solución del sistema lineal en y
    
Aplicar la Condición de Frontera seleccionada a u_proxima. Nuevamente, vienen 3 diferentes opciones de condiciones de frontera, Dirichlet, Neumann y Robin. 
        
Opción 1: Dirichlet (bordes fijos en 0); se fija la temperatura en los bordes a un valor constante y representa un material en contacto con un medio frio constante. 

        u_proxima[0, :] = 0; u_proxima[-1, :] = 0; u_proxima[:, 0] = 0; u_proxima[:, -1] = 0
        
Opción 2: Neumann (flujo de calor nulo en los bordes); la derivada normal de la temperatura es 0, no hay flujo a traves del borde, es decir, no entra ni sale calor por los bordes. 

        u_proxima[0, :] = u_proxima[1, :]
        u_proxima[-1, :] = u_proxima[-2, :]
        u_proxima[:, 0] = u_proxima[:, 1]
        u_proxima[:, -1] = u_proxima[:, -2]
        
Opción 3: Robin (convección en los bordes); es un caso intermedio entre Dirichlet y Neumann, representa la transferencia de calor por conveccion entre el sistema y un medio exterior.

        beta = 3.0 # Se puede modificar el valor  (Coeficiente de transferencia de calor en fronteras [W/m²K])
        u_proxima[0, :] = u_proxima[1, :] / (1 + beta * dx)
        u_proxima[-1, :] = u_proxima[-2, :] / (1 + beta * dx)
        u_proxima[:, 0] = u_proxima[:, 1] / (1 + beta * dy)
        u_proxima[:, -1] = u_proxima[:, -2] / (1 + beta * dy)
    
    
Actualizamos la solución completa para el siguiente paso de tiempo; se actualiza u para el siguiente paso y se guarda el estado si corresponde. 

        u[:, :] = u_proxima[:, :]  # Se copia el resultado actualizado de u_proxima a u
    
        # Guardar estado actual periodicamente
        if n % intervalo_guardado == 0:
            historial_temperaturas.append(u.copy())
    
Visualización; muestra la primera distribucion de temperatura en un grafico de contornos de color. 
    
    print(f"Total de estados guardados: {len(historial_temperaturas)}")
    
    fig, ax = plt.subplots(figsize=(7, 6))
    mapa_calor = ax.contourf(X, Y, historial_temperaturas[0], 20, cmap='hot')
    cbar = fig.colorbar(mapa_calor, ax=ax) 
    ax.set_title(f"Crank-Nicolson ADI 2D (t={0:.3f} s)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
Función para actualizar la animación; pdate(frame) es la funcion que actualiza el grafico en cada paso y FuncAnimation crea la animacion completa de la evolucion temporal. 

    def update(frame):
        ax.clear() # Limpiar el gráfico anterior
        mapa_calor = ax.contourf(X, Y, historial_temperaturas[frame], 20, cmap='hot') # Crear nuevo mapa de calor
        ax.set_title(f"t = {frame * intervalo_guardado * dt:.3f} s")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    
    anim = FuncAnimation(fig, update, frames=len(historial_temperaturas), interval=400, blit=False) # 400ms entre frames
    plt.close(fig) # Cerrar la figura estática

Guardar animacion; guarda la animacion en formato .mp4 

    mp4_filename = "difusion_calor_Onda-Robin.mp4"
    
    print(f"Saving animation to {mp4_filename}...")
    try:
        anim.save(mp4_filename, writer="ffmpeg", fps=10)
        print("Animation saved successfully.")
    except Exception as e:
        print(f"Error saving animation: {e}")
    
    
Mapa de calor 2D del resultado final, es decir, la visualizacion final de la distribucion de temperatura en t = T. 
    
    plt.figure(figsize=(7, 6))
    cp_final = plt.contourf(X, Y, u, 20, cmap='hot')
    plt.colorbar(cp_final)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.show()

  

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

















