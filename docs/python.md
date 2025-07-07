# Soluión en Python 

En el presente apartado explicaremos a detalle cada linea de codigo de la resolucion de la Ecuacion de calor realizada en Python.

Iniciamos importando las librerias, en donde numpy nos ayuda a realizar calculos numericos y el manejo de arreglos/matirces; matplotlib.pyplot, grafica y visualiza datos; 
scipy.sparse.diags, permite crear matrices dispersas que ahorran memoria y son mas eficientes; scipy.sparse.linalg.spsolve, resuelve sistemas lineales que usan matrices dispersas; FuncAnimation, permite crear animaciones y Image, sirve para mostrar imagenes en entornos como Jupyter Notebook. 

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    from matplotlib.animation import FuncAnimation
    from IPython.display import Image

Se define la clase principal para simular la ecuacion de calor en 2D. En donde el constructor de la clase es __init__ y alpha es la difusividad termica, Lx, Ly son las dimensiones del dominio, Nx, Ny son las divisiones del espacio y dt, T es el paso y tiempo total de simulacion. 

    class SimuladorDifusionCalor2D:
        
        def __init__(self, alpha, Lx, Ly, Nx, Ny, dt, T):

            # Parámetros físicos y de discretización espacial
            self.alpha = alpha # Coeficiente de difusión térmica
            self.Lx, self.Ly = Lx, Ly # Dimensiones del dominio (largo y ancho)
            self.Nx, self.Ny = Nx, Ny # Cantidad de divisiones en las direcciones x e y
            self.dx = Lx / Nx # Tamaño del paso espacial en x
            self.dy = Ly / Ny # Tamaño del paso espacial en y

Define las propiedades fisicas y espaciales del dominio. 
            
            self.x = np.linspace(0, Lx, Nx + 1) # Crear un array con las coordenadas x de la malla
            self.y = np.linspace(0, Ly, Ny + 1) # Crear un array con las coordenadas y de la malla
            # malla de coordenadas 2D
            self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')

Crea la malla espacial 2D.

            # Parámetros temporales
            self.dt = dt # Tamaño del paso de tiempo
            self.T = T # Tiempo total de simulación
            self.nt = int(T / dt) # Número total de pasos de tiempo

Configura la discretizacion temporal.
    
            # Parámetros usados en el método ADI 
            self.parametro_x = alpha * dt / (2 * self.dx ** 2) # Parámetro de difusión en dirección x
            self.parametro_y = alpha * dt / (2 * self.dy ** 2) # Parámetro de difusión en dirección y

Coeficientes del metodo ADI.
    
            # Construir las matrices dispersas necesarias para el método ADI
            self._construir_matrices()

Llama al metodo para construir las matrices tridiagonales. En donde se inicializa las matrices de temperatura y el historial para animacion.             
    
            # Inicializar las matrices de temperatura
            self.u = np.zeros((Nx + 1, Ny + 1)) # Matriz de temperatura actual, inicializada a ceros
            self.u_proxima = np.zeros_like(self.u) # Matriz para la temperatura en el siguiente paso de tiempo, inicializada de forma similar a u
            self.historial_temperaturas = [] # Lista para almacenar el estado de la temperatura en diferentes momentos
            self.intervalo_guardado = 10  # Guarda el estado de la temperatura cada 10 pasos de tiempo

Construccion de matrices tridiagonales
    
        # Método para construir las matrices dispersas 
        def _construir_matrices(self):
        
            Nx, Ny = self.Nx, self.Ny # Obtener el número de puntos de discretización
            px, py = self.parametro_x, self.parametro_y # Obtener los parámetros de difusión

EL self.Ax cre matrices implicitas y el self.Bx crea mastrices explicitas para la direccion x. 

            # Construir la matriz Ax para la difusión en x (implícita)
            # Es una matriz tridiagonal (diagonal principal, subdiagonal y superdiagonal)
            self.Ax = diags([[-px] * (Nx - 2), [1 + 2 * px] * (Nx - 1), [-px] * (Nx - 2)], [-1, 0, 1], shape=(Nx-1, Nx-1))
            # Construir la matriz Bx para la difusión en x (parte explícita)
            # Es una matriz tridiagonal
            self.Bx = diags([[px] * (Nx - 2), [1 - 2 * px] * (Nx - 1), [px] * (Nx - 2)], [-1, 0, 1], shape=(Nx-1, Nx-1))

EL self.Ay cre matrices implicitas y el self.By crea mastrices explicitas para la direccion y. 
    
            # Construir la matriz Ay para la difusión en y (implícita)
            self.Ay = diags([[-py] * (Ny - 2), [1 + 2 * py] * (Ny - 1), [-py] * (Ny - 2)], [-1, 0, 1], shape=(Ny-1, Ny-1))
            # Construir la matriz By para la difusión en y (explícita)
            self.By = diags([[py] * (Ny - 2), [1 - 2 * py] * (Ny - 1), [py] * (Ny - 2)], [-1, 0, 1], shape=(Ny-1, Ny-1))

Condiciones iniciales; permite elegir entre el pulso gaussiano centrado, paraboloide invertido y la ondsa senoidal. 
    
        # Establecer la condición inicial de temperatura
        def condicion_inicial(self, tipo):
            
            # Pulso gaussiano centrado en (0.5, 0.5)
            if tipo == 'gauss':
                self.u[:, :] = np.exp(-100 * ((self.X - 0.5) ** 2 + (self.Y - 0.5) ** 2))
            # Paraboloide centrado
            elif tipo == 'paraboloide':
                self.u[:, :] = 10 * ((self.X - 0.5) ** 2 + (self.Y - 0.5) ** 2)
            # Onda Senosoidal
            elif tipo == 'senoidal':
                self.u[:, :] = np.sin(2 * np.pi * self.X) * np.sin(2 * np.pi * self.Y)
            # Si no funciona, lanzar error de condición incial no conocida
            else:
                raise ValueError(f"Condición inicial desconocida: {tipo}")

Simulacion; ejecuta la simulacion usando ADI, en donde el primer medio paso es implicito en x, el segundo medio paso es implicito en y, y por ultimo, se aplican las condiciones de frontera. 
    
        # Método para ejecutar la simulación
        def simular(self, frontera):
            # Iterar sobre el número total de pasos de tiempo
            for n in range(self.nt):
                # Matriz temporal para almacenar resultados intermedios
                u_intermedia = np.zeros_like(self.u)
                
                # Recorremos cada fila fija (j) y resolvemos en x (columnas)
                for j in range(1, self.Ny):
                    # Calcular el lado derecho del sistema lineal en la dirección x
                    rhs = self.Bx.dot(self.u[1:self.Nx, j])
                    # Resolver el sistema lineal en x
                    u_intermedia[1:self.Nx, j] = spsolve(self.Ax, rhs)
    
                # Segundo medio paso (explícito en x, implícito en y)
                # resolvemos en la dirección y, manteniendo x fijo
                for i in range(1, self.Nx):
                    # Calcular el lado derecho del sistema lineal en la dirección y
                    rhs = self.By.dot(u_intermedia[i, 1:self.Ny])
                    # Resolver el sistema lineal en y
                    self.u_proxima[i, 1:self.Ny] = spsolve(self.Ay, rhs)
    
                # Aplicar condiciones de frontera 
                if frontera == 'dirichlet':
                    # Condición de frontera de Dirichlet (temperatura fija en los bordes)
                    self.u_proxima[0, :] = 0 # Borde izquierdo
                    self.u_proxima[-1, :] = 0 # Borde derecho
                    self.u_proxima[:, 0] = 0 # Borde inferior
                    self.u_proxima[:, -1] = 0 # Borde superior
                elif frontera == 'neumann':
                    # Condición de frontera de Neumann (flujo de calor cero en los bordes)
                    self.u_proxima[0, :] = self.u_proxima[1, :] # Borde izquierdo
                    self.u_proxima[-1, :] = self.u_proxima[-2, :] # Borde derecho
                    self.u_proxima[:, 0] = self.u_proxima[:, 1] # Borde inferior
                    self.u_proxima[:, -1] = self.u_proxima[:, -2] # Borde superior
                elif frontera == 'robin':
                    # Condición de frontera de Robin 
                    beta = 3.0 # Coeficiente de convección 
                    self.u_proxima[0, :] = self.u_proxima[1, :] / (1 + beta * self.dx) # Borde izquierdo
                    self.u_proxima[-1, :] = self.u_proxima[-2, :] / (1 + beta * self.dx) # Borde derecho
                    self.u_proxima[:, 0] = self.u_proxima[:, 1] / (1 + beta * self.dy) # Borde inferior
                    self.u_proxima[:, -1] = self.u_proxima[:, -2] / (1 + beta * self.dy) # Borde superior
                # Lanzar un error si el tipo de frontera no es válido
                else:
                    raise ValueError(f"Condición de frontera no válida: {frontera}")
    
Actualiza la matriz de temperatura para el siguiente paso de tiempo.

                self.u[:, :] = self.u_proxima[:, :]
    
Guarda el estado actual periodicamente.

                if n % self.intervalo_guardado == 0:
                    self.historial_temperaturas.append(self.u.copy()) #.copy() para guardar una copia y no una referencia
    
Muestra el mapa de calor en el estado final de la simulación.

        def mostrar_estado_final(self):
            plt.figure(figsize=(7, 6))
            mapa = plt.contourf(self.X, self.Y, self.historial_temperaturas[-1], 20, cmap='hot') # Crea un gráfico de mapa de calor del estado final de la temperatura
            plt.colorbar(mapa)  # Agrega barra de color para indicar los valores de temperatura
            plt.title(f"Mapa de Calor")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.tight_layout()
            plt.show()

Crea una animacion de la evolucion de la temperatura. 
    
        def animar(self, intervalo_ms=400, fps=10):
    
            fig, ax = plt.subplots(figsize=(7, 6))
            cont = ax.contourf(self.X, self.Y, self.historial_temperaturas[0], 20, cmap='hot')
            fig.colorbar(cont, ax=ax)
    
Función que se llama en cada frame de la animación.

            def update(frame):
                # Limpiar los ejes para dibujar el nuevo frame
                ax.clear()
                # Crear el gráfico de contorno para el estado de temperatura del frame actual
                cont = ax.contourf(self.X, self.Y, self.historial_temperaturas[frame], 20, cmap='hot')
                # Actualizar el título con el tiempo actual de la simulación
                ax.set_title(f"t = {frame * self.intervalo_guardado * self.dt:.3f} s")
                ax.set_xlabel("x")
                ax.set_ylabel("y")
    
Crea el objeto FuncAnimation.

            anim = FuncAnimation(fig, update, frames=len(self.historial_temperaturas), interval=intervalo_ms, blit=False)
            # Cerrar la figura para evitar que se muestre como un gráfico estático
            plt.close(fig)
            # Mostrar la animación como un video HTML 
            display(HTML(anim.to_jshtml()))
    
Llama la clase simulador y establecer las condiciones iniciales y de frontera. En donde sim crea el objeto de simulacion, sim.condicion_inicial('paraboloide')  establece la condicion inicial tipo paraboloide y sim.simular(frontera='neumann') corre la simulacion con condiciones de Neumann.

    # Caso: Condicion inicial de Paraboloide y condición de frontera de Neumann
    sim = SimuladorDifusionCalor2D(alpha=1.0, Lx=1.0, Ly=1.0, Nx=40, Ny=40, dt=0.0005, T=0.1)
    sim.condicion_inicial('paraboloide')  
    sim.simular(frontera='neumann') 
    
Mapa de calor para condición inicial de Paraboloide y condición de frontera de Neumann; muestra el mapa de calor de la ultima distribucion de la temperatura. 

    sim.mostrar_estado_final()
    
Animación para condición inicial de Paraboloide y condición de frontera de Neumann.

    sim.animar()
    
Caso: Condicion inicial de Paraboloide y condición de frontera de Neumann.

    sim_gaussiano = SimuladorDifusionCalor2D(alpha=1.0, Lx=1.0, Ly=1.0, Nx=40, Ny=40, dt=0.0005, T=0.1)
    sim_gaussiano.condicion_inicial('gauss')  
    sim_gaussiano.simular(frontera='dirichlet')  
    
Mapa de calor para condición inicial de Gauss y condición de frontera de Dirichlet.

    sim_gaussiano.mostrar_estado_final()
    
Animación para condición inicial de Gauss y condición de frontera de Dirichlet.

    sim_gaussiano.animar()
    
Caso: Condicion inicial de onda senoidal y condición de frontera de Dirichlet.

    sim_senoidal = SimuladorDifusionCalor2D(alpha=1.0, Lx=1.0, Ly=1.0, Nx=40, Ny=40, dt=0.0005, T=0.1)
    sim_senoidal.condicion_inicial('senoidal')  
    sim_senoidal.simular(frontera='robin')  
    
Mapa de calor para condición inicial de  y condición de frontera de Dirichlet.

    sim_senoidal.mostrar_estado_final()
    
Animación para condición inicial de  y condición de frontera de Dirichlet.

    sim_senoidal.animar()
    
        
        
        
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

















