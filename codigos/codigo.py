#Este codigo incluye 3 opciones de condiciones de frontera y 3 opciones de condiciones iniciales, por lo que para graficar, solo se debe de dejar una de ambas condiciones.  

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from matplotlib.animation import FuncAnimation
from IPython.display import Image

# Parámetros
alpha = 1.0 #Constante de difusión térmica

#Configuración del dominio
Longitud_x, Longitud_y = 1.0, 1.0  # Estos son los largos del dominio en x e y respectivamente 
num_divisiones_x, num_divisiones_y = 20, 20  # Cantidad de divisiones en las direcciones x e y

#Tamaño de cada celda (Resolución espacial)
dx, dy = Longitud_x / num_divisiones_x, Longitud_y / num_divisiones_y  # Distancia entre nodos consecutivos en cada eje

#Malla de coordenadas
x = np.linspace(0, Longitud_x, num_divisiones_x+1)  # Puntos equiespaciados de 0 a Lx (incluyendo extremos)
y = np.linspace(0, Longitud_y, num_divisiones_y+1)  # Puntos equiespaciados de 0 a Ly (incluyendo extremos)
X, Y = np.meshgrid(x, y, indexing='ij')

#Configuración temporal
T = 0.2      # Tiempo total de simulación
dt = 0.001   # Paso temporal
nt = int(T / dt)  # Número total de pasos de tiempo

# Cálculo de parámetros auxiliares para Crank-Nicolson en x e y
parametro_x = alpha * dt / (2 * dx**2)  # Parámetro de difusión en x, dividido entre 2 por Crank-Nicolson
parametro_y = alpha * dt / (2 * dy**2)  # Igual pero en y

# Construcción de matrices tridiagonales para el método implícito en x
Ax = diags([[-parametro_x]*(num_divisiones_x-1), [1+2*parametro_x]*(num_divisiones_x-1), [-parametro_x]*(num_divisiones_x-1)], [-1,0,1], shape=(num_divisiones_x-1, num_divisiones_x-1))  # Matriz del lado izquierdo (implícito) en x
Bx = diags([[parametro_x]*(num_divisiones_x-1), [1-2*parametro_x]*(num_divisiones_x-1), [parametro_x]*(num_divisiones_x-1)], [-1,0,1], shape=(num_divisiones_x-1, num_divisiones_x-1))    # Matriz del lado derecho (explícito) en x

# Construcción de matrices tridiagonales para el método implícito en y
Ay = diags([[-parametro_y]*(num_divisiones_y-1), [1+2*parametro_y]*(num_divisiones_y-1), [-parametro_y]*(num_divisiones_y-1)], [-1,0,1], shape=(num_divisiones_y-1, num_divisiones_y-1))  # Implícito en y
By = diags([[parametro_y]*(num_divisiones_y-1), [1-2*parametro_y]*(num_divisiones_y-1), [parametro_y]*(num_divisiones_y-1)], [-1,0,1], shape=(num_divisiones_y-1, num_divisiones_y-1))    # Explícito en y

# Condición Inicial
u = np.zeros((num_divisiones_x+1, num_divisiones_y+1))  # Inicialización de la matriz de temperatura u en todo el dominio

# Opción 1: Pulso gaussiano centrado en (0.5, 0.5)
#u[:, :] = np.exp(-100 * ((X - 0.5)**2 + (Y - 0.5)**2))

# Opción 2: Paraboloide centrado (Alternativa físicamente consistente)
#u[:, :] = 10 * ((X - 0.5)**2 + (Y - 0.5)**2)

# Opción 3: Onda senosoidal suave (Para patrones periódicos)
u[:, :] = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)

# Matriz para almacenar la próxima temperatura
u_proxima = np.zeros_like(u) # Matriz que contendrá la temperatura actualizada


# SIMULACIÓN

# Lista para guardar estados intermedios
historial_temperaturas = []
intervalo_guardado = 10 # Guarda cada 10 pasos de tiempo

# Comienza el lazo de integración temporal
for n in range(nt):
    # Paso intermedio: resolvemos en la dirección x, manteniendo y fijo
    u_intermedia = np.zeros_like(u)  # Matriz temporal para almacenar resultados intermedios
    # Recorremos cada fila fija (j) y resolvemos en x (columnas)
    for j in range(1, num_divisiones_y):
        rhs = Bx.dot(u[1:num_divisiones_x, j])              # Lado derecho del sistema: combinación explícita
        u_intermedia[1:num_divisiones_x, j] = spsolve(Ax, rhs)    # Resolvemos el sistema lineal en x

    # Paso final: resolvemos en la dirección y, manteniendo x fijo
    for i in range(1, num_divisiones_x):
        rhs = By.dot(u_intermedia[i, 1:num_divisiones_y])  # Lado derecho para dirección y
        u_proxima[i, 1:num_divisiones_y] = spsolve(Ay, rhs)   # Solución del sistema lineal en y

    # Aplicar la Condición de Frontera seleccionada a u_proxima 
    
    # Opción 1: Dirichlet (bordes fijos en 0)
    #u_proxima[0, :] = 0; u_proxima[-1, :] = 0; u_proxima[:, 0] = 0; u_proxima[:, -1] = 0
    
    # Opción 2: Neumann (flujo de calor nulo en los bordes) 
    #u_proxima[0, :] = u_proxima[1, :]
    #u_proxima[-1, :] = u_proxima[-2, :]
    #u_proxima[:, 0] = u_proxima[:, 1]
    #u_proxima[:, -1] = u_proxima[:, -2]
    
    #Opción 3: Robin (convección en los bordes) 
    #beta = 3.0 # Se puede modificar el valor  (Coeficiente de transferencia de calor en fronteras [W/m²K])
    #u_proxima[0, :] = u_proxima[1, :] / (1 + beta * dx)
    #u_proxima[-1, :] = u_proxima[-2, :] / (1 + beta * dx)
    #u_proxima[:, 0] = u_proxima[:, 1] / (1 + beta * dy)
    #u_proxima[:, -1] = u_proxima[:, -2] / (1 + beta * dy)


    # Actualizamos la solución completa para el siguiente paso de tiempo
    u[:, :] = u_proxima[:, :]  # Se copia el resultado actualizado de u_proxima a u

    # Guardar estado actual periodicamente
    if n % intervalo_guardado == 0:
        historial_temperaturas.append(u.copy())

# Visualización

print(f"Total de estados guardados: {len(historial_temperaturas)}")

fig, ax = plt.subplots(figsize=(7, 6))
mapa_calor = ax.contourf(X, Y, historial_temperaturas[0], 20, cmap='hot')
cbar = fig.colorbar(mapa_calor, ax=ax) 
ax.set_title(f"Crank-Nicolson ADI 2D (t={0:.3f} s)")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Función para actualizar la animación
def update(frame):
    ax.clear() # Limpiar el gráfico anterior
    mapa_calor = ax.contourf(X, Y, historial_temperaturas[frame], 20, cmap='hot') # Crear nuevo mapa de calor
    ax.set_title(f"t = {frame * intervalo_guardado * dt:.3f} s")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

  # No es necesario retornar objetos cuando blit=False
    pass 

# Crear animación

anim = FuncAnimation(fig, update, frames=len(historial_temperaturas), interval=400, blit=False) # 400ms entre frames
plt.close(fig) # Cerrar la figura estática

# Para guardar 
mp4_filename = "difusion_calor_Onda-Robin.mp4"

print(f"Saving animation to {mp4_filename}...")
try:
    anim.save(mp4_filename, writer="ffmpeg", fps=10)
    print("Animation saved successfully.")
except Exception as e:
    print(f"Error saving animation: {e}")


# Mapa de calor 2D del resultado final

plt.figure(figsize=(7, 6))
cp_final = plt.contourf(X, Y, u, 20, cmap='hot')
plt.colorbar(cp_final)
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()
