# Solución en C++
El presente código en C++ tiene como objetivo resolver la ecuación de calor en dos dimensiones. Para resolverla de forma eficiente, el código implementa el método Crank–Nicolson con Alternating Direction Implicit (ADI). Este método permite dividir cada paso temporal en dos subpasos alternados en dirección horizontal y vertical, facilitando la resolución de matrices tridiagonales mediante el eficiente método de Thomas.

Además, se integran técnicas de paralelismo de memoria compartida usando OpenMP, lo que permite reducir significativamente los tiempos de cálculo. El código está diseñado para permitir la elección entre distintas condiciones iniciales

---------------------------------------------------------

La siguiente función crea una matriz representada como un vector de 1D:\
Se utilizaron las siguientes bibliotecas:
 
* `<iostream>`: Para poder utilizar inputs y outputs con std.
* `<cmath>`: Para poder agregar funciones matemáticas
* `<vector>`: Para poder usar std::vector
* `<omp.h>`: Para paralelizar con OpenMP

```cpp
#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
```

**Función que crea una matriz llena de ceros, y la almacena en 1 vector de 1D**\
`n` Número de filas.\
`m` Número de columnas.\
`return std::vector<double>` Se retorna un vector `nxm` que contiene los arreglos de ceros

```cpp
std::vector<double> crearMatrizCeros(int n, int m) {
    return std::vector<double>(n * m, 0.0);
}
```

**Función que imprime la matriz creada anteriormente**\
Se imprime en la consola la matriz nxn contenida en un vector de 1D. Para esto se utilzan 3 cifras decimales por cada entrada de la matriz.\
`matriz` Vector que contiene la matriz creada previamente de forma lineal.\
`n` Número de filas de la matriz.

```cpp
void imprimirMatriz(const std::vector<double>& matriz, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%.3f ", matriz[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");
}
```

**Función que multiplica la matriz anteior por un vector.**\
Se hace una multiplicación entre la matriz guardada en el vector 1D y un vector de las mismas dimensiones. 
El resultado obtenido de la multiplicación sobreescribe el vector original utilizado.\
`matriz` Vector que contiene la matriz creada previamente de forma lineal.\
`vec` Vector por el cual se multiplica la matriz, y en el cual se guardan los resultados de la multiplicación.\
`n` Número de filas de la matriz.\

```cpp
void multiplicarMatrizVector(const std::vector<double>& matriz, std::vector<double>& vec, int n) {
    std::vector<double> temp(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            temp[i] += matriz[i * n + j] * vec[j];
    vec = temp; 
}
```

**Función que resuelve el sistema tridiagonal usando el método de Thomas**\
Se resuelve el sistema de ecuaciones lineales obtenido, el cual incluye una matriz tridiagonal, utilizando el método de Thomas. La solución se almacena en el vector u, sobrescribiendo así su contenido.\
`T` Matriz tridiagonal resuelta por medio del método de Thomas, y almacenada de forma 1D.\
`u` Vector con los términos independientes del sistema tridiagonal.\
`n` Dimensión de la matriz.\
El resultado final se guarda en el vector `u`.

```cpp
void resolverTridiagonal(const std::vector<double>& T, std::vector<double>& u, int n) {
    std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), x(n, 0.0);
    for (int i = 0; i < n; i++) b[i] = T[i * n + i];
    for (int i = 0; i < n-1; i++) {
        a[i+1] = T[(i+1) * n + i];
        c[i] = T[i * n + (i + 1)];
    }
    c[0] /= b[0];
    x[0] = u[0] / b[0];
    for (int i = 1; i < n; i++) {
        double m = b[i] - a[i] * c[i-1];
        c[i] /= m;
        x[i] = (u[i] - a[i] * x[i-1]) / m;
    }
    for (int i = n - 2; i >= 0; i--)
        x[i] -= c[i] * x[i + 1];
    u = x;
}
```

 **Función que crea una matriz tridiagonal implicita**\
 Se crea una matriz tridiagonal implícita con valores 1+2r en la diagonal principal y -r en las otras dos diagonales.\
 `n` Dimensión de la matriz.\
 `r` Valor de la variable discretizada que relaciona la parte espacial y temporal de la ecuación diferencial.\
`return Laplaciano` Se retorna un Laplaciano tridiagonal implícito nxn.

```cpp
std::vector<double> crearLaplacianoImplicito(int n, double r) {
    std::vector<double> Laplaciano(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (i == 0 || i == n - 1) {
            Laplaciano[i * n + i] = 1.0;
        } else {
            Laplaciano[i * n + (i - 1)] = -r;
            Laplaciano[i * n + i] = 1 + 2 * r;
            Laplaciano[i * n + (i + 1)] = -r;
        }
    }
    return Laplaciano;
}
```

**Función que crea una matriz tridiagonal explicita**\
Se crea una matriz tridiagonal explícita con valores 1-2r en la diagonal principal y r en las otras dos diagonales.\
`n` Dimensión de la matriz.\
`r` Valor de la variable discretizada que relaciona la parte espacial y temporal de la ecuación diferencial.\
`return Laplaciano` Se retorna un Laplaciano tridiagonal explícito nxn.

```cpp
std::vector<double> crearLaplacianoExplicito(int n, double r) {
    std::vector<double> Laplaciano(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (i == 0 || i == n - 1) {
            Laplaciano[i * n + i] = 1;
        } else {
            Laplaciano[i * n + (i - 1)] = r;
            Laplaciano[i * n + i] = 1 - 2 * r;
            Laplaciano[i * n + (i + 1)] = r;
        }
    }
    return  Laplaciano;
}
```

**Función que calcula el método de Crank Nicholson con Alternating Direction Implicit (ADI)**\
La función resuelve la ecuación de calor en 2D dividiendo la parte temporal de problema en 2 subpasos:
Una parte para la dirección y y otra para la dirección x. Esto facilita poder utilizar el método de Crank
Nicholson al reducirlo con matrices tridiagonales. Además en cada subpaso se calcula la multiplicación de de matrices y la resolución de la matriz tridiagonal por medio del método de Thomas.\
Para la parte de la paralelización se utilizó estratégicamente en los bucles for que involucran filas y columnas  de la malla, ya que al agregarlo en dichos `for`, la velocidad del código mejora considerablemente.\
`matriz` Matriz almacenada como un vector 1D.\
`n` Número de puntos a utilizar\
`r` Valor de la variable discretizada que relaciona la parte espacial y temporal de la ecuación diferencial.\
`pasos` Número de pasos a utilizar.\
`bordeIzq` Valor de frontera en el borde izquierdo.\
`bordeDer` Valor de frontera en el borde derecho.\
`bordeInf` Valor de frontera del borde inferior.\
`bordeSup` Valor de frontera del borde superior.

```cpp
void CN_2D_ADI_Advance(std::vector<double>& matriz, int n, double r, int pasos,
                       double bordeIzq, double bordeDer, double bordeInf, double bordeSup) {

    r = r / 2.0;
    std::vector<double> S = crearLaplacianoExplicito(n, r);
    std::vector<double> T = crearLaplacianoImplicito(n, r);
    std::vector<double> temp(n * n, 0.0);
    std::vector<double> fila(n);
    std::vector<double> columna(n);

    for (int t = 0; t < pasos; ++t) {

        # pragma omp parallel for
        for (int i = 0; i < n; ++i) {
          std::vector<double> fila(n);
          for (int j = 0; j < n; ++j)
            fila[j] = matriz[i * n +j];
          multiplicarMatrizVector(S, fila, n);
          for (int j = 0; j < n; ++j)
            matriz[i * n +j] = fila[j];
        }


        # pragma omp parallel for
        for (int j = 0; j < n; ++j) {
            std::vector<double> columna(n);
            for (int i = 0; i < n; ++i)
                columna[i] = matriz[i * n + j];
            resolverTridiagonal(T, columna, n);
            for (int i = 0; i < n; ++i)
                matriz[i * n + j] = columna[i];
        }



        # pragma omp parallel for
        for (int j = 0; j < n; ++j) {
          std::vector<double> columna(n);
          for (int i = 0; i < n; ++i)
            columna[i] = matriz[i * n +j];
          multiplicarMatrizVector(S, columna, n);
          for (int i = 0; i < n; ++i)
            matriz[i * n + j] = columna[i];
        }


   
        # pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            std::vector<double> fila(n);
            for (int j = 0; j < n; ++j)
                fila[j] = matriz[i * n +j];
            resolverTridiagonal(T, fila, n);
            for (int j = 0; j < n; ++j)
                matriz[i * n + j] = fila[j];
        }
    }
}
```

**Clase** 
`solucion_ecuacion_calor` Consiste en los procesos necesarios para resolver el problema.
En esta se configuran las condiciones iniciales y de frontera, además de que se crean las funciones con las cuáles se calcula la solución del 
problema.

**Atributos**
- `n, dt, t, alpha2, ds, r, pasos`: parámetros numéricos.
- `opcion`: condición inicial seleccionada.
- `condicion_frontera`: tipo de borde (1: Dirichlet, 2: Neumann, 3: Robin).
- `matriz`: contiene los valores de temperatura.

**Funciones públicas**

```cpp
solucion_ecuacion_calor(int n, double dt, double t, double alpha2 = 1.0);
void opcion_escogida(int opcion_variable);
void condicion_front_escog(int condicion);
void aplicar_condicion_inicial();
void evaluar_condicion_frontera();
void resolver();
void imprimir();
```

**Funciones privadas**

```cpp
double evaluar_condicion_ini(double x, double y);
```
Evaluá la condición inicial dependiendo de la opción escogida.

**Función main**\
Inicializa todo el código para resolver la ecuación de calor en 2-D. Para esto el usuario primero debe de escoger una de las tres condiciones iniciales y una de las tres condiciones de frontera. El programa evalúa dicha condición incial usando el método de Crank Nicholson y la reducción de la matriz tridiagonal con el método de Thomas, y luego imprime el resultado obtenido.\
`return 0` si el programa se ejecuta sin ningún problema, y 1 si la opción elegida no era válida.

```cpp
int main() {
    int ns = 50;
    double dt = 0.0005, t = 0.1;

    int opcion_cond_ini;
    std::cout << "Escoja una de las siguientes condiciones iniciales al ingresar el número correspondiente: \n";
    std::cout << " 1. Pulso Gaussiano centrado \n2. Paraboloide centrado \n3. Onda senoidal suave \n";
    std::cin >> opcion_cond_ini;

    int condicion_frontera;
    std::cout << "Escoja una de las siguientes condiciones de frontera al ingresar el número correspondiente: \n";
    std::cout << " 1. Dirichlet \n2. Neuman \n3. Robin \n";
    std::cin >> condicion_frontera;

    if (opcion_cond_ini < 1 || opcion_cond_ini > 3) {
        std::cerr << "La opción elegida no es válida. \n";
        return 1;
    }

    // Se crea el constructor y se usa para llamar las funciones.
    solucion_ecuacion_calor constructor(ns, dt, t);
    constructor.opcion_escogida(opcion_cond_ini);
    constructor.condicion_front_escog(condicion_frontera);
    constructor.aplicar_condicion_inicial();
    constructor.resolver();
    constructor.imprimir();

    return 0;
}
```
