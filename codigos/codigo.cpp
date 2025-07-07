#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>

//Función que crea una matriz llena de ceros, y la almacena en 1 vector de 1D.
std::vector<double> crearMatrizCeros(int n, int m) {
    return std::vector<double>(n * m, 0.0);
}

// Función que imprime una matriz creada anteriormente.
void imprimirMatriz(const std::vector<double>& matriz, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%.4f ", matriz[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");
}

//Función que multiplica la matriz anterior por un vector.
void multiplicarMatrizVector(const std::vector<double>& matriz, std::vector<double>& vec, int n) {
    std::vector<double> temp(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            temp[i] += matriz[i * n + j] * vec[j];
    vec = temp;
}

// Función que resuelve el sistema tridiagonal usando el método de Thomas.
void resolverTridiagonal(const std::vector<double>& T, std::vector<double>& u, int n) {
    std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), x(n, 0.0);
    for (int i = 0; i < n; i++) b[i] = T[i * n + i];
    for (int i = 0; i < n - 1; i++) {
        a[i + 1] = T[(i + 1) * n + i];
        c[i] = T[i * n + (i + 1)];
    }
    c[0] /= b[0];
    x[0] = u[0] / b[0];
    for (int i = 1; i < n; i++) {
        double m = b[i] - a[i] * c[i - 1];
        c[i] /= m;
        x[i] = (u[i] - a[i] * x[i - 1]) / m;
    }
    for (int i = n - 2; i >= 0; i--)
        x[i] -= c[i] * x[i + 1];
    u = x;
}

//Función que crea una matriz tridiagonal implicita.
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

// Función que crea una matriz tridiagonal explicita.
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

// Función que calcula el método de Crank Nicholson con Alternating Direction Implicit (ADI)
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
                fila[j] = matriz[i * n + j];
            multiplicarMatrizVector(S, fila, n);
            for (int j = 0; j < n; ++j)
                matriz[i * n + j] = fila[j];
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
                columna[i] = matriz[i * n + j];
            multiplicarMatrizVector(S, columna, n);
            for (int i = 0; i < n; ++i)
                matriz[i * n + j] = columna[i];
        }

        # pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            std::vector<double> fila(n);
            for (int j = 0; j < n; ++j)
                fila[j] = matriz[i * n + j];
            resolverTridiagonal(T, fila, n);
            for (int j = 0; j < n; ++j)
                matriz[i * n + j] = fila[j];
        }
    }
}

// Clase que se encarga del proceso de la solución de la ecuación de calor.
class solucion_ecuacion_calor {
private:
    int n;
    int pasos;
    double dt, t;
    double ds;
    double r;
    double alpha2;
    int opcion;
    int condicion_frontera;
    std::vector<double> matriz;

    double evaluar_condicion_ini(double x, double y) {
        if (opcion == 1) {
          //Condición inicial de pulso gaussiano centrado.
            double dx = x - 0.5;
            double dy = y - 0.5;
            return exp(-100 * (dx * dx + dy * dy));
        } else if (opcion == 2) {
          // Condición inicial de paraboloide centrado
            double dx = x - 0.5;
            double dy = y - 0.5;
            return 10.0 * (dx * dx + dy * dy);
        } else if (opcion == 3) {
          // Condición inicial de la onda senoidal suave
            const double Pi = 3.1415926535;
            return sin(2 * Pi * x) * sin(2 * Pi * y);
        }
        return 0.0;
    }

public:
    //Se crea un constructor:
    solucion_ecuacion_calor(int n_variable, double dt_variable, double t_variable, double alpha2_variable = 1.0) {
        n = n_variable;
        dt = dt_variable;
        t = t_variable;
        alpha2 = alpha2_variable;
        ds = 1.0 / (n - 1);
        r = alpha2 * dt / (ds * ds);
        pasos = static_cast<int>(t / dt);
        matriz = crearMatrizCeros(n, n);
    }

    // Se guarda la opcion escogida para la condición de frontera en la variable condicion_frontera.
    void condicion_front_escog(int condicion) {
        condicion_frontera = condicion;
    }

    // Se guarda la opcion escogida para la condición inicial en la variable opcion
    void opcion_escogida(int opcion_variable) {
        opcion = opcion_variable;
    }

    // Se aplica la condición inicial escogida
    void aplicar_condicion_inicial() {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double x = i * ds;
                double y = j * ds;
                matriz[j * n + i] = evaluar_condicion_ini(x, y);
            }
        }
    }

    // Se evalúa la condición de frontera correspondiente
    void evaluar_condicion_frontera() {
        double beta = 3.0;
        for (int i = 0; i < n; ++i) {
            if (condicion_frontera == 1) {
                matriz[i * n + 0] = 0;
                matriz[i * n + (n - 1)] = 0;
                matriz[0 * n + i] = 0;
                matriz[(n - 1) * n + i] = 0;
            } else if (condicion_frontera == 2) {
                matriz[i * n + 0] = matriz[i * n + 1];
                matriz[i * n + (n - 1)] = matriz[i * n + (n - 2)];
                matriz[0 * n + i] = matriz[1 * n + i];
                matriz[(n - 1) * n + i] = matriz[(n - 2) * n + i];
            } else if (condicion_frontera == 3) {
                matriz[i * n + 0] = matriz[i * n + 1] / (1 + beta * ds);
                matriz[i * n + (n - 1)] = matriz[i * n + (n - 2)] / (1 + beta * ds);
                matriz[0 * n + i] = matriz[1 * n + i] / (1 + beta * ds);
                matriz[(n - 1) * n + i] = matriz[(n - 2) * n + i] / (1 + beta * ds);
            }
        }
    }

    void resolver() {
        for (int paso = 0; paso < pasos; ++paso) {
            evaluar_condicion_frontera();
            CN_2D_ADI_Advance(matriz, n, r, 1, 0, 0, 0, 0);
        }
    }

    void imprimir() {
        imprimirMatriz(matriz, n);
    }
};

// Función main.
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
