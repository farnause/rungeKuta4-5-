
# Objectivo
El objetivo de este programa es implementar el método Runge-Kutta de orden 4 y 5 de 6 etapas con paso variable para resolver una ecuación diferencial con las condiciones iniciales dadas. Los coeficientes son

| $i$ | $c_i$ | $a_{i1}$ | $a_{i2}$ | $a_{i3}$ | $a_{i4}$ | $a_{i5}$ | $a_{i6}$ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | $0$ | $0$ | $0$ | $0$ | $0$ | $0$ | $0$ |
| 2 | $\frac14$ | $\frac14$ | $0$ | $0$ | $0$ | $0$ | $0$ |
| 3 | $\frac38$ | $\frac{3}{32}$ | $\frac{9}{32}$ | $0$ | $0$ | $0$ | $0$ |
| 4 | $\frac{12}{13}$ | $\frac{1932}{2197}$ | $-\frac{7200}{2197}$ | $\frac{7296}{2197}$ | $0$ | $0$ | $0$ |
| 5 | $1$ | $\frac{439}{216}$ | $-8$ | $\frac{3680}{513}$ | $-\frac{845}{4104}$ | $0$ | $0$ |
| 6 | $\frac12$ | $-\frac{8}{27}$ | $2$ | $-\frac{3544}{2565}$ | $\frac{1859}{4104}$ | $-\frac{11}{40}$ | $0$ |

Coeficientes de las etapas de cada orden

| $i$ | $b_i$ (orden 4) | $\hat b_i$ (orden 5) |
|---:|---:|---:|
| 1 | $\frac{25}{216}$ | $\frac{16}{135}$ |
| 2 | $0$ | $0$ |
| 3 | $\frac{1408}{2565}$ | $\frac{6656}{12825}$ |
| 4 | $\frac{2197}{4104}$ | $\frac{28561}{56430}$ |
| 5 | $-\frac15$ | $-\frac9{50}$ |
| 6 | $0$ | $\frac2{55}$ |

La gracia de un Runge-Kutta de orden 4 y 5 es que se puede estimar el error cometido en cada paso comparando las soluciones de orden 4 y 5:

$$
e_k \approx \left\|x_{k}^{(5)}-x_{k}^{(4)}\right\|.
$$

La EDO que se quiere resolver es

$$
\begin{cases}
\dot{x}=ax-y,\\
\dot{y}=x+ay,
\end{cases}
\qquad
x(0)=x_0,\quad y(0)=y_0.
$$

La solución exacta es

$$
\begin{aligned}
x(t)&=e^{at}\bigl(x_0\cos t-y_0\sin t\bigr),\\
y(t)&=e^{at}\bigl(x_0\sin t+y_0\cos t\bigr).
\end{aligned}
$$


# Cómo funciona
El programa cuenta con 4 archivos .c y sus respectivos archivos con los headers. En `Funcio.c` está declarado el campo de la EDO, en `rkf45.c` están las funciones necesarias para implementar el método Runge-Kutta 4(5), en `otros.c` hay funciones auxiliares y, por último, en `main.c` están declaradas las constantes relevantes para poder iniciar el Runge-Kutta. Concretamente, en este archivo se especifica
- **Condiciones iniciales**: El valor de la solución de la EDO en `t0`, `(x0,y0)`.
- **`n`**: dimensión del campo.
- **`h`**: variable que controla la longitud del paso actual dentro de `rk45`. 
- **Rango permitido de `h`**: donde `h` es el valor de la longitud del paso en cada iteración del Runge-Kutta. Concretamente, se define `hmin=tol` y `hmax=0.1`, que son el valor mínimo y máximo de `h` respectivamente.
- **`t` y `T`**: La variable `t` es la que controla el tiempo que hemos integrado, y su valor inicial es 0. Por otro lado, `T` corresponde al `t` final, o hasta qué tiempo queremos integrar. 
- **`a`**: parámetros que puede tener el campo de la EDO. Está declarado como un double con valor 0, pero sería sencillo adaptarlo a un vector que pudiera contener más de un valor. 
- **Tolerancia**: error máximo permitido del método Runge-Kutta, que corresponde a la variable `tol`. Por defecto su valor es `tol=1e-10`.
- **`e` y `ev`**: `e` corresponde al módulo en norma 1 del error estimado total. Por otro lado, `ev` es un vector, en que cada coordenada representa el error estimado cometido en la coordenada correspondiente de la solución. 
- **`npasmax`**: número de pasos que le permitimos hacer al programa.

La función `flux` es la que hace la integración numérica hasta el tiempo final especificado, utilizando `rk45` para hacer cada paso de la integración. La función `flux` se puede parar en caso de encontrar algún error, como que el número de pasos hecho sea mayor que `npasmax` o en caso de que `rk45` encuentre algún problema. 

# Descripción
### Funciones presentes en `Funcio.c`
- **camp**: función que hace el papel del campo de la EDO. 

Argumentos: `int n`: dimensión del campo,`double t`: variable que controla el tiempo de integración, `double *x`: pointer que apunta al vector donde se va almacenando el valor de la solución de la EDO, `double *f`: pointer que apunta a un vector que almacena el valor del campo después de recibir el valor de la solución en el tiempo correspondiente, `double *prm`: parámetros que pueda tener el campo.


### Funciones presentes en `otros.c`
- **Abs**: función que te devuelve el valor absoluto del input. Debido a experiencias pasadas, prefiero declarar yo mismo la función valor absoluto y no usar la de `math.h`.  

Argumentos: `double x`: valor del que se quiere el valor absoluto.  

- **minomax**: función que, dependiendo del input, te devuelve el mínimo o máximo entre dos valores.  

Argumentos: `double x, double y`: valores que queremos comparar, `int a`: si `a=-1` la función devuelve el mínimo, si `a=1` te devuelve el máximo.  

- **Start**: función void que acepta un vector y lo convierte en el vector nulo.

Argumentos: `double *v`: vector que queremos convertir en el vector nulo, `int n`: dimensión del vector.

### Funciones presentes en `rkf45.c`
- **rk45**: función que efectúa cada paso de la integración numérica, controlando el valor de `h` para forzar que el error se mantenga por debajo de la tolerancia. 

Argumentos: `double *ev`: vector que almacenará el error estimado en cada coordenada de la solución, `double *e`: error en norma 1 del error estimado, `double *t`: variable que controla internamente el tiempo integrado, `double *x`: vector en el que se almacenará la solución,  `double *h`: variable que controla internamente la longitud del paso de la iteración actual, `double hmin` i `double hmax`: rango de valores permitidos para `h`, `double tol`: tolerancia permitida, `int n`: dimensión del campo de la EDO, `int (*camp)(...)`: campo de la EDO, `double *prm`: parámetros del campo de la EDO.

- **flux**: función que integra hacia el tiempo final indicado, o hasta que se superen el número máximo de pasos indicado. Según el signo de `T`, integra hacia atrás o hacia adelante en el tiempo, respecto el valor de `t` inicial. 

Argumentos: hereda los argumentos de `rk45` y añade `int npasmax`: número de pasos máximos que permitimos. 

## Observaciones
Es posible mostrar el error estimado utilizando `e` y `ev`, aunque no está programado que se muestre por defecto.
# Comentarios personales respecto el programa
Esta fue una de las prácticas que hice en la optativa de cálculo numérico de 4º curso del grado de Matemáticas en la UAB. Esta práctica me hizo sentir mucho más cómodo trabajando con memoria dinámica. 
También me obligó a organizar mis programas en varios archivos, cosa que permite no tener un main extremadamente largo. La organización de mis códigos siempre sufre de este problema y esta práctica fue un paso hacia mejorar este aspecto en particular.    
Como mejora del programa, se podría dar más libertad al usuario para introducir las condiciones iniciales o incluso la tolerancia.  
