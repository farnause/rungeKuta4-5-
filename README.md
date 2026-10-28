# Objectivo
El objetivo de este programa es implementar el método Runge-Kutta de orden 4 y 5 de 6 etapas con paso variable para resolver una equación diferencial con las condiciones iniciales dadas. Los coeficientes son

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

La gracia de un Runge-Kutta de orden 4 y 5 se puede estimar el error cometido en cada paso comparando las soluciones de orden 4 y 5:

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
El programa cuenta con 4 archivos .c y sus respectivos archivos con los headers. En `Funcio.c` esta declarado el campo de la EDO, en `rkf45.c` esta programado el métoddo Runge-Kutta 4(5), en `otros.c` hay funciones auxiliares y por último en `main.c` estan declaradas las constantes relevantes para poder iniciar el Runge-Kutta. Concretamente, en este archivo se especifica
- **Condiciones iniciales**: El valor de la solución de la EDO en `t0`, `(x0,y0)`.
- **`n`**: dimensión del campo.
- **`h`**: variable que controla la longitud del paso actual dentro de `rk45`. 
- **Rango permitido de `h`**: donde `h` es el valor de la longitud del paso en cada iteración del Runge-Kutta. Concretamente se define `hmin=tol` y `hmax=0.1`, que son el valor mínimo y máximo de `h` respectivamente.
- **`t` y `T`**: La variable `t` es la que controla el tiempo que hemos integrado, y su valor inicial es 0. Por otro lado `T` corresponde al `t` final, o hasta que tiempo queremos integrar. 
- **`a`**: parámetros que puede tener el campo de la EDO. Esta declarado como un double con valor 0, pero sería sencillo adaptar-lo a un vector que pudiera contener más de un valor. 
- **Tolerancia**: error màximo permitido del mèetodo Runge-Kutta, que coresponde a la variable `tol`. Por defecto esta puesta en `tol=1e-10`.
- **`e` i `ev`**: `e` corresponde al módulo en norma 1 del error estimado total. Por otro lado `ev` es un vector, en que cada coordenada representa el error estimado cometido en la coordenada correspondente de la solución.  

Por último el valor retornado por la función `rk45` nos dice si ha ocurrido algún problema durante el proceso de integración, donde un retorno de 1 indica que se ha integrado hasta el final sin ningún problema, y un retorno negativo indica que ha ocurrido algún error en el proceso, y se muestra el mensaje correspondiente.  
Comentar que mientras es posible mostrar el error estimado con `e` y `ev`, no esta porgramado que se muestre por defecto.
# Descripción
##### Funciones presentes en `Funcio.c`
- **camp**: función que hace el papel del campo de la EDO. 

Argumentos: `int n`: dimensión del campo,`double t`: variable que controla el tiempo de integración, `double *x`: pointer que apunta al vector donde se va almacenando el valor de la solución de la EDO, `double *f`: pointer que apunta a un vector que almacena el valor del campo después de recivir el valor de la solución en el tiempo correspondiente, `double *prm`: parámetros que pueda tener el campo.


##### Funciones presentes en `otros.c`
-**Abs**: función que te devuelve el valor absoluto del input. Debido a experiencias pasadas prefiero declarar yo mismo la función valor absoluto y no usar la de `math.h`.  

Argumentos: `double x`: valor del que se quiere el valor absoluto.  

-**minomax**: función que dependiendo del input, te devuelve el mínimo o máximo entre dos valores.  

Arguemntos: `double x, double y`: valores que queremos comparar, `int a`: si `a=-1` la función devuelve el mínimo, si `a=1` te devuelve el máximo.  

-**Start**: función void que acepta un vector y lo convierte en el vector nulo.

Argumentos: `double *v`: vector que queremos convertir en el vector nulo, `int n`: dimensión del vecvtor.

##### Funciones presentes en `rk45.c`
-**rk45**: función que implementa el Runge-Kutta 4(5). 

Argumentos: `double *ev`: vector que almacenará el error estimado en cada coordenada de la solución, `double *e`: error en norma 1 del error estiamdo, `double *t`: variable que controla internamente el tiempo integrado, `double *x`: vector en el que se almacenará la solución,  `double *h`: varaible que controla internamente la longitud del paso de la iteración actual, `double hmin` i `double hmax`: rango de valores permitidos para `h`, `double tol`: tolerancia permitida, `int n`: dimensión del campo de la EDO, `int (*camp)(...)`: campo de la EDO, `double *prm`: paràmetros del campo de la EDO.

## Observaciones
Esta fue una de las prácticas que hice en la optativa de cálculo numérico de 4to curso del grado de Matemáticas en la UAB. 
# Comentarios personales respecto el programa
Esta práctica me hizo sentir mucho más cómodo trabajando con memória dinámica.  
También me obligó a organizar mis programas en varios archivos, cosa que permite no tener un main extremadamente largo. La organización de mis codigos siempre sufre de este problema y esta práctica fue un paso hacia mejorar este aspecto en particular.  
Como mejora del programa, se podria dar más libertad al usuario para introducir la las conduiciones iniciales o incluso la tolerancia.  
