# Simulacion de Frenado Magnetico en Torres de Caida

## 🎢 Descripcion del Proyecto

Este proyecto simula computacionalmente la desaceleracion por frenado magnetico (corrientes de Eddy) de una gondola en una torre de caida de un parque de diversiones. La simulacion permite visualizar y analizar el comportamiento de la velocidad y aceleracion de la gondola durante su descenso, considerando diferentes fases del movimiento asi como parametros de la torre.

### Caracteristicas Principales

- Simulación completa del movimiento de la gondola en diferentes fases
- Visualización de velocidad y aceleración en función del tiempo y altura
- Analisis de seguridad con advertencias sobre aceleraciones y velocidades peligrosas
- Implementación en MATLAB con unidades del Sistema Internacional ( SI )
- Personalizacion de parámetros fisicos y geometricos de la torre

## Estructura del Proyecto

El proyecto esta organizado en tres entregables incrementales, cada uno iterando sobre el anterior.

1. **Entregable 1**: Simulacion basica del frenado magnetico
2. **Entregable 2**: Simulacion con caida libre seguida de frenado magnetico
3. **Entregable 3**: Simulacion completa con parametros detallados y advertencias de seguridad

## Como Ejecutar el Codigo

### Requisitos Previos

- MATLAB ( R2019b o mayor )
- Paquete de Ecuaciones Diferenciales de MATLAB ( para la función `ode45` )

### Ejecucion de los Scripts

#### Entregable 1: Simulación Basica

```matlab
% En la consola de MATLAB
entregable1
```

Esta parte simula el movimiento de la gondola bajo la influencia del frenado magnetico, partiendo del reposo. Genera graficas de velocidad y aceleración contra tiempo y altura.

#### Entregable 2: Simulacion con Caida Libre

```matlab
% En la consola de MATLAB
entregable2_corregido_v3
```

Esta parte simula una fase inicial de caida libre seguida por el frenado magnetico. La velocidad final de la caida libre se utiliza como velocidad inicial para la fase de frenado.

#### Entregable 3: Simulacion Completa

```matlab
% En la consola de MATLAB
entregable3_corregido
```

Esta parte implementa la simulacion completa con parametros detallados de la torre, analisis de seguridad y advertencias sobre aceleraciones excesivas o velocidades de impacto que puedan ser peligrosas.

## Parametros Configurables

Todos los scripts permiten la modificación de varios parametros para personalizar la simulación:

```matlab
% Parámetros físicos
m = 500;          % Masa de la góndola (kg)
g_accel = 9.81;   % Aceleración debida a la gravedad (m/s²)

% Parámetros de la torre (Entregable 3)
H_total_torre = 100;        % Altura total de la torre (m)
H_conductora_inicio = 70;   % Altura donde comienza la parte conductora (m)
H_conductora_fin = 10;      % Altura donde termina la parte conductora (m)

% Parámetros de seguridad (Entregable 3)
max_accel_segura_g = 3;     % Máxima aceleración segura (en unidades de g)
max_vel_final_segura = 1.0; % Máxima velocidad segura al llegar al piso (m/s)
```

## Explicación de la Fisica

### Fuerzas Involucradas

En este sistema actuan principalmente dos fuerzas:

1. **Fuerza Gravitacional (Fg)**: La fuerza con la que la Tierra atrae a la gondola.
   ```
   Fg = m * g
   ```
   Donde `m` es la masa de la gondola y `g` es la aceleración debida a la gravedad.

2. **Fuerza de Frenado por Corrientes de Eddy (Fz)**: La fuerza que se opone al movimiento de la gondola, generada por el sistema de frenado magnético. Esta fuerza es proporcional a la velocidad de la gondola.
   ```
   Fz = k * v
   ```
   Donde `k` es una constante de proporcionalidad y `v` es la velocidad de la gondola.

### Ecuacion Diferencial del Movimiento

Aplicando la Segunda Ley de Newton ( F = m·a ), obtenemos:

```
m * (dv/dt) = Fg - Fz
m * (dv/dt) = m*g - k*v
dv/dt = g - (k/m)*v
```

Esta ecuación diferencial describe como cambia la velocidad de la góndola con el tiempo bajo la influencia de la gravedad y el frenado magnético.

### Velocidad Terminal

La velocidad terminal es la velocidad constante que alcanza la gondola cuando la fuerza de frenado magnetico ( Fz ) se iguala en magnitud a la fuerza gravitacional ( Fg ):

```
Fg = Fz
m*g = k*v_terminal
v_terminal = (m*g)/k
```

A esta velocidad, la aceleración es cero y la gondola continúa descendiendo a velocidad constante.

### Fases del Movimiento

La simulacion completa considera tres fases distintas:

1. **Caida Libre Inicial**: La gondola cae bajo la influencia de la gravedad sin resistencia.
   ```
   a = g
   v(t) = g*t
   h(t) = h_inicial - 0.5*g*t²
   ```

2. **Frenado Magnetico**: La gondola experimenta tanto la fuerza gravitacional como la fuerza de frenado magnetico.
   ```
   a(t) = g - (k/m)*v(t)
   ```
   Esta fase se resuelve numericamente usando `ode45` en MATLAB.

3. **Caida Libre Final**: Si la gondola sale de la zona conductora antes de llegar al suelo, continua en caida libre.

## Resultados y Visualizaciones

Los scripts generan varias graficas que muestran:

- Velocidad vs. Tiempo
- Aceleracion vs. Tiempo
- Velocidad y Aceleracion vs. Altura

el Entregable 3 proporciona advertencias de seguridad cuando:
     
- La aceleración maxima excede un limite seguro ( por defecto, 3g )
- La velocidad final al llegar al piso excede un valor seguro ( por defecto, 1.0 m/s )

---

<p align="center">
</p>
