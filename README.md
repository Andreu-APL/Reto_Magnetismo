# Simulación de Frenado Magnético en Torres de Caída

![Banner](https://i.imgur.com/JQXvmZs.jpg)

## 🎢 Descripción del Proyecto

Este proyecto simula computacionalmente la desaceleración por frenado magnético (corrientes de Eddy) de una góndola en una torre de caída de un parque de diversiones. La simulación permite visualizar y analizar el comportamiento de la velocidad y aceleración de la góndola durante su descenso, considerando diferentes fases del movimiento y parámetros de la torre.

### Características Principales

- ✅ Simulación completa del movimiento de la góndola en diferentes fases
- ✅ Visualización de velocidad y aceleración en función del tiempo y altura
- ✅ Análisis de seguridad con advertencias sobre aceleraciones y velocidades peligrosas
- ✅ Implementación en MATLAB con unidades del Sistema Internacional (SI)
- ✅ Personalización de parámetros físicos y geométricos de la torre

## 📊 Estructura del Proyecto

El proyecto está organizado en tres entregables incrementales, cada uno añadiendo complejidad y funcionalidades:

1. **Entregable 1**: Simulación básica del frenado magnético
2. **Entregable 2**: Simulación con caída libre seguida de frenado magnético
3. **Entregable 3**: Simulación completa con parámetros detallados y advertencias de seguridad

## 🚀 Cómo Ejecutar el Código

### Requisitos Previos

- MATLAB (recomendado R2019b o superior)
- Paquete de Ecuaciones Diferenciales de MATLAB (para la función `ode45`)

### Ejecución de los Scripts

#### Entregable 1: Simulación Básica

```matlab
% En la consola de MATLAB
entregable1
```

Este script simula el movimiento de la góndola bajo la influencia del frenado magnético, partiendo del reposo. Genera gráficas de velocidad y aceleración contra tiempo y altura.

#### Entregable 2: Simulación con Caída Libre

```matlab
% En la consola de MATLAB
entregable2_corregido_v3
```

Este script simula una fase inicial de caída libre seguida por el frenado magnético. La velocidad final de la caída libre se utiliza como velocidad inicial para la fase de frenado.

#### Entregable 3: Simulación Completa

```matlab
% En la consola de MATLAB
entregable3_corregido
```

Este script implementa la simulación completa con parámetros detallados de la torre, análisis de seguridad y advertencias sobre aceleraciones excesivas o velocidades de impacto peligrosas.

## 📝 Parámetros Configurables

Todos los scripts permiten la modificación de varios parámetros para personalizar la simulación:

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

## 🧲 Explicación de la Física

### Fuerzas Involucradas

En este sistema actúan principalmente dos fuerzas:

1. **Fuerza Gravitacional (Fg)**: La fuerza con la que la Tierra atrae a la góndola.
   ```
   Fg = m * g
   ```
   Donde `m` es la masa de la góndola y `g` es la aceleración debida a la gravedad.

2. **Fuerza de Frenado por Corrientes de Eddy (Fz)**: La fuerza que se opone al movimiento de la góndola, generada por el sistema de frenado magnético. Esta fuerza es proporcional a la velocidad de la góndola.
   ```
   Fz = k * v
   ```
   Donde `k` es una constante de proporcionalidad y `v` es la velocidad de la góndola.

### Ecuación Diferencial del Movimiento

Aplicando la Segunda Ley de Newton (F = m·a), obtenemos:

```
m * (dv/dt) = Fg - Fz
m * (dv/dt) = m*g - k*v
dv/dt = g - (k/m)*v
```

Esta ecuación diferencial describe cómo cambia la velocidad de la góndola con el tiempo bajo la influencia de la gravedad y el frenado magnético.

### Velocidad Terminal

La velocidad terminal es la velocidad constante que alcanza la góndola cuando la fuerza de frenado magnético (Fz) se iguala en magnitud a la fuerza gravitacional (Fg):

```
Fg = Fz
m*g = k*v_terminal
v_terminal = (m*g)/k
```

A esta velocidad, la aceleración es cero y la góndola continúa descendiendo a velocidad constante.

### Fases del Movimiento

La simulación completa considera tres fases distintas:

1. **Caída Libre Inicial**: La góndola cae bajo la influencia de la gravedad sin resistencia.
   ```
   a = g
   v(t) = g*t
   h(t) = h_inicial - 0.5*g*t²
   ```

2. **Frenado Magnético**: La góndola experimenta tanto la fuerza gravitacional como la fuerza de frenado magnético.
   ```
   a(t) = g - (k/m)*v(t)
   ```
   Esta fase se resuelve numéricamente usando `ode45` en MATLAB.

3. **Caída Libre Final** (si aplica): Si la góndola sale de la zona conductora antes de llegar al suelo, continúa en caída libre.

## 📈 Resultados y Visualizaciones

Los scripts generan varias gráficas que muestran:

- Velocidad vs. Tiempo
- Aceleración vs. Tiempo
- Velocidad y Aceleración vs. Altura

Además, el Entregable 3 proporciona advertencias de seguridad cuando:
- La aceleración máxima excede un umbral seguro (por defecto, 3g)
- La velocidad final al llegar al piso excede un valor seguro (por defecto, 1.0 m/s)

## 🔍 Análisis de Seguridad

El análisis de seguridad implementado en el Entregable 3 es crucial para el diseño de torres de caída reales. Las aceleraciones excesivas pueden ser peligrosas para los pasajeros, y una velocidad de impacto demasiado alta al final del recorrido podría causar daños estructurales o lesiones.

## 📚 Referencias

- Serway, R. A., & Jewett, J. W. (2018). Physics for Scientists and Engineers. Cengage Learning.
- Griffiths, D. J. (2017). Introduction to Electrodynamics. Cambridge University Press.
- MATLAB Documentation: [ode45](https://www.mathworks.com/help/matlab/ref/ode45.html)

---

<p align="center">
<i>Desarrollado como proyecto de simulación física para aplicaciones de ingeniería</i>
</p>
