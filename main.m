%%% main.m

% -------------------- PARTE 1: GRAFICACIÓN DE LAS ESPIRAS --------------------
nl = 5;               % Número de espiras
N = 20;               % Puntos por espira
R = 1.5;              % Radio de cada espira
sz = 1;               % Separación entre espiras (paso en Z)
I = 300;              % Corriente (A)
mo = 4*pi*1e-7;       % Permeabilidad magnética del vacío
km = mo * I / (4*pi); % Constante de Biot–Savart
rw = 0.2;             % Grosor efectivo del alambre
ds = 0.1;             % Resolución para graficación de campo (más grande = más rápido)

% Generar las espiras
[Px, Py, Pz, dx, dy, dz] = espiras(nl, N, R, sz);

% -------------------- PARTE 2: VISUALIZACIÓN DEL CAMPO MAGNÉTICO --------------------
plot_option = true;  % Solo para graficar el campo (sin devolver z)
campoB(ds, km, Px, Py, Pz, dx, dy, nl, N, rw, plot_option);

% -------------------- PARTE 3: SIMULACIÓN DE LA TRAYECTORIA --------------------
% Ahora sí calculamos Bz y z para usarlos en la simulación
plot_option = false;
ds = 0.005;  % Resolución más fina para el cálculo del gradiente

[Bz, z] = campoB(ds, km, Px, Py, Pz, dx, dy, nl, N, rw, plot_option);

% Parámetros del dipolo magnético
mag = 1000;     % Magnitud del momento magnético
m = 0.009;      % Masa (kg)
zo = 4.9;       % Posición inicial
dt = 0.05;      % Paso de tiempo (s)
vz = 0.7;       % Velocidad inicial (m/s)
gamma = 0.08;   % Coeficiente de fricción

% Ejecutar la simulación
trayectoria(Bz, z, mag, m, zo, dt, vz, gamma);
