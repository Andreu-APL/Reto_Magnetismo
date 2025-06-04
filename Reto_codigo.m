clear; clc; close all;

%% 1. PARÁMETROS DE ENTRADA (Definidos por el Usuario)
% ------------------------------------------------------------------------
disp('--- Configuración de la Simulación ---');

% Parámetros Físicos (Unidades del Sistema Internacional: kg, m, s, N)
m = 1200;      % Masa de la góndola (kg)
g = 9.81;      % Aceleración debido a la gravedad (m/s^2)
k = 3500;      % Constante de frenado magnético (N*s/m)
               % ¡Este es un valor CRÍTICO para ajustar!
               % k baja = frenado suave (peligro de impacto fuerte)
               % k alta = frenado brusco (peligro de aceleración alta)

% Dimensiones de la Torre (m)
H_torre = 65;             % Altura total desde donde se suelta la góndola
H_inicio_frenado = 25;    % Altura (desde el suelo) donde COMIENZA el frenado
                          % La góndola cae libremente desde H_torre hasta H_inicio_frenado
                          % Luego frena desde H_inicio_frenado hasta el suelo (0 m)

% Condiciones Iniciales
y0_inicial = H_torre;       % Posición inicial (altura)
v0_inicial = 0;             % Velocidad inicial (m/s)

% Límites de Seguridad
max_G_peligrosas = 4.5;     % Límite de aceleración/desaceleración en unidades de 'g'
max_v_impacto = 1.5;        % Límite de velocidad al llegar al suelo (m/s)

% Parámetros para el solver ODE (ode45)
tiempo_max_simulacion_fase = 20; % Tiempo máximo para cada fase de la simulación (s)
                                 % Aumentar si la simulación se corta prematuramente

fprintf('Masa (m): %.1f kg\n', m);
fprintf('Constante de Frenado (k): %.1f Ns/m\n', k);
fprintf('Altura Torre: %.1f m\n', H_torre);
fprintf('Altura Inicio Frenado: %.1f m\n', H_inicio_frenado);
fprintf('Distancia Caída Libre: %.1f m\n', H_torre - H_inicio_frenado);
fprintf('Distancia de Frenado: %.1f m\n\n', H_inicio_frenado);

if H_inicio_frenado >= H_torre
    error('La altura de inicio de frenado debe ser MENOR que la altura total de la torre.');
end
if H_inicio_frenado < 0
    error('La altura de inicio de frenado no puede ser negativa.');
end


%% 2. DEFINICIÓN DE LAS ECUACIONES DE MOVIMIENTO (Funciones para ODE45)
% ------------------------------------------------------------------------
% El solver ode45 necesita una función que describa cómo cambian el estado
% (posición 'y' y velocidad 'v') con el tiempo.
% El estado Y es un vector columna: Y(1) es posición 'y', Y(2) es velocidad 'v'.
% La función debe devolver dY/dt, que es [dy/dt; dv/dt].
% dy/dt = v  (esto es Y(2))
% dv/dt = aceleración (calculada según la fase)

% Ecuación para la Fase 1: Caída Libre (sin frenado magnético)
% dv/dt = -g
eq_caida_libre = @(t, Y) [
    Y(2);       % dy/dt = v
    -g          % dv/dt = -g (aceleración constante hacia abajo)
];

% Ecuación para la Fase 2: Con Frenado Magnético
% dv/dt = -g - (k/m)*v
% Aquí, Y(2) es la velocidad 'v'. La fuerza de frenado es -k*v (opuesta a v)
eq_con_frenado = @(t, Y, m_pasada, g_pasada, k_pasada) [
    Y(2);                                   % dy/dt = v
    -g_pasada - (k_pasada/m_pasada) * Y(2)  % dv/dt = -g - (k/m)*v
];

%% 3. FASE 1: SIMULACIÓN DE LA CAÍDA LIBRE
% ------------------------------------------------------------------------
disp('--- Iniciando Fase 1: Caída Libre ---');

% Condiciones iniciales para la Fase 1
estado_inicial_f1 = [y0_inicial; v0_inicial]; % [altura_inicial; velocidad_inicial]
t_span_f1 = [0, tiempo_max_simulacion_fase];   % Intervalo de tiempo para simular

% Opciones para ode45: Queremos detener la simulación de esta fase
% cuando la góndola alcance la altura H_inicio_frenado.
% Esto se hace con una 'función de evento'.
evento_alcanzar_altura_frenado = @(t, Y) event_stop_at_height(t, Y, H_inicio_frenado);
opciones_f1 = odeset('Events', evento_alcanzar_altura_frenado);

% Llamada al solver ODE45 para la Fase 1
% T1: Vector de tiempos
% SOL1: Matriz de soluciones (col1=posición 'y', col2=velocidad 'v')
% te1, ye1, ie1: Información sobre el evento que detuvo la simulación
[T1, SOL1, te1, ye1, ie1] = ode45(eq_caida_libre, t_span_f1, estado_inicial_f1, opciones_f1);

% Extraer resultados de la Fase 1
Y_pos_f1 = SOL1(:,1); % Posiciones (alturas) en la Fase 1
V_vel_f1 = SOL1(:,2); % Velocidades en la Fase 1
A_acc_f1 = -g * ones(size(T1)); % Aceleración en Fase 1 es constante (-g)

% Estado al FINAL de la Fase 1 (será el INICIO de la Fase 2)
t_fin_f1 = T1(end);
y_fin_f1 = Y_pos_f1(end);
v_fin_f1 = V_vel_f1(end);

fprintf('Fase 1 completada.\n');
fprintf('  Tiempo transcurrido: %.2f s\n', t_fin_f1);
fprintf('  Altura al final de Fase 1: %.2f m (objetivo: %.2f m)\n', y_fin_f1, H_inicio_frenado);
fprintf('  Velocidad al final de Fase 1: %.2f m/s (%.2f km/h)\n\n', v_fin_f1, abs(v_fin_f1)*3.6);

% Inicializar variables combinadas con los resultados de la Fase 1
% Estas se completarán con los resultados de la Fase 2
T_combinado = T1;
Y_pos_combinado = Y_pos_f1;
V_vel_combinado = V_vel_f1;
A_acc_combinado = A_acc_f1;

% Variables para almacenar los resultados de la Fase 2
T2_abs = []; % Tiempos absolutos de la Fase 2
Y_pos_f2 = []; V_vel_f2 = []; A_acc_f2 = [];

%% 4. FASE 2: SIMULACIÓN CON FRENADO MAGNÉTICO
% ------------------------------------------------------------------------
% Esta fase solo se ejecuta si la Fase 1 realmente alcanzó la altura de frenado.
% ie1 (índice de evento) estará vacío si el evento no ocurrió.
if isempty(ie1)
    warning('ADVERTENCIA: La Fase 1 no alcanzó la altura de inicio de frenado.');
    fprintf('La góndola pudo haber golpeado el suelo en caída libre o el tiempo_max_simulacion_fase es muy corto.\n');
    % En este caso, la velocidad final es la última velocidad de la Fase 1
    velocidad_final_impacto = V_vel_f1(end);
    altura_final_impacto = Y_pos_f1(end);
else
    disp('--- Iniciando Fase 2: Frenado Magnético ---');
    
    % Condiciones iniciales para la Fase 2 (vienen del final de la Fase 1)
    estado_inicial_f2 = [y_fin_f1; v_fin_f1];
    % El tiempo para ode45 siempre empieza en 0 para su simulación interna.
    % Luego ajustaremos T2 para que sea absoluto.
    t_span_f2 = [0, tiempo_max_simulacion_fase];

    % Opciones para ode45: Detener cuando la góndola llegue al suelo (y=0)
    evento_llegar_al_suelo = @(t,Y) event_stop_at_height(t, Y, 0); % Detener a altura 0
    opciones_f2 = odeset('Events', evento_llegar_al_suelo);

    % Llamada al solver ODE45 para la Fase 2
    % Nota: Pasamos m, g, k como parámetros adicionales a eq_con_frenado
    [T2_rel, SOL2, te2, ye2, ie2] = ode45(@(t,Y) eq_con_frenado(t, Y, m, g, k), ...
                                          t_span_f2, estado_inicial_f2, opciones_f2);

    % Extraer resultados de la Fase 2
    Y_pos_f2 = SOL2(:,1);
    V_vel_f2 = SOL2(:,2);
    A_acc_f2 = -g - (k/m) * V_vel_f2; % Aceleración calculada con la ecuación de frenado

    % Ajustar los tiempos de la Fase 2 para que sean absolutos (continúen desde T1)
    T2_abs = t_fin_f1 + T2_rel;

    fprintf('Fase 2 completada.\n');
    if isempty(ie2)
        warning('ADVERTENCIA: La Fase 2 no terminó llegando al suelo (y=0).');
        fprintf('Puede que la góndola se haya detenido antes o el tiempo_max_simulacion_fase es muy corto.\n');
        velocidad_final_impacto = V_vel_f2(end); % Última velocidad registrada
        altura_final_impacto = Y_pos_f2(end);    % Última altura registrada
    else
        fprintf('  Góndola llegó al suelo.\n');
        velocidad_final_impacto = V_vel_f2(end);
        altura_final_impacto = Y_pos_f2(end);
    end
    fprintf('  Tiempo total del viaje: %.2f s\n', T2_abs(end));
    fprintf('  Altura final: %.2f m\n', altura_final_impacto);
    fprintf('  Velocidad final (al llegar al suelo/parar): %.2f m/s (%.2f km/h)\n\n', ...
            velocidad_final_impacto, abs(velocidad_final_impacto)*3.6);
    
    % Combinar resultados de Fase 1 y Fase 2
    % Evitar duplicar el punto de transición usando (2:end) para la Fase 2
    if ~isempty(T2_abs) % Solo combinar si la Fase 2 produjo resultados
        T_combinado = [T1; T2_abs(2:end)];
        Y_pos_combinado = [Y_pos_f1; Y_pos_f2(2:end)];
        V_vel_combinado = [V_vel_f1; V_vel_f2(2:end)];
        A_acc_combinado = [A_acc_f1; A_acc_f2(2:end)];
    end
end

%% 5. ANÁLISIS DE SEGURIDAD
% ------------------------------------------------------------------------
disp('--- Análisis de Seguridad ---');

% Encontrar la aceleración máxima positiva (mayor frenado)
% y la mínima (más negativa, usualmente -g en caída libre)
max_aceleracion = max(A_acc_combinado);
min_aceleracion = min(A_acc_combinado); % Será -g durante la caída libre

% Convertir a unidades de 'g'
% La aceleración "sentida" más fuerte es usualmente el frenado (aceleración positiva)
max_aceleracion_Gs = max_aceleracion / g;

fprintf('Velocidad máxima (magnitud) alcanzada: %.2f m/s (%.2f km/h)\n', ...
        max(abs(V_vel_combinado)), max(abs(V_vel_combinado))*3.6);
fprintf('Aceleración máxima (positiva, durante frenado): %.2f m/s^2 (%.2f G)\n', ...
        max_aceleracion, max_aceleracion_Gs);
fprintf('Aceleración mínima (más negativa): %.2f m/s^2 (%.2f G)\n', ...
        min_aceleracion, min_aceleracion/g);

% Advertencias de Seguridad
advertencias = {};
if max_aceleracion_Gs > max_G_peligrosas
    msg = sprintf('¡PELIGRO! Aceleración máxima de frenado (%.2f G) SUPERA el límite de %.1f G.', ...
                  max_aceleracion_Gs, max_G_peligrosas);
    advertencias{end+1} = msg;
    disp(msg);
end

if abs(velocidad_final_impacto) > max_v_impacto && altura_final_impacto < 0.1 % Solo si realmente tocó el suelo
    msg = sprintf('¡PELIGRO! Velocidad de impacto (%.2f m/s) SUPERA el límite de %.1f m/s.', ...
                  abs(velocidad_final_impacto), max_v_impacto);
    advertencias{end+1} = msg;
    disp(msg);
end

if altura_final_impacto > 0.1 && ~(isempty(ie1) && isempty(ie2)) % Si no se detuvo por tiempo agotado en ambas fases
     % (es decir, si al menos una fase completó su simulación ODE)
    if isempty(ie2) && ~isempty(T2_abs) % Si Fase 2 corrió pero no llegó al suelo
        msg = sprintf('ADVERTENCIA: La góndola NO LLEGÓ AL SUELO. Se detuvo a %.2f m de altura.', altura_final_impacto);
        advertencias{end+1} = msg;
        disp(msg);
    end
end

if isempty(advertencias)
    disp('Simulación DENTRO de los límites de seguridad establecidos.');
end
fprintf('\n');

%% 6. GRÁFICAS
% ------------------------------------------------------------------------
disp('--- Generando Gráficas ---');

% Figura 1: Velocidad y Aceleración vs. Altura (como solicitado)
figure('Name', 'Resultados Principales: vs. Altura', 'NumberTitle', 'off');

% Subplot 1: Rapidez vs. Altura
subplot(1,2,1);
plot(Y_pos_combinado, abs(V_vel_combinado), 'b-', 'LineWidth', 1.5);
hold on;
% Marcar punto de inicio de frenado si la Fase 2 ocurrió
if ~isempty(ie1)
    plot(y_fin_f1, abs(v_fin_f1), 'ko', 'MarkerFaceColor', 'y', 'DisplayName', 'Inicio Frenado');
    xline(H_inicio_frenado, 'k--', 'Inicio Frenado');
end
hold off;
xlabel('Altura (m)');
ylabel('Rapidez (|v|, m/s)');
title('Rapidez vs. Altura');
grid on;
xlim([0, H_torre * 1.05]); % Un poco más de la altura de la torre
ylim_rapidez = ylim; % Guardar límites actuales
ylim([0, ylim_rapidez(2)]); % Asegurar que el eje Y de rapidez empiece en 0

% Subplot 2: Aceleración vs. Altura (en Gs)
subplot(1,2,2);
plot(Y_pos_combinado, A_acc_combinado / g, 'r-', 'LineWidth', 1.5);
hold on;
% Marcar punto de inicio de frenado
if ~isempty(ie1)
    % Encontrar la aceleración justo al inicio del frenado
    idx_frenado_inicio = find(T_combinado >= t_fin_f1, 1, 'first');
    if ~isempty(idx_frenado_inicio)
        plot(Y_pos_combinado(idx_frenado_inicio), A_acc_combinado(idx_frenado_inicio)/g, ...
            'ko', 'MarkerFaceColor', 'y', 'DisplayName', 'Inicio Frenado');
    end
    xline(H_inicio_frenado, 'k--', 'Inicio Frenado');
end
% Líneas de referencia para aceleración
yline(0, 'k:', '0 G');
yline(-1, 'b:', '-1 G (Caída Libre)');
yline(max_G_peligrosas, 'm--', sprintf('Límite %.1f G', max_G_peligrosas), 'LabelVerticalAlignment', 'bottom');
yline(-max_G_peligrosas, 'm--', sprintf('Límite -%.1f G', max_G_peligrosas), 'LabelVerticalAlignment', 'top'); % Límite simétrico por si acaso
hold off;
xlabel('Altura (m)');
ylabel('Aceleración (unidades de g)');
title('Aceleración vs. Altura');
grid on;
xlim([0, H_torre * 1.05]);

sgtitle(sprintf('Simulación Torre de Caída (m=%.0fkg, k=%.0fNs/m, H_{freno}=%.0fm)', m, k, H_inicio_frenado), 'FontWeight', 'bold');

% Figura 2: Velocidad y Aceleración vs. Tiempo (útil para entender la dinámica)
figure('Name', 'Perfiles Temporales', 'NumberTitle', 'off');

% Subplot 1: Velocidad vs. Tiempo
subplot(2,1,1);
plot(T_combinado, V_vel_combinado, 'b-', 'LineWidth', 1.5);
hold on;
if ~isempty(ie1)
    plot(t_fin_f1, v_fin_f1, 'ko', 'MarkerFaceColor', 'y');
    xline(t_fin_f1, 'k--', 'Inicio Frenado');
end
hold off;
xlabel('Tiempo (s)');
ylabel('Velocidad (m/s)');
title('Velocidad vs. Tiempo');
grid on;

% Subplot 2: Aceleración vs. Tiempo (en Gs)
subplot(2,1,2);
plot(T_combinado, A_acc_combinado / g, 'r-', 'LineWidth', 1.5);
hold on;
if ~isempty(ie1)
    idx_frenado_inicio = find(T_combinado >= t_fin_f1, 1, 'first');
    if ~isempty(idx_frenado_inicio)
        plot(T_combinado(idx_frenado_inicio), A_acc_combinado(idx_frenado_inicio)/g, ...
            'ko', 'MarkerFaceColor', 'y');
    end
    xline(t_fin_f1, 'k--', 'Inicio Frenado');
end
yline(0, 'k:', '0 G');
yline(-1, 'b:', '-1 G (Caída Libre)');
yline(max_G_peligrosas, 'm--', sprintf('Límite %.1f G', max_G_peligrosas), 'LabelVerticalAlignment', 'bottom');
hold off;
xlabel('Tiempo (s)');
ylabel('Aceleración (unidades de g)');
title('Aceleración vs. Tiempo');
grid on;

sgtitle(sprintf('Perfiles Temporales (m=%.0fkg, k=%.0fNs/m, H_{freno}=%.0fm)', m, k, H_inicio_frenado), 'FontWeight', 'bold');

disp('--- Simulación Finalizada ---');

%% 7. FUNCIÓN DE EVENTO (Usada por ODE45 para detener la simulación)
% ------------------------------------------------------------------------
% Esta función se llama en cada paso de la integración de ode45.
% Debe devolver:
%   value: El valor de una función que queremos que sea cero.
%   isterminal: 1 si la integración debe detenerse cuando value=0, 0 si no.
%   direction: 0 para detectar ceros sin importar la dirección,
%              1 para detectar ceros cuando la función 'value' está aumentando,
%             -1 para detectar ceros cuando la función 'value' está disminuyendo.

function [value, isterminal, direction] = event_stop_at_height(t, Y, target_height)
    % Y(1) es la posición actual (altura)
    value = Y(1) - target_height; % Queremos que (altura_actual - altura_objetivo) sea 0
    isterminal = 1;               % Sí, detener la integración
    direction = -1;               % Detener cuando la altura ESTÁ DISMINUYENDO y cruza target_height
                                  % (o cuando Y(1) cruza target_height desde arriba)
                                  % Si target_height fuera mayor que y(0), se usaría direction = 1.
end