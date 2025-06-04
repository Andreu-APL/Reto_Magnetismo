clear; clc; close all;

m = 1200;
g = 9.81;
k = 3500;

H_torre = 65;
H_inicio_frenado = 25;

y0_inicial = H_torre;
v0_inicial = 0;

max_G_peligrosas = 4.5;
max_v_impacto = 1.5;

tiempo_max_simulacion_fase = 20;

disp('--- Configuracion de la Simulacion ---');
fprintf('Masa (m): %.1f kg\n', m);
fprintf('Constante de Frenado (k): %.1f Ns/m\n', k);
fprintf('Altura Torre: %.1f m\n', H_torre);
fprintf('Altura Inicio Frenado: %.1f m\n', H_inicio_frenado);
fprintf('Distancia Caida Libre: %.1f m\n', H_torre - H_inicio_frenado);
fprintf('Distancia de Frenado: %.1f m\n\n', H_inicio_frenado);

if H_inicio_frenado >= H_torre
    error('La altura de inicio de frenado debe ser MENOR que la altura total de la torre.');
end
if H_inicio_frenado < 0
    error('La altura de inicio de frenado no puede ser negativa.');
end

eq_caida_libre = @(t, Y) [
    Y(2);
    -g
];

eq_con_frenado = @(t, Y, m_pasada, g_pasada, k_pasada) [
    Y(2);
    -g_pasada - (k_pasada/m_pasada) * Y(2)
];

disp('--- Iniciando Fase 1: Caida Libre ---');

estado_inicial_f1 = [y0_inicial; v0_inicial];
t_span_f1 = [0, tiempo_max_simulacion_fase];

evento_alcanzar_altura_frenado = @(t, Y) event_stop_at_height(t, Y, H_inicio_frenado);
opciones_f1 = odeset('Events', evento_alcanzar_altura_frenado);

[T1, SOL1, te1, ye1, ie1] = ode45(eq_caida_libre, t_span_f1, estado_inicial_f1, opciones_f1);

Y_pos_f1 = SOL1(:,1);
V_vel_f1 = SOL1(:,2);
A_acc_f1 = -g * ones(size(T1));

t_fin_f1 = T1(end);
y_fin_f1 = Y_pos_f1(end);
v_fin_f1 = V_vel_f1(end);

fprintf('Fase 1 completada.\n');
fprintf('  Tiempo transcurrido: %.2f s\n', t_fin_f1);
fprintf('  Altura al final de Fase 1: %.2f m (objetivo: %.2f m)\n', y_fin_f1, H_inicio_frenado);
fprintf('  Velocidad al final de Fase 1: %.2f m/s (%.2f km/h)\n\n', v_fin_f1, abs(v_fin_f1)*3.6);

T_combinado = T1;
Y_pos_combinado = Y_pos_f1;
V_vel_combinado = V_vel_f1;
A_acc_combinado = A_acc_f1;

T2_abs = [];
Y_pos_f2 = []; V_vel_f2 = []; A_acc_f2 = [];
velocidad_final_impacto = 0; % Inicializar
altura_final_impacto = 0; % Inicializar

if isempty(ie1)
    warning('ADVERTENCIA: La Fase 1 no alcanzo la altura de inicio de frenado.');
    fprintf('La gondola pudo haber golpeado el suelo en caida libre o el tiempo_max_simulacion_fase es muy corto.\n');
    velocidad_final_impacto = V_vel_f1(end);
    altura_final_impacto = Y_pos_f1(end);
else
    disp('--- Iniciando Fase 2: Frenado Magnetico ---');
    
    estado_inicial_f2 = [y_fin_f1; v_fin_f1];
    t_span_f2 = [0, tiempo_max_simulacion_fase];

    evento_llegar_al_suelo = @(t,Y) event_stop_at_height(t, Y, 0);
    opciones_f2 = odeset('Events', evento_llegar_al_suelo);

    [T2_rel, SOL2, te2, ye2, ie2] = ode45(@(t,Y) eq_con_frenado(t, Y, m, g, k), ...
                                          t_span_f2, estado_inicial_f2, opciones_f2);

    Y_pos_f2 = SOL2(:,1);
    V_vel_f2 = SOL2(:,2);
    A_acc_f2 = -g - (k/m) * V_vel_f2;

    T2_abs = t_fin_f1 + T2_rel;

    fprintf('Fase 2 completada.\n');
    if isempty(ie2)
        warning('ADVERTENCIA: La Fase 2 no termino llegando al suelo (y=0).');
        fprintf('Puede que la gondola se haya detenido antes o el tiempo_max_simulacion_fase es muy corto.\n');
        velocidad_final_impacto = V_vel_f2(end);
        altura_final_impacto = Y_pos_f2(end);
    else
        fprintf('  Gondola llego al suelo.\n');
        velocidad_final_impacto = V_vel_f2(end);
        altura_final_impacto = Y_pos_f2(end);
    end
    fprintf('  Tiempo total del viaje: %.2f s\n', T2_abs(end));
    fprintf('  Altura final: %.2f m\n', altura_final_impacto);
    fprintf('  Velocidad final (al llegar al suelo/parar): %.2f m/s (%.2f km/h)\n\n', ...
            velocidad_final_impacto, abs(velocidad_final_impacto)*3.6);
    
    if ~isempty(T2_abs)
        T_combinado = [T1; T2_abs(2:end)];
        Y_pos_combinado = [Y_pos_f1; Y_pos_f2(2:end)];
        V_vel_combinado = [V_vel_f1; V_vel_f2(2:end)];
        A_acc_combinado = [A_acc_f1; A_acc_f2(2:end)];
    end
end

disp('--- Analisis de Seguridad ---');

max_aceleracion = 0; % Inicializar
if ~isempty(A_acc_combinado)
    max_aceleracion = max(A_acc_combinado);
    min_aceleracion = min(A_acc_combinado);
else % Caso donde solo hubo Fase 1 y no se completó
    max_aceleracion = -g;
    min_aceleracion = -g;
end


max_aceleracion_Gs = max_aceleracion / g;

fprintf('Velocidad maxima (magnitud) alcanzada: %.2f m/s (%.2f km/h)\n', ...
        max(abs(V_vel_combinado)), max(abs(V_vel_combinado))*3.6);
fprintf('Aceleracion maxima (positiva, durante frenado): %.2f m/s^2 (%.2f G)\n', ...
        max_aceleracion, max_aceleracion_Gs);
fprintf('Aceleracion minima (mas negativa): %.2f m/s^2 (%.2f G)\n', ...
        min_aceleracion, min_aceleracion/g);

advertencias = {};
if max_aceleracion_Gs > max_G_peligrosas
    msg = sprintf('PELIGRO! Aceleracion maxima de frenado (%.2f G) SUPERA el limite de %.1f G.', ...
                  max_aceleracion_Gs, max_G_peligrosas);
    advertencias{end+1} = msg;
    disp(msg);
end

if abs(velocidad_final_impacto) > max_v_impacto && altura_final_impacto < 0.1
    msg = sprintf('PELIGRO! Velocidad de impacto (%.2f m/s) SUPERA el limite de %.1f m/s.', ...
                  abs(velocidad_final_impacto), max_v_impacto);
    advertencias{end+1} = msg;
    disp(msg);
end

if altura_final_impacto > 0.1 && ~(isempty(ie1) && isempty(ie2))
    if isempty(ie2) && ~isempty(T2_abs)
        msg = sprintf('ADVERTENCIA: La gondola NO LLEGO AL SUELO. Se detuvo a %.2f m de altura.', altura_final_impacto);
        advertencias{end+1} = msg;
        disp(msg);
    end
end

if isempty(advertencias)
    disp('Simulacion DENTRO de los limites de seguridad establecidos.');
end
fprintf('\n');

disp('--- Generando Graficas ---');

figure('Name', 'Resultados Principales: vs. Altura', 'NumberTitle', 'off');

subplot(1,2,1);
plot(Y_pos_combinado, abs(V_vel_combinado), 'b-', 'LineWidth', 1.5);
hold on;
if ~isempty(ie1)
    plot(y_fin_f1, abs(v_fin_f1), 'ko', 'MarkerFaceColor', 'y');
    xline(H_inicio_frenado, 'k--', 'Inicio Frenado');
end
hold off;
xlabel('Altura (m)');
ylabel('Rapidez (|v|, m/s)');
title('Rapidez vs. Altura');
grid on;
xlim([0, H_torre * 1.05]);
ylim_rapidez = ylim;
ylim([0, ylim_rapidez(2)]);

subplot(1,2,2);
plot(Y_pos_combinado, A_acc_combinado / g, 'r-', 'LineWidth', 1.5);
hold on;
if ~isempty(ie1)
    idx_frenado_inicio = find(T_combinado >= t_fin_f1, 1, 'first');
    if ~isempty(idx_frenado_inicio)
        plot(Y_pos_combinado(idx_frenado_inicio), A_acc_combinado(idx_frenado_inicio)/g, ...
            'ko', 'MarkerFaceColor', 'y');
    end
    xline(H_inicio_frenado, 'k--', 'Inicio Frenado');
end
yline(0, 'k:', '0 G');
yline(-1, 'b:', '-1 G (Caida Libre)');
yline(max_G_peligrosas, 'm--', sprintf('Limite %.1f G', max_G_peligrosas), 'LabelVerticalAlignment', 'bottom');
yline(-max_G_peligrosas, 'm--', sprintf('Limite -%.1f G', max_G_peligrosas), 'LabelVerticalAlignment', 'top');
hold off;
xlabel('Altura (m)');
ylabel('Aceleracion (unidades de g)');
title('Aceleracion vs. Altura');
grid on;
xlim([0, H_torre * 1.05]);

sgtitle(sprintf('Simulacion Torre de Caida (m=%.0fkg, k=%.0fNs/m, H_freno=%.0fm)', m, k, H_inicio_frenado), 'FontWeight', 'bold');

figure('Name', 'Perfiles Temporales', 'NumberTitle', 'off');

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
yline(-1, 'b:', '-1 G (Caida Libre)');
yline(max_G_peligrosas, 'm--', sprintf('Limite %.1f G', max_G_peligrosas), 'LabelVerticalAlignment', 'bottom');
hold off;
xlabel('Tiempo (s)');
ylabel('Aceleracion (unidades de g)');
title('Aceleracion vs. Tiempo');
grid on;

sgtitle(sprintf('Perfiles Temporales (m=%.0fkg, k=%.0fNs/m, H_freno=%.0fm)', m, k, H_inicio_frenado), 'FontWeight', 'bold');

disp('--- Simulacion Finalizada ---');

function [value, isterminal, direction] = event_stop_at_height(t, Y, target_height)
    value = Y(1) - target_height;
    isterminal = 1;
    direction = -1;
end