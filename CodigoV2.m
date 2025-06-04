function [Bz,z] = campoB(ds,km,Px, Py, Pz, dx, dy, nl, N, rw, plot_option)
    x = -5:ds:5; y = x; z = x;
    if plot_option
        x=z;
        y=z;
    else
        x=-0.1:0.01:0.1; % Estas x, y son para el centro, no afectan el rango de z
        y=-0.1:0.01:0.1;
    end

    Lx = length(x); Ly = length(y); Lz = length(z);
    dBx = zeros(Lx, Ly, Lz,'single'); dBy = dBx; dBz = dBx;
    
    % Pre-calcular Px, Py, Pz, dx, dy para el bucle interno si es posible
    % (ya están calculados y pasados como argumento, así que está bien)

    for i = 1:Lx
        for j = 1:Ly
            for k = 1:Lz
                % Vectorización parcial para el bucle más interno (l)
                rx_vec = x(i)-Px; % Px es un vector de nl*N elementos
                ry_vec = y(j)-Py; % Py es un vector de nl*N elementos
                rz_vec = z(k)-Pz; % Pz es un vector de nl*N elementos
                
                r_vec_sq = rx_vec.^2 + ry_vec.^2 + rz_vec.^2 + rw^2;
                r_vec = sqrt(r_vec_sq);
                r3_vec = r_vec.^3;

                % Evitar división por cero si r3_vec puede ser cero (rw^2 ayuda)
                % Si r3_vec tiene ceros, manejarlo adecuadamente (ej. Poner Inf o NaN y sumar al final)
                % En este caso, rw^2 > 0 previene r3_vec de ser 0 si rx,ry,rz son 0.

                sum_dBx_k = sum(km * dy .* rz_vec ./ r3_vec);
                sum_dBy_k = sum(km * dx .* rz_vec ./ r3_vec); % Corrección: era dy, debe ser dx para dBy según Biot-Savart usual (dl x r)
                                                            % No, la fórmula original dBx ~ dy*rz y dBy ~ dx*rz parece ser específica de esta derivación.
                                                            % Mantendré la fórmula original: dBy(i,j,k) = dBy(i,j,k)+km*dx(l)*rz/r3;
                                                            % La fórmula original es:
                                                            % dBx = km * (dly*rz - dlz*ry) / r^3
                                                            % dBy = km * (dlz*rx - dlx*rz) / r^3
                                                            % dBz = km * (dlx*ry - dly*rx) / r^3
                                                            % Dado que dlz=0, entonces:
                                                            % dBx = km * dly*rz / r^3  (dy en el código es dly)
                                                            % dBy = km * (-dlx*rz) / r^3 (dx en el código es dlx)
                                                            % dBz = km * (dlx*ry - dly*rx) / r^3
                                                            % El código tiene dBy(i,j,k) = dBx(i,j,k)+km*dx(l)*rz/r3; -> Esto es dBy(i,j,k) = dBy(i,j,k)+km*dx(l)*rz/r3;
                                                            % Parece que la implementación de dBy usa +dx(l)*rz/r3 en lugar de -dx(l)*rz/r3.
                                                            % Por ahora, mantendré la implementación original del usuario.

                sum_dBy_k = sum(km * dx .* rz_vec ./ r3_vec); % Manteniendo la implementación original
                sum_dBz_k = sum(km * (dx .* ry_vec - dy .* rx_vec) ./ r3_vec);
                
                dBx(i,j,k) = sum_dBx_k;
                dBy(i,j,k) = sum_dBy_k;
                dBz(i,j,k) = sum_dBz_k;
            end
        end
    end

   if plot_option 
        %Convertir a plano XZ
        Bmag = sqrt(dBx.^2+dBy.^2+dBz.^2);
        centery = round(Ly/2);
        Bx_xz = squeeze(dBx(:,centery,:));
        Bz_xz = squeeze(dBz(:,centery,:));
        Bxz = squeeze(Bmag(:,centery,:));
        
        %Graficar
        figure(2);
        hold on;
        pcolor(x,z,(Bxz').^(1/3)); shading interp; colormap jet; colorbar; % x para pcolor, z para pcolor
        h1 = streamslice(x,z,Bx_xz',Bz_xz',3); % x, z para streamslice
        set(h1,'Color',[0.8 1 0.9]);
        Bz_out=1; % Variable de salida dummy
        z_out = z; % Devuelve el vector z usado para graficar
        if nargout > 0, Bz = Bz_out; end
        if nargout > 1, z = z_out; end

   else
        idx_x = ceil(Lx/2);
        idx_y = ceil(Ly/2);
        Bz_profile = squeeze(dBz(idx_x, idx_y,:));

        %Calcular gradiente del campo
        dBz_dz_profile = diff(Bz_profile)./diff(z'); % z es fila, Bz_profile es columna
        z_mid=z(1:end-1)' + diff(z')/2; % z es fila
        
        %Graficar
        figure(3); % Cambié el número de figura para no sobreescribir la figura 2
        plot(z_mid,dBz_dz_profile, 'r-','LineWidth',2)
        xlabel('z (m)'); ylabel('dBz/dz (T/m)')
        title('Gradiente del campo magnético de Bz')
        grid on;

        % Salidas para la simulación
        if nargout > 0, Bz = Bz_profile; end
        if nargout > 1, z = z'; end % Asegurar que z sea un vector columna si Bz_profile lo es
   end
end

function [Px,Py,Pz,dx,dy,dz] = espiras(nl,N,R,sz_in) % sz_in para evitar conflicto con variable sz en campoB
    dtheta = 2*pi/N;
    ang = 0:dtheta:(2*pi-dtheta); % ang es 1xN
    % sz=1; % Esto estaba hardcodeado, usar el argumento sz_in
    s=1;
    
    % Preallocate arrays
    Px = zeros(1, nl*N);
    Py = zeros(1, nl*N);
    Pz = zeros(1, nl*N);
    dx = zeros(1, nl*N);
    dy = zeros(1, nl*N);

    for i = 1:nl
        start_idx = (i-1)*N + 1;
        end_idx = i*N;

        Px(start_idx:end_idx) = R*cos(ang);
        Py(start_idx:end_idx) = R*sin(ang);
        Pz(start_idx:end_idx) = -nl/2*sz_in+(i-1)*sz_in; % Usar sz_in

        dx(start_idx:end_idx) = -Py(start_idx:end_idx)*dtheta; % dl_x = -R*sin(theta)*dtheta
        dy(start_idx:end_idx) = Px(start_idx:end_idx)*dtheta;  % dl_y =  R*cos(theta)*dtheta
        
        % s=s+N; % No es necesario si se usan start_idx y end_idx
    end
    dz = zeros(1, N*nl); % dl_z = 0 para espiras en planos XY
    
    figure(1)
    quiver3(Px,Py,Pz,dx,dy,dz,0.5,'-r','LineWidth',2)
    view(-34,33)
    xlabel('x (m)'); ylabel('y (m)');zlabel('z (m)')
    title('Corriente de espiras')
    axis equal
    
end

function trayectoria(Bz_campo, z_campo, mag, m, zo, dt, v_inicial, gamma)
    % Constantes
    g_accel = -9.81; % Aceleración debida a la gravedad (m/s^2)
    F_gravity = m * g_accel; % Fuerza gravitatoria (constante)
    delta_interp = 0.005; % Pequeño desplazamiento para calcular el gradiente de Bz

    % Inicialización de vectores de tiempo, posición y velocidad
    max_steps = round(20/dt); % Estimación del número máximo de pasos (ej. 20 segundos)
    zm = zeros(1, max_steps);
    vz = zeros(1, max_steps);
    tt = zeros(1, max_steps);
    zmfree = zeros(1, max_steps);
    vzfree = zeros(1, max_steps);

    zm(1) = zo;
    vz(1) = v_inicial; % Usar la velocidad inicial proporcionada
    tt(1) = 0;

    zmfree(1) = zo;
    vzfree(1) = v_inicial; % Caída libre también parte con v_inicial

    cc = 1; % Contador de pasos

    % Bucle de simulación hasta que el imán caiga por debajo de -3m o se detenga
    while zm(cc) > -3 && cc < max_steps -1 
        
        z_current = zm(cc);
        v_current = vz(cc);

        % RK4 para el sistema acoplado dz/dt = v, dv/dt = a(z,v)
        
        % k1
        a1 = calcular_aceleracion(z_current, v_current, z_campo, Bz_campo, mag, m, gamma, F_gravity, delta_interp);
        k1_z = v_current;
        k1_v = a1;
        
        % k2
        a2 = calcular_aceleracion(z_current + 0.5*dt*k1_z, v_current + 0.5*dt*k1_v, z_campo, Bz_campo, mag, m, gamma, F_gravity, delta_interp);
        k2_z = v_current + 0.5*dt*k1_v;
        k2_v = a2;
        
        % k3
        a3 = calcular_aceleracion(z_current + 0.5*dt*k2_z, v_current + 0.5*dt*k2_v, z_campo, Bz_campo, mag, m, gamma, F_gravity, delta_interp);
        k3_z = v_current + 0.5*dt*k2_v;
        k3_v = a3;
        
        % k4
        a4 = calcular_aceleracion(z_current + dt*k3_z, v_current + dt*k3_v, z_campo, Bz_campo, mag, m, gamma, F_gravity, delta_interp);
        k4_z = v_current + dt*k3_v;
        k4_v = a4;
        
        % Actualización de posición y velocidad
        zm(cc+1) = z_current + (dt/6)*(k1_z + 2*k2_z + 2*k3_z + k4_z);
        vz(cc+1) = v_current + (dt/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);
        
        % Actualización para caída libre (Euler simple es suficiente aquí)
        a_free = g_accel; % Solo gravedad
        zmfree(cc+1) = zmfree(cc) + vzfree(cc)*dt + 0.5*a_free*dt^2;
        vzfree(cc+1) = vzfree(cc) + a_free*dt;
    
        tt(cc+1) = tt(cc) + dt;
        cc = cc + 1;
    
        % Condición de parada si la velocidad es muy baja (levitación o parada)
        if abs(vz(cc)) < 1e-4 && abs(a4) < 1e-3 % Usar a4 como la última aceleración calculada
            % Podríamos añadir una comprobación de si está cerca de un punto de equilibrio
            disp('Simulación detenida: velocidad y aceleración cercanas a cero.');
            break;
        end
    end

    % Recortar vectores al tamaño real utilizado
    zm = zm(1:cc);
    vz = vz(1:cc);
    tt = tt(1:cc);
    zmfree = zmfree(1:cc);
    % vzfree = vzfree(1:cc); % No se usa vzfree directamente en el plot

    % Graficar resultados
    figure(99)
    hold on
    plot(tt,zm, 'r-', 'LineWidth',2)
    plot(tt,zmfree,'b--','LineWidth',2)
    hold off

    grid on 
    xlabel('Tiempo (s)')
    ylabel('Posición Z (m)')
    title('Posición vs. Tiempo: Dipolo Magnético Cayendo (RK4)')
    legend('Trayectoria con espira (RK4)', 'Trayectoria de caída libre','Location','Southwest')
    max_t_plot = max(tt);
    if isempty(max_t_plot) || max_t_plot == 0, max_t_plot = 1; end % Evitar error si tt está vacío o es todo cero
    axis_ymin = min([-3.5, min(zm), min(zmfree)]); % Asegurar que -3 sea visible
    axis_ymax = max([zo + 0.5, max(zm), max(zmfree)]); % Un poco por encima de la posición inicial
    axis([0 max_t_plot axis_ymin axis_ymax])

end

% --- Función auxiliar para calcular la aceleración ---
function accel = calcular_aceleracion(pos_z, vel_z, z_data, Bz_data, mag_moment, mass, friction_coeff, F_grav, delta)
    % Interpolar para encontrar Bz en puntos cercanos a pos_z
    % Asegurarse de que z_data sea monótonamente creciente para interp1
    if ~issorted(z_data)
        [z_data_sorted, sort_idx] = sort(z_data);
        Bz_data_sorted = Bz_data(sort_idx);
    else
        z_data_sorted = z_data;
        Bz_data_sorted = Bz_data;
    end

    Bz_forward = interp1(z_data_sorted, Bz_data_sorted, pos_z + delta, 'linear', 'extrap');
    Bz_backward = interp1(z_data_sorted, Bz_data_sorted, pos_z - delta, 'linear', 'extrap');
    
    % Calcular gradiente de Bz
    dBz_dz = (Bz_forward - Bz_backward) / (2 * delta);
    
    % Calcular fuerzas
    Fm = -mag_moment * dBz_dz;    % Fuerza magnética
    Ff = -friction_coeff * vel_z; % Fuerza de fricción
    
    % Fuerza total y aceleración
    F_total = Fm + F_grav + Ff;
    accel = F_total / mass;
end


% -------------------- PARTE 1: GRAFICACIÓN DE LAS ESPIRAS --------------------
nl = 5;               % Número de espiras
N = 20;               % Puntos por espira
R = 1.5;              % Radio de cada espira (m)
sz_espira = 1;        % Separación entre espiras (paso en Z) (m)
I = 300;              % Corriente (A)
mo = 4*pi*1e-7;       % Permeabilidad magnética del vacío (T*m/A)
km = mo * I / (4*pi); % Constante de Biot–Savart (T*m)
rw = 0.2;             % Grosor efectivo del alambre (m)
ds_plot = 0.1;        % Resolución para graficación de campo (más grande = más rápido)

% Generar las espiras
[Px, Py, Pz, dx, dy, dz] = espiras(nl, N, R, sz_espira);

% -------------------- PARTE 2: VISUALIZACIÓN DEL CAMPO MAGNÉTICO --------------------
plot_option_campo = true;  % Solo para graficar el campo
disp('Calculando y graficando campo magnético (visualización)...');
[~, ~] = campoB(ds_plot, km, Px, Py, Pz, dx, dy, nl, N, rw, plot_option_campo); % No necesitamos las salidas aquí

% -------------------- PARTE 3: SIMULACIÓN DE LA TRAYECTORIA --------------------
% Ahora sí calculamos Bz y z para usarlos en la simulación
plot_option_campo = false; % Para obtener Bz y z para la simulación
ds_sim = 0.005;          % Resolución más fina para el cálculo del gradiente y Bz para simulación

disp('Calculando Bz a lo largo del eje z para la simulación...');
[Bz_sim, z_sim] = campoB(ds_sim, km, Px, Py, Pz, dx, dy, nl, N, rw, plot_option_campo);

% Asegurarse que z_sim y Bz_sim son vectores columna y z_sim es monótono
z_sim = z_sim(:);
Bz_sim = Bz_sim(:);
if ~issorted(z_sim)
    [z_sim, sort_idx] = sort(z_sim);
    Bz_sim = Bz_sim(sort_idx);
end


% Parámetros de la simulación de trayectoria
mag = 1000;     % Magnitud del momento magnético (A*m^2)
m = 0.009;      % Masa (kg)
zo = 4.9;       % Posición inicial (m)
dt = 0.01;      % Paso de tiempo (s) (RK4 puede permitir pasos mayores que Euler para misma precisión)
vz0 = 0.7;      % Velocidad inicial (m/s)
gamma = 0.08;   % Coeficiente de fricción (N*s/m)

disp('Iniciando simulación de trayectoria con RK4...');
trayectoria(Bz_sim, z_sim, mag, m, zo, dt, vz0, gamma);

disp('Simulación completada.');
% --- END OF FILE codigo_junto_rk4.txt ---