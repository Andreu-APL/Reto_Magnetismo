%%% Main.m

% -------------------- PARTE 1: GRAFICACIÓN DE LAS ESPIRAS --------------------

nl = 5;               % Número total de espiras del solenoide: define cuántas vueltas de alambre hay
N = 20;               % Número de segmentos por espira: usado para discretizar la espira en elementos de corriente
R = 1.5;              % Radio de cada espira (en metros): distancia del centro al alambre
sz = 1;               % Separación entre espiras (paso en el eje Z): determina la longitud del solenoide
I = 300;              % Corriente que circula por las espiras (Amperes): fuente del campo magnético según Biot–Savart
mo = 4*pi*1e-7;       % Permeabilidad magnética del vacío (H/m): constante fundamental que aparece en la ley de Biot–Savart
km = mo * I / (4*pi); % Factor constante para Biot–Savart, simplifica los cálculos posteriores
rw = 0.2;             % Grosor efectivo del alambre (usado para evitar singularidades cuando el punto de cálculo se acerca mucho)
ds = 0.1;             % Resolución espacial para el cálculo del campo: menor = mayor precisión y mayor tiempo de cálculo

% Generar las espiras y elementos de corriente usando función Espiras
% Px, Py, Pz son posiciones; dx, dy, dz son los diferenciales de línea dl
[Px, Py, Pz, dx, dy, dz] = Espiras(nl, N, R, sz);

% -------------------- PARTE 2: VISUALIZACIÓN DEL CAMPO MAGNÉTICO --------------------

plot_option = true;  % Activamos la opción de graficar el campo sin devolver datos
CampoB(ds, km, Px, Py, Pz, dx, dy, nl, N, rw, plot_option);  % Cálculo del campo magnético con ley de Biot–Savart

% -------------------- PARTE 3: SIMULACIÓN DE LA TRAYECTORIA --------------------

% Se desactiva la visualización y se activa el cálculo real de Bz vs z
plot_option = false;
ds = 0.005;  % Mayor resolución para capturar bien el gradiente de campo magnético

% Calculamos la componente Bz del campo y la posición z
[Bz, z] = CampoB(ds, km, Px, Py, Pz, dx, dy, nl, N, rw, plot_option);

% Parámetros del dipolo magnético
mag = 1000;     % Magnitud del momento magnético del objeto (Am^2)
m = 0.009;      % Masa del objeto (kg)
zo = 4.9;       % Posición inicial del objeto sobre el eje z (m)
dt = 0.05;      % Paso de tiempo para la simulación (s)
vz = 0.7;       % Velocidad inicial del objeto (m/s) hacia abajo
gamma = 0.08;   % Coeficiente de fricción (rozamiento viscoso)

% Ejecutar la simulación del movimiento del dipolo magnético bajo fuerzas magnéticas, gravitatorias y fricción
Trayectoria(Bz, z, mag, m, zo, dt, vz, gamma);



%%
%%Función Espiras
%%
function [Px,Py,Pz,dx,dy,dz] = Espiras(nl,N,R,sz)
    dtheta = 2*pi/N;                        % Paso angular para discretizar la espira
    ang = 0:dtheta:(2*pi-dtheta);          % Ángulos alrededor del círculo (evita duplicado en 2π)
    sz=1;
    s=1;
    for i = 1:nl
        % Coordenadas de los puntos de la espira i (proyectada en XY)
        Px(s:s+N-1) = R*cos(ang);          % Coordenada X de los puntos de la espira
        Py(s:s+N-1) = R*sin(ang);          % Coordenada Y
        Pz(s:s+N-1) = -nl/2*sz+(i-1)*sz;   % Coordenada Z (cada espira se coloca a sz metros de la anterior)

        % Componentes del diferencial de línea dl = (dx, dy, 0) para cada segmento de corriente
        dx(s:s+N-1) = -Py(s:s+N-1)*dtheta; % dx = -y dθ, componente tangente
        dy(s:s+N-1) = Px(s:s+N-1)*dtheta;  % dy = x dθ

        s=s+N;
    end
    dz = zeros(1, N*nl);  % No hay componente vertical del diferencial de línea (la corriente circula en planos XY)

    % Graficar vectores de corriente
    figure(1)
    quiver3(Px,Py,Pz,dx,dy,dz,0.5,'-r','LineWidth',2)
    view(-34,33)
    xlabel('x'); ylabel('y'); zlabel('z')
    title('Corriente de espiras')
    axis equal
end

%%Función Campo Magnético

function [Bz,z] = CampoB(ds,km,Px, Py, Pz, dx, dy, nl, N, rw, plot_option)

    x = -5:ds:5; y = x; z = x;         % Definimos espacio tridimensional para calcular campo

    if plot_option
        x=z; y=z;                      % Si solo vamos a graficar, tomamos un plano simétrico XZ
    else
        x=-0.1:0.01:0.1; y=-0.1:0.01:0.1;  % Si vamos a simular trayectoria, nos enfocamos solo en el eje Z
    end

    % Inicialización de arreglos para el campo magnético
    Lx = length(x); Ly = length(y); Lz = length(z);
    dBx = zeros(Lx, Ly, Lz,'single'); dBy = dBx; dBz = dBx;

    % Cálculo de campo con la ley de Biot–Savart
    for i = 1:Lx
        for j = 1:Ly
            for k = 1:Lz
                for l = 1:nl*N
                    % Vector desde elemento de corriente hasta punto de evaluación
                    rx = x(i)-Px(l);
                    ry = y(j)-Py(l);
                    rz = z(k)-Pz(l);

                    % r al cubo (añadiendo rw² para evitar singularidades)
                    r = sqrt(rx^2+ry^2+rz^2+rw^2);
                    r3 = r^3;

                    % Componentes de campo usando Biot–Savart: dB = (μ₀ I / 4π) * (dl × r̂) / r²
                    dBx(i,j,k) = dBx(i,j,k)+km*dy(l)*rz/r3;  % componente X del campo
                    dBy(i,j,k) = dBy(i,j,k)+km*dx(l)*rz/r3;  % componente Y
                    dBz(i,j,k) = dBz(i,j,k)+km*(dx(l)*ry-dy(l)*rx)/r3; % componente Z
                end
            end
        end
    end

    if plot_option 
        % Visualización del campo en plano XZ
        Bmag = sqrt(dBx.^2+dBy.^2+dBz.^2);
        centery = round(Ly/2);
        Bx_xz = squeeze(dBx(:,centery,:));
        Bz_xz = squeeze(dBz(:,centery,:));
        Bxz = squeeze(Bmag(:,centery,:));

        % Mapa de color y líneas de campo
        figure(2);
        hold on;
        pcolor(x,y,(Bxz').^(1/3)); shading interp; colormap jet; colorbar;
        h1 = streamslice(x,z,Bx_xz',Bz_xz',3);  % Líneas de campo
        set(h1,'Color',[0.8 1 0.9]);
        Bz=1;
    else
        % Si estamos en modo simulación, extraemos Bz en el eje Z
        idx_x = ceil(Lx/2);
        idx_y = ceil(Ly/2);
        Bz = squeeze(dBz(idx_x, idx_y,:));

        % Calculamos el gradiente del campo dBz/dz
        dBz_dz_profile = diff(Bz)./diff(z);
        z_mid=z(1:end-1) + diff(z)/2;

        % Graficar gradiente
        figure
        plot(z_mid,dBz_dz_profile, 'r-','LineWidth',2)
        xlabel('z'); ylabel('dBz/dz')
        title('Gradiente del campo magnético de Bz')
    end
end


%%Función Trayectoria

function Trayectoria(Bz,z,mag,m,zo,dt,vz,gamma)

    w=-m*9.81;                      % Fuerza gravitacional hacia abajo (N)
    zm(1)=zo;                      % Posición inicial del dipolo
    zmfree(1)=zo;                  % Posición de objeto en caída libre
    vz(1)=0.7;                     % Velocidad inicial
    vzfree(1)=0;                   % Velocidad inicial en caída libre
    tt(1)=0;                       % Tiempo inicial
    cc =1;

    while zm(cc) > -3              % Simular hasta que el dipolo caiga por debajo de z = -3
        delta = 0.005;
        % Estimamos gradiente del campo en la posición actual
        Bz_forward = interp1(z, Bz, zm(cc) + delta, 'linear', 'extrap');
        Bz_backward = interp1(z, Bz, zm(cc) - delta, 'linear', 'extrap');
        dBz_dz = (Bz_forward - Bz_backward) / (2 * delta);

        Fm(cc) = -mag * dBz_dz;         % Fuerza magnética = -m·∇B (atracción o repulsión según gradiente)
        Ff = -gamma * vz(cc);           % Fuerza de fricción (lineal en la velocidad)
        F(cc) = Fm(cc) + w + Ff;        % Fuerza total: magnética + peso + fricción
        a = F(cc) / m;                  % Aceleración instantánea

        % Actualizamos posición y velocidad con cinemática clásica
        zm(cc+1) = zm(cc) + vz(cc)*dt + 0.5*a*dt^2;
        vz(cc+1) = vz(cc) + a*dt;

        % También calculamos la trayectoria sin campo (solo gravedad)
        zmfree(cc+1) = zmfree(cc) + vzfree(cc)*dt + 0.5*(-9.81)*dt^2;
        vzfree(cc+1) = vzfree(cc) + (-9.81)*dt;

        tt(cc+1) = tt(cc) + dt;
        cc = cc + 1;

        % Condición de parada si el dipolo se "estaciona"
        if abs(vz(cc)) < 1e-3
            break
        end
    end

    % Graficar resultados
    figure(99)
    hold on
    plot(tt,zm, 'r-', 'LineWidth',2)
    plot(tt,zmfree,'b--','LineWidth',2)

    grid on 
    xlabel('Time(s)')
    ylabel('Z position (m)')
    title(' Posicion vs tiempo en un Dipolo magnetico cayendo a través de una espira cargada ')
    legend('Trayectoria con la espira cargada', 'Trayectoria de caída libre','Southwest')
    axis([0 max(tt) -6 6])
end

% Resúmen físico