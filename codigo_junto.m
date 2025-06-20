
function [Bz,z] = campoB(ds,km,Px, Py, Pz, dx, dy, nl, N, rw, plot_option)
    x = -5:ds:5; y = x; z = x;
    if plot_option
        x=z;
        y=z;
    else
        x=-0.1:0.01:0.1;
        y=-0.1:0.01:0.1;
    end

    Lx = length(x); Ly = length(y); Lz = length(z);
    dBx = zeros(Lx, Ly, Lz,'single'); dBy = dBx; dBz = dBx;
    
    for i = 1:Lx
        for j = 1:Ly
            for k = 1:Lz
                for l = 1:nl*N
                    rx = x(i)-Px(l);
                    ry = y(j)-Py(l);
                    rz = z(k)-Pz(l);
                    
                    r = sqrt(rx^2+ry^2+rz^2+rw^2);
                    r3 = r^3;
    
                    dBx(i,j,k) = dBx(i,j,k)+km*dy(l)*rz/r3;
                    dBy(i,j,k) = dBy(i,j,k)+km*dx(l)*rz/r3;
                    dBz(i,j,k) = dBz(i,j,k)+km*(dx(l)*ry-dy(l)*rx)/r3;
    
                end
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
        pcolor(x,y,(Bxz').^(1/3)); shading interp; colormap jet; colorbar;
        h1 = streamslice(x,z,Bx_xz',Bz_xz',3);
        set(h1,'Color',[0.8 1 0.9]);
        Bz=1;
   else
        idx_x = ceil(Lx/2);
        idx_y = ceil(Ly/2);
        Bz = squeeze(dBz(idx_x, idx_y,:));

        %Calcular gradiente del campo
        dBz_dz_profile = diff(Bz)./diff(z);
        z_mid=z(1:end-1) + diff(z)/2;
        
        %Graficar
        figure
        plot(z_mid,dBz_dz_profile, 'r-','LineWidth',2)
        xlabel('z'); ylabel('dBz/dz')
        title('Gradiente del campo magnético de Bz')
   end
end

function [Px,Py,Pz,dx,dy,dz] = espiras(nl,N,R,sz)
    dtheta = 2*pi/N;
    ang = 0:dtheta:(2*pi-dtheta);
    sz=1;
    s=1;
    for i = 1:nl
    Px(s:s+N-1) = R*cos(ang);
    Py(s:s+N-1) = R*sin(ang);
    Pz(s:s+N-1) = -nl/2*sz+(i-1)*sz;

    dx(s:s+N-1) = -Py(s:s+N-1)*dtheta;
    dy(s:s+N-1) = Px(s:s+N-1)*dtheta;
    
    s=s+N;
    end
    dz = zeros(1, N*nl);
    figure(1)
    quiver3(Px,Py,Pz,dx,dy,dz,0.5,'-r','LineWidth',2)
    view(-34,33)
    xlabel('x'); ylabel('y');zlabel('z')
    title('Corriente de espiras')
    axis equal
    
end
    function trayectoria(Bz,z,mag,m,zo,dt,vz,gamma)
    w=-m*9.81;
    zm(1)=zo;
    zmfree(1)=zo;
    vz(1)=0.7;
    vzfree(1)=0;
    tt(1)=0;
    cc =1;
   while zm(cc) > -3
        delta = 0.005;
        Bz_forward = interp1(z, Bz, zm(cc) + delta, 'linear', 'extrap');
        Bz_backward = interp1(z, Bz, zm(cc) - delta, 'linear', 'extrap');
        dBz_dz = (Bz_forward - Bz_backward) / (2 * delta);
    
        Fm(cc) = -mag * dBz_dz;
        Ff = -gamma * vz(cc);
        F(cc) = Fm(cc) + w + Ff;
        a = F(cc) / m;
    
        zm(cc+1) = zm(cc) + vz(cc)*dt + 0.5*a*dt^2;
        vz(cc+1) = vz(cc) + a*dt;
    
        zmfree(cc+1) = zmfree(cc) + vzfree(cc)*dt + 0.5*(-9.81)*dt^2;
        vzfree(cc+1) = vzfree(cc) + (-9.81)*dt;
    
        tt(cc+1) = tt(cc) + dt;
        cc = cc + 1;
    
        if abs(vz(cc)) < 1e-3
        break
        end
    end


    figure(99)
    hold on
    plot(tt,zm, 'r-', 'LineWidth',2)
    plot(tt,zmfree,'b--','LineWidth',2)

    grid on 
    xlabel('Time(s)')
    ylabel('Z position (m)')
    title(' Posicion vs tiempo en un Dipolo magnetico callendo a traves de una espira cargada ')
    legend('Trayectoria con la espira cargada', 'Trayectoria de caida libre','Southwest')
    axis([0 max(tt) -6 6])

    end

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

mag = 1000;     % Magnitud del momento magnético
m = 0.009;      % Masa (kg)
zo = 4.9;       % Posición inicial
dt = 0.05;      % Paso de tiempo (s)
vz = 0.7;       % Velocidad inicial (m/s)
gamma = 0.08;   % Coeficiente de fricción


trayectoria(Bz, z, mag, m, zo, dt, vz, gamma);
