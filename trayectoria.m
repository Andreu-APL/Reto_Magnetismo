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