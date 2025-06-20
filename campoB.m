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