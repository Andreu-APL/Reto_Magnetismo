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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        