function [X,j] = qp_admm(H,b,R,r,ep,ed, jmax,rho)
    if(jmax < 0 )
        disp('Erro! jmax deve ser positivo');
        X = [];
        return
    end
    
    if(ep <0 )
        disp('Erro! epsf deve ser >0 ');
        X = [];
        return
    end
    
    if(ed <0 )
        disp('Erro! epso deve ser >0 ');
        X = [];
        return
    end
    
    if(rho <0 )
        disp('Erro! rho deve ser >0 ');
        X = [];
        return
    end
    
    nrin = size(R,1);

    z = zeros(nrin,1);
    phi = zeros(nrin,1);
    
    Htil = H+R'*R*rho;
    H1 = -inv(Htil);

    
    for j=1:round(jmax)
        z1 = z;
        
        btil = b+rho*R'*(phi-z-r);
        X = H1*btil;
        
        z = min(R*X-r+phi,zeros(nrin,1));
        phi = phi + (R*X-z-r);
        
        if(norm(R*X-z-r)<= ep && norm(rho*R'*(z-z1)) <= ed)
            return
        end
        
    end

end