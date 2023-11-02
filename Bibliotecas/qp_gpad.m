function [X,j] = qp_gpad(H,b,R,r,epsf,epso, jmax)
    if(jmax < 0 )
        disp('Erro! jmax deve ser positivo');
        X = [];
        return
    end
    
    if(epsf <0 )
        disp('Erro! epsf deve ser >0 ');
        X = [];
        return
    end
    
    if(epso <0 )
        disp('Erro! epso deve ser >0 ');
        X = [];
        return
    end
    
    J = @(du) (1/2)*du'*H*du + b'*du;
    
    H1 = -inv(H);
    
    nrin = size(R,1);
    
    phi_j = zeros(nrin,1);
    phi_j1 = zeros(nrin,1);
    
%     L = sqrt(max(eig(H'*H)));
%     L = norm(H,'fro');

    L = norm(R*H1*R','fro');

    
    for j=1:round(jmax)
        if(j==1)
            beta = 0;
        else
            beta = (j-1)/(j+2);
        end
        
        a = phi_j + beta*(phi_j - phi_j1);
        
        X = H1*(R'*a+b);
        s = (R*X-r);
        
        temp = sum(s >= epsf);
        
        if(~temp)
            if(-a'*s < epso)
                return
            end
        end
        
        phi_j1 = phi_j;
        
        phi_j = max(a+(1/L)*s, zeros(nrin,1));
        
    end

end