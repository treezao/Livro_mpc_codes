function [X,j] = qp_gpm(H,b,R,r,epsilon, jmax)
    if(jmax < 0 )
        disp('Erro! jmax deve ser positivo');
        X = [];
        return
    end
    
    if(epsilon <0 )
        disp('Erro! epsf deve ser >0 ');
        X = [];
        return
    end
    
    H1 = inv(H);
    
    nrin = size(R,1);
    
    phi_j = 1*ones(nrin,1);
    phi_j1 = zeros(nrin,1);
    
    L = norm(R*H1*R','fro');
%     L = norm(R*H1*R');

    
    for j=1:round(jmax)
        X = -H1*(R'*phi_j+b);
        
        phi_j1 = phi_j;
        
        temp = phi_j + (1/L)*(R*X-r);
        phi_j = max(temp, zeros(nrin,1));
        
        if(norm(phi_j-phi_j1)<epsilon)
            break;
        end
        
    end

end