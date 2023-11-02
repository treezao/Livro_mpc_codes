function [X,it] = qp_barreira(H,b,R,r,Lnt,epsilon,mu, beta,w)
    it = 0;
    if(beta > 1 || beta < 0 )
        disp('Erro! Beta deve ser entre 0 e 1');
        X = [];
        return
    end
    
    if(w <0.01 || w > 0.5)
        disp('Erro! w deve ser entre 0 e 0.5');
        X = [];
        return
    end
    
    if(Lnt < 0 )
        disp('Erro! Lnt deve ser positivo');
        X = [];
        return
    end
    
    if(epsilon <0 )
        disp('Erro! epsilon deve ser >0 ');
        X = [];
        return
    end
    
    if( mu < 1)
        disp('Erro! mu deve ser > 1');
        X = [];
        return
    end
    
    if (size(H,1) ~= size(b,1))
        disp('Erro! Matrizes H e b com dimensões incompatíveis')
        X = [];
        return
    end
    
    if (size(H,1) ~= size(R,2))
        disp('Erro! Matrizes H e R com dimensões incompatíveis')
        X = [];
        return
    end
    
    if (size(R,1) ~= size(r,1))
        disp('Erro! Matrizes R e r com dimensões incompatíveis')
        X = [];
        return
    end  
    kb = mu*epsilon;
    Lnt = round(Lnt);
    
    J = @(du) (1/2)*du'*H*du + b'*du;
    Eta = @(du,kb) sum(-(1/kb)*log(-(R*du-r)));
    
    nrin = size(R,1);      
    
    X = zeros(size(H,1),1);
    
    xnt = Inf;
    
    while nrin/kb > epsilon
        i = 0;
        while i <= Lnt && norm(xnt)>epsilon
            Go = H*X+b;
            Ho = H;
            for j=1:nrin
                Go = Go - (1/kb)/(R(j,:)*X-r(j))*R(j,:)';
                Ho = Ho + (1/kb)/(R(j,:)*X-r(j))^2*(R(j,:)'*R(j,:));
            end

            xnt = -Ho\Go;
            
            l = 1;
            while max((R*(X+l*xnt)-r))>=0
                l = beta*l;
            end
            
            Jdu = J(X);
            Etadu = Eta(X,kb);
            
            while ( J(X+l*xnt)+Eta(X+l*xnt,kb) > (Jdu +Etadu + l*w*Go'*xnt) )
                l = beta * l;
            end
            
            X = X + l*xnt;
            
            i = i+1;
        end
        
        kb = mu * kb;
        it = it +1;
    end
    
    it = it*round(Lnt);

end