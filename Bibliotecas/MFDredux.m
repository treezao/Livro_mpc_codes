function [Ao,Bo] = MFDredux(A,B)
n = size(A,1);
m = size(A,2);
tol = 1e-6; % tolerância para indicar raízes iguais.

Ao = {};
Bo = {};
for i=1:n
    raizes = [];
    for j=1:m
        r = roots(A{i,j});
        while(length(r)>0)
            r0 = r(1); % raíz do polinômio
            mult = 0; % multiplicidade
            I = find(abs(r-r0) < tol); % encontrando os índices de raízes iguais com certa tolerância

            mult = size(I,1); % multiplicidade é o número de raízes encontradas
            
            r(I) = []; % removendo as raízes iguais processadas
            
            % verifica se raíz r0 já existia em outras FTs e se a
            % multiplicidade é maior ou menor
            if(length(raizes)>0)
                I = find(abs(raizes(:,1)-r0) < tol);

                if(size(I,1)==0) % se não existia
                    raizes = [raizes; r0 mult];
                else % se existia
                    if(raizes(I,2) < mult) % verifica multiplicidade
                        raizes(I,2) = mult;
                    end % se multiplicade já era maior, não faz nada
                end
            else
                raizes = [r0 mult];
            end
        end
    end
    
    Ao{i} = poly(repelem(raizes(:,1), raizes(:,2)));    
    
    
    for j=1:m
        [q,r] = deconv(Ao{i},A{i,j});
        
        Bo{i,j} = conv(B{i,j},q);
        
    end
end


end