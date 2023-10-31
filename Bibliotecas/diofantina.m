%% Função para o cálculo da equação diofantina

function [E,F] = diofantina(A,N1,N2)

% Cálculo dos polinômios E(z) e F(z)

nA = size(A,2); % ordem de A(z)

% Cálculo de F(z)
f(1,:) = [1 zeros(1,nA-2)]; % inicializa F

for j=1:N2
    for i=1:nA-2
        f(j+1,i) = f(j,i+1)-f(j,1)*A(i+1);
    end
    f(j+1,nA-1) = -f(j,1)*A(nA);
end

F = f(1+N1:1+N2,:);

% Cálculo de E(z)
E = zeros(N2);
e(1) = 1;
E(1,1) = e(1);

for i=2:N2
    e(i) = f(i,1);
    E(i,1:i) = e;
end

E = E(N1:N2,:);