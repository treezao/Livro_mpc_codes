%% Função para o cálculo da equação diofantina

function [E,F] = diofantinaC(A,C,N1,N2)

% Cálculo dos polinômios E(z) e F(z)
% delta = [1 -1]; % delta = 1-z^{-1}
% AD = conv(A,delta) ; % AD = A(z)*Delta(z)

nA1 = size(A,2); % numero de elemnetos de A(z)
nC = size(C,2); % numero de elementos de C(z)


%%% calculo da ordem máxima do polinômio F
nF1 = max(nC,N1-1+nA1)-N1;
nF2 = max(nC,N2-1+nA1)-N2;
nF = max(nF1,nF2);

% Cálculo de F(z)

if(nC>nA1)
    A = [A zeros(1,nC-nA1)];
    nA = size(A,2);

else
    nA = nA1;
end

f(1,:) = [C zeros(1,nA-nC)]; % inicializa F



for j=1:N2
    for i=1:nA-1
        f(j+1,i) = f(j,i+1)-f(j,1)*A(i+1);
    end
%     f(j+1,nA-1) = -f(j,1)*A(nA);
end

F = f(1+N1:1+N2,1:nF);



% Cálculo de E(z)
E = zeros(N2);
e(1) = C(1);
E(1,1) = e(1);

for i=2:N2
    e(i) = f(i,1);
    E(i,1:i) = e;
end

E = E(N1:N2,:);


 
