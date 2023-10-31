clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
Ts = 1; % período de amostragem
A = [1 -1.1]; % denominador do modelo do processo
B = [0.1 0.2]; % numerador do modelo do processo
Bq = [0.3]; % numerador do modelo da perturbação do processo

d=1; % atraso em relação à manipulada
dq = 0; % atraso em relação à perturbação

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A
nbq = size(Bq,2)-1; % ordem do polinômio Bq



%% parâmetros de ajuste

N1 = 2; %horizonte de predição inicial
N2 = 4; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 2; % horizonte de controle
Nq = 1; % horizonte da ação antecipativa

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle


%% montagem das matrizes
Btil = conv(B,[zeros(1,d) 1]); % incorporação do atraso no numerador B
Bqtil = conv(Bq,[zeros(1,dq) 1]); % incorporação do atraso no numerador Bq


G = zeros(N,N); % matriz dinâmica G
H = zeros(N,nb+d); % matriz dos incrementos passados de controle 
Gq = zeros(N,N-1); % matris dos incrementos da perturbação futura
Hq = zeros(N,nbq+dq); % matriz dos incrementos passados da perturbação 

[E,F] = diofantina(conv(A,[1 -1]),N1,N2); % obtenção dos polinômios Ej, Fj

for i=N1:N2
    EjB = conv(E(i-N1+1,1:i),Btil);
    EjBq = conv(E(i-N1+1,1:i),Bqtil);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

    Gq(i-N1+1,1:i) = EjBq(i:-1:1);    
    Hq(i-N1+1,:) = EjBq(i+1:end);
    
end
G = G(:,1:Nu);
Gq = Gq(:,1:Nq);


G,F,H,Gq,Hq

%% obtenção do ganho irrestrito
Qu = lambda*eye(Nu); % matriz de ponderação dos incrementos de controle
Qe = delta*eye(N);  % matriz de ponderaão dos erros futuros

Kmpc = (G'*Qe*G+Qu)\G'*Qe;

Kmpc1 = Kmpc(1,:);

Kmpc1*H
Kmpc1*F
Kmpc1*Gq
sum(Kmpc1)


%% obtenção das funções de transferências

delta = tf([1 -1],[1 0],-1)
den1 = tf([1 Kmpc1*H],[1 zeros(1,size(Kmpc1*H,2))],-1)
C1num = tf(sum(Kmpc1),1,-1)
C2num = tf(Kmpc1*F,[1 zeros(1,size(Kmpc1*F,2)-1)],-1)

%%% controlador de realimentação
Cr = minreal(C2num/den1/delta)

%%% filtro de referência
Fr = minreal(C1num/C2num)

%%% controlador feed-forward
Ca = -tf([Kmpc1*Gq Kmpc1*Hq],1,-1)/den1



