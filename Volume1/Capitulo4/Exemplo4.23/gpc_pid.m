clear all, 
% close all, 
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%
Ts = 1; % período de amostragem
A = [1 -0.8187];; % denominador do modelo do processo
B = [0.3625]; % numerador do modelo do processo

d=0; % atraso em relação à manipulada

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A



%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 15; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle


%% montagem das matrizes
Btil = conv(B,[zeros(1,d) 1]); % incorporação do atraso no numerador B


G = zeros(N,N); % matriz dinâmica G0
H = zeros(N,nb+d); % matriz dos incrementos passados de controle 

[E,F] = diofantina(conv(A,[1 -1]),N1,N2); % obtenção dos polinômios Ej, Fj

for i=N1:N2
    EjB = conv(E(i-N1+1,1:i),Btil);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

end
G = G(:,1:Nu);


G,F,H

%% obtenção do ganho irrestrito
Qu = lambda*eye(Nu); % matriz de ponderação dos incrementos de controle
Qe = delta*eye(N);  % matriz de ponderaão dos erros futuros

Kmpc = (G'*Qe*G+Qu)\G'*Qe;

Kmpc1 = Kmpc(1,:);

%% obtenção das funções de transferências

delta = tf([1 -1],[1 0],-1)
den1 = tf([1 Kmpc1*H],[1 zeros(1,size(Kmpc1*H,2))],-1)
C1num = tf(sum(Kmpc1),1,-1)
C2num = tf(Kmpc1*F,[1 zeros(1,size(Kmpc1*F,2)-1)],-1)

%%% controlador de realimentação
Cr = minreal(C2num/den1/delta)

%%% filtro de referência
Fr = minreal(C1num/C2num)

%%% definição da planta
Gz = tf(B,A,-1)

%%% FT de malha fechada para a referência
Hr = minreal(feedback(Gz*Cr,1)*Fr)

