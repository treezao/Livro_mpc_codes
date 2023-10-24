clear all,
close all,
clc

%% obtenção da lei de controle do dmc

g = [0.2 0.6 1.4 1.8 2 2 2 2 2 2 2 2 2 2 2 2 2]'; %%% resposta ao degrau do processo



Nss = 5; % horizonte de modelo
Nu = 2; % horizonte de controle
N = 3; % horizonte de predição
lambda = 0.5; % ponderação do esforço de controle
Ts = 1; % período de amostragem (em segundos)



%%% montagem da matriz dinâmica G
G = [];

for i=1:Nu
    G(i:N,i) = g(1:N-i+1);
end


%%% obtenção da matriz de ganhos K
x = (G'*G+eye(Nu)*lambda)
xi = inv(G'*G+eye(Nu)*lambda)
Kdmc = xi*G'




%%% obtenção da matriz H que multiplica os incrementos de controle passados
%%% e Hq que multiplica os incrementos passados da perturbação
H1 = [];
H2 = [];
for i=1:Nss
    H1(1:N,i) = -g(i);
    H2(1:N,i) = g(i+1:i+N);
end

H = H2+H1

Hq1(1:N,1) = 0;
Hq2(1:N,1) = g(1:N);
for i=1:Nss
    Hq1(1:N,i+1) = -g(i);
    Hq2(1:N,i+1) = g(i+1:i+N);
end

Hq = Hq2+Hq1

%%% obtenção dos coeficientes do controlador equivalente

KH = Kdmc(1,:)*H
KHq = Kdmc(1,:)*Hq %

f0 = sum(Kdmc(1,:))

%% controlador equivalente
C = tf(f0,conv([1 -1],[1 KH]),Ts,'variable','z^-1')

pole(C)



