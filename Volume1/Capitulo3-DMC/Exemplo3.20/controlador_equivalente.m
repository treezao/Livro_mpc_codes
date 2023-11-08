clear all
close all
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% modelo
modDegrau = [0.2,0.6,1.4,1.8,2,2,2,2,2,2,2]';

Ts = 1;
Gz = tf([0 0.2 0.4 0.8 0.4 0.2],[1],Ts,'Variable','z^-1'); % modelo discreto do sistema
na = size(Gz.den{1},2)-1;

figure
step(Gz)

%% PARAMETROS DE AJUSTE DO CONTROLADOR
N1 = 1; %horizonte de predição inicial
N2 = 3; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 2; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = .5; % ponderação do esforço de controle

Nss=5; % horizonte de modelo
Gcoef = step(Gz,Ts:Ts:2*Nss*Ts); % coeficientes da resposta ao degrau
Gqcoef = step(Gz,0:Ts:2*Nss*Ts); % coeficientes da resposta ao degrau da perturbação

%% Monta as matrizes offline que definem o problema DMC
G = zeros(N2,Nu);
G(:,1) = Gcoef(1:N2,1);


for i=2:Nu
    G(i:end,i) = G(1:end-(i-1),1);    
end

G = G(N1:end,:);

Qy = delta*eye(N2-N1+1);
Qu = lambda*eye(Nu);

Kdmc = inv(G'*Qy*G+Qu)*G'*Qy

Kdmc1 = Kdmc(1,:);

%%% calcula a matriz H para o cálculo da resposta livre no caso DMC
H1dmc = [];
H2dmc = [];

H1qdmc = [];
H2qdmc = [];

for i=N1(1):N2(1)
    H1dmc = [H1dmc;Gcoef(i+1:i+Nss)'];
    H2dmc = [H2dmc;Gcoef(1:Nss)'];

    H1qdmc = [H1qdmc;Gqcoef(i+1:i+Nss+1)'];
    H2qdmc = [H2qdmc;Gqcoef(1:Nss+1)'];    
    
end
H = H1dmc-H2dmc
Hq = H1qdmc-H2qdmc


%% cálculo do controlador equivalente
f0 = sum(Kdmc1);
Hu = conv([1 -1],[1 Kdmc1*H])
Huq = Kdmc1*Hq;

Cr = tf(f0,Hu,Ts,'Variable','z^-1') % controlador
Cff = -tf(Huq,[1 Kdmc1*H],Ts,'Variable','z^-1') % controlador antecipativo

mf = feedback(Cr*Gz,1)
