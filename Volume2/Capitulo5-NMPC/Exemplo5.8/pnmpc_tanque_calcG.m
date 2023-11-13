clear all, close all, clc

%% parametros da simulação
nit = 100; % numero de iterações

h0 = 50;
u0 = 50;
q0 = 20;

global K1;
global K2;
global g;
K1 = 2;
K2 = 0.1;
g = 9.81;

N1 = 1;
N2 = 3;
Nu = N2;

%%% calculo de y0
y0 = zeros(N2,1);
y0(1) = modeloTanque(h0,u0,q0);
for i=2:N2
    y0(i) = modeloTanque(y0(i-1),u0,q0);
end

%%% calculo de yi para obter G
epsilon = u0/1000;
y = zeros(N2,N2);
for j=1:N2
    u = h0*ones(N2,1);
    u(j:end) = u(j:end)+epsilon;
    
    y(1,j) = modeloTanque(h0,u(1),q0);
    for i=2:N2
        y(i,j) = modeloTanque(y(i-1,j),u(i),q0);
    end
end

%%% cálculo de G
G = (y-repmat(y0,1,N2))./epsilon


%%
figure
plot(G(:,1))

%% função do modelo do sistema de tanque
function hk = modeloTanque(hk1,ak1,qk1)
    
    %%% parâmetros do modelo
    global K1;
    global K2;
    global g;
    %%% modelo discreto do tanque
    hk = hk1 + K1*ak1 -K2*qk1*sqrt(2*g*hk1);
    
    if(hk<0)
        hk = 0;
    end
end


