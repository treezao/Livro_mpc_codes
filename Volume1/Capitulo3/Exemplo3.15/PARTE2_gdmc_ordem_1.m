% clear all, 
% close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
global Kdmc1 Nf H N
%%
s = tf('s');
Ts = 1; % perído de amostragem em minutos

z = tf('z',Ts);
Gz = 0.1/(z-1.1)/z^2; % modelo discretizado
Gs = d2c(0.1/(z-1.1),'zoh')*exp(-2*Ts*s); % modelo contínuo


%% parâmetros de ajuste

N1 = 3; %horizonte de predição inicial
N2 = 25; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

Nss=80; % horizonte de modelo
Nf = 40; % horizonte de modelo filtrado
betaf = 0.8; % polo do filtro do gdmc
Gcoef = step(Gz,Ts:Ts:Nss*Ts);

tsim = 150*Ts;

%% montando as matrizes do DMC

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

%%% cálculo dos filtros dos erros de predição (SISO) para o GDMC

F = tf(0,1,Ts);
nf = 1;

pz = 1.1; % obtem o polo indesejado (instável) de malha aberta

for i=N1(1):N2(1)
    %%% monta sistema para obtenção dos parâmetros do filtro
    %%% primeira equação (z^i - F_i(z) = 0 -> para z=pz
    %%% segunda equação F_i(1) = 0 -> força ganho unitário
    indf = i-N1(1)+1;
    Af = [pz 1;
          1 1];
    bf = [pz^i*(pz-betaf);
          (1-betaf)];
    X = Af\bf;
    F(indf,1) = (X(1)*z+X(2))/(z-betaf);
    %%% teste da condição
%     pz^i-(X(1)*pz+X(2))/(pz-betaf)^2

    %%% armazena coeficientes gtil
    modDegrauUF{i} = filter(F(indf,1).num{1},F(indf,1).den{1},Gcoef);

end


%%% calcula a matriz H para o cálculo da resposta livre no caso GDMC
H1 = [];
H2 = [];

for i=N1(1):N2(1)
    H1 = [H1;Gcoef(i+1:i+Nf)'];
    H2 = [H2;modDegrauUF{i}(1:Nf)'];
    
end
H = H1-H2

%%
clear gdmc_simulink
sim_res = sim('PARTE2_sim_gdmc.slx')

%%
tout = sim_res.sim_output.time;
stout= size(tout,1)
np = 500;
if(stout>np)
    ind = 1:round(stout/np):stout;
end
tout1 = tout(ind);
ref = sim_res.sim_output.signals(1).values(ind);
yGDMC1 = sim_res.sim_output.signals(2).values(ind);

tin1 = sim_res.sim_in.time;
uGDMC1 = sim_res.sim_in.signals(1).values;


%% plota figura
cores = gray(3);
cores = cores(1:end-1,:);

hf = figure
h=subplot(2,1,1)
plot(tout1,yGDMC1,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tout1,ref,'--','LineWidth',tamlinha,'Color',cores(2,:))
ylim([0 1.6])
h.YTick = [0 0.5 1 1.5];
hl = legend('GDMC Ordem 1','Referência','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(tin1,uGDMC1,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
h.YTick = [-2 -1 0 1 2]
ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.6720 0.4895 0.2625 0.1242];




