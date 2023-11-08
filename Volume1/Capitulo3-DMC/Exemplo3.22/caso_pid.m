% clear all, 
% close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');
Ts = 1; % perído de amostragem em minutos

z = tf('z',Ts);
Ts = 1;
Gz = tf([0 0.2 0.4 0.8 0.4 0.2],[1],Ts,'Variable','z^-1'); % modelo discreto do sistema
num = Gz.num{1}; % numerador do modelo discreto
den = Gz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Gz.inputdelay; % atraso discreto


%% controlador PID
Cpid = tf([0.7457, -0.8722, +0.305],conv([1 -1],[1 -0.1624]),Ts)
Fpid = tf(0.65582*[1 -0.9 0.265],[1 -1.17 0.409],Ts)
Fpid = Fpid/dcgain(Fpid)

Cffpid = -tf([0.7251,0.5322,0.4166],[1,0.3706,0.4397],Ts)

mfpid = feedback(Cpid*Gz,1)

%% inicialização vetores
%%% inicializacao dos estados e variaveis da simulação
nin = 5+1;
nit = 40 + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+round(25/Ts):end) = 2;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+round(4/Ts):end) = 1;

rf = zeros(nit,1);
erro = zeros(nit,1);
entradas = zeros(nit,1); % vetor do sinal de controle
du = zeros(nit,1); % vetor dos incrementos de controle
ua = zeros(nit,1);
uf = zeros(nit,1);

for k = nin:nit
    %% -- simulador
    saidas(k) = -den(2:end)*saidas(k-1:-1:k-na) + num*(entradas(k-dd:-1:k-nb-dd-1) + perts(k-dd:-1:k-dd-nb-1));
    
    %% -- Controlador PID
    rf(k) = -Fpid.den{1}(2:end)*rf(k-1:-1:k-2) + Fpid.num{1}*refs(k:-1:k-2);
    
    erro(k) = rf(k) - saidas(k);
    
    %%% realimentacao
    uf(k) = -Cpid.den{1}(2:end)*uf(k-1:-1:k-2)+ Cpid.num{1}*erro(k:-1:k-2);
    
    %%% antecipativo
    ua(k) = -Cffpid.den{1}(2:end)*ua(k-1:-1:k-2)+ Cffpid.num{1}*perts(k:-1:k-2);

    entradas(k) = ua(k) + uf(k);
    du(k) = entradas(k)-entradas(k-1);
end

%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;

cores = gray(3);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(t,saidas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx),'--','LineWidth',tamlinha,'Color',cores(2,:))
% ylim([0 1.6])
% h.YTick = [0 0.5 1 1.5];
hl = legend('PID','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7202 0.4960 0.2054 0.1242]




