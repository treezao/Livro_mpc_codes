clear all, 
close all, 
clc
close_system


run('../../../Bibliotecas/parametrosFiguras.m')
%%

Ts = 1;
z =tf('z',Ts);
s = tf('s')
%% caso 1- sem erro de modelagem
Kv = 0.4;
G = Kv/s
Gq = 0.2/s;
L1 = 4*Ts;

%%% planta
P = G*exp(-L1*s);
Pz = c2d(P,Ts,'zoh');

Pq = Gq*exp(-L1*s);
Pqz = c2d(Pq,Ts,'zoh');

d = L1/Ts;

Gnz = c2d(G,Ts,'zoh')
Lnz = z^-d;
Gqz = c2d(Gq,Ts,'zoh')*Lnz


t5d = 15
pd = 3/t5d
pdz = exp(-pd*Ts)

Kc = (-pdz+1)/0.4

Czps = Kc;

Kmf = (1-pdz)

a = 1- Kmf/((5*(1-pdz)+1))

Fe = (z-a)/(1-a)/z

Spsf = Gnz*(1-Lnz*Fe)
pole(Spsf)
zero(Spsf)

Spsf = minreal(Spsf)

mf = feedback(Czps*Gnz,1);

mfq1 = 1- mf*Lnz*Fe

%%% cenário de simulação
tsim = 80;
tref = 5;
ampref = 1;
tpert = 40;
amppert = 0.2;

out = sim('boiler_PSF_sim');

deltat = 10;
t1 = out.tout(1:deltat:end);
ref1 = out.simout(1:deltat:end,1);
y1 = out.simout(1:deltat:end,2);
u1 = out.simout(1:deltat:end,3);

%% caso 2- com erro de modelagem
L1r = 4*Ts*1.065;

%%% planta
P = G*exp(-L1r*s);
Pz = c2d(P,Ts,'zoh');

Pq = Gq*exp(-L1r*s);
Pqz = c2d(Pq,Ts,'zoh');

out = sim('boiler_PSF_sim');

deltat = 10;
t2 = out.tout(1:deltat:end);
ref2 = out.simout(1:deltat:end,1);
y2 = out.simout(1:deltat:end,2);
u2 = out.simout(1:deltat:end,3);




%%
cores = gray(4);
cores = cores(1:end-1,:);

fator = 10;
ind1 = 1:fator:size(y1,1);

%%% figura caso 1
hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t1(ind1),y1(ind1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1(ind1),ref1(ind1),'-.','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 2])
hl = legend('Controlada','Referência','Location','SouthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

% h.YTickLabel = trocaponto(h.YTickLabel)

h = subplot(2,1,2)
plot(t1(ind1),u1(ind1),'LineWidth',tamlinha,'Color',cores(1,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

xlabel('Tempo (minutos)','FontSize',tamletra)

% h.YTickLabel = trocaponto(h.YTickLabel)

% hl.Position = [0.6720 0.4928 0.2054 0.1806]; 

% print('boiler_PSF_1','-depsc')

%%
%%% figura caso 2
fator = 10;
ind2 = 1:fator:size(y2,1);

hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t2(ind2),y2(ind2),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2(ind2),ref2(ind2),'-.','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 2])
hl = legend('Controlada','Referência','Location','SouthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

% h.YTickLabel = trocaponto(h.YTickLabel)

h = subplot(2,1,2)
plot(t2(ind2),u2(ind2),'LineWidth',tamlinha,'Color',cores(1,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

xlabel('Tempo (minutos)','FontSize',tamletra)

% h.YTickLabel = trocaponto(h.YTickLabel)
% hl.Position = [0.6785 0.7018 0.2368 0.1136];

% print('boiler_PSF_2','-depsc')



