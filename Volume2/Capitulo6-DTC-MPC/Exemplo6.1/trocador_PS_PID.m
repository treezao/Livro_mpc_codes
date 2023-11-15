clear all, 
close all, 
clc
close_system

run('../../../Bibliotecas/parametrosFiguras.m')
%%


Ts = 0.15;
z =tf('z',Ts);
s = tf('s');

%% caso 1 - atraso de 10 amostras
Ke = 1;
tau = 1.5;
G = Ke/(tau*s+1)
L1 = 10*Ts;

P = G*exp(-L1*s);

%%% ajuste PID
alfa1=0.3;
T0=(sqrt(alfa1^2+alfa1)+alfa1)*L1/2;
kx=(2*tau)/((4*T0+L1)*Ke);

C = kx*(tau*s+1)*(L1/2*s+1)/(tau*s)/(alfa1*L1/2*s+1)
Cz = c2d(C,Ts,'tustin')


%%% ajuste PS - controlador primário
Ti = tau;
Kc = 5;

Cps = Kc*(1 + 1/Ti/s);
Czps = c2d(Cps,Ts,'tustin');

Gnz = c2d(G,Ts,'zoh')
Lnz = z^-(L1/Ts);

%%% cenário de simulação
tsim = 30;
tref = 4;
ampref = 1;
tpert = 20;
amppert = 0.5;



out = sim('trocador_PS_PID_sim');

deltat = 5;
t1 = out.tout(1:deltat:end);
ref1 = out.simout(1:deltat:end,1);
y1 = out.simout(1:deltat:end,2:3);
u1 = out.simout(1:deltat:end,4:5);


%% caso 2 - atraso de 50 amostras
L2 = 50*Ts;

P = G*exp(-L2*s);

%%% ajuste PID
alfa1=0.3;
T0=(sqrt(alfa1^2+alfa1)+alfa1)*L2/2;
kx=(2*tau)/((4*T0+L2)*Ke);

C = kx*(tau*s+1)*(L2/2*s+1)/(tau*s)/(alfa1*L2/2*s+1)
Cz = c2d(C,Ts,'tustin')


%%% ajuste PS - controlador primário idêntico ao caso 1

Lnz = z^-(L2/Ts);

%%% cenário de simulação
tsim = 100;
tref = 4;
ampref = 1;
tpert = 50;
amppert = 0.5;



out = sim('trocador_PS_PID_sim');


t2 = out.tout(1:deltat:end);
ref2 = out.simout(1:deltat:end,1);
y2 = out.simout(1:deltat:end,2:3);
u2 = out.simout(1:deltat:end,4:5);

%%
cores = gray(4);
cores = cores(1:end-1,:);


%%% figura caso 1
hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t1,y1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,y1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,ref1,'-.','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 2])
hl = legend('PID','PS','Referência','Location','SouthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h.YTickLabel = trocaponto(h.YTickLabel)

h = subplot(2,1,2)
plot(t1,u1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,u1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)


hl.Position = [0.6810 0.6411 0.2054 0.1806]; 

% print('trocador_PS_PID_1','-depsc')


%%
%%% figura caso 2

hf = figure
hf.Position = tamfigura;

h=subplot(2,1,1)
plot(t2,y2(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,y2(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,ref2,'-.','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 2])
hl = legend('PID','PS','Referência','Location','SouthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h.YTickLabel = trocaponto(h.YTickLabel)

h = subplot(2,1,2)
plot(t2,u2(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,u2(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)

hl.Position = [0.6738 0.5863 0.2054 0.1806]; 
% print('trocador_PS_PID_2','-depsc')