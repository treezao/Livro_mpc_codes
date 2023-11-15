clear all, 
close all, 
clc
close_system

run('../../../Bibliotecas/parametrosFiguras.m')
%%

Ts = 0.15;
z =tf('z',Ts);
s = tf('s');

%% caso 1- ajuste inicial
Ke = 1;
tau = 1.5;
G = Ke/(tau*s+1)
L1 = 10*Ts;
Lr = 12*Ts;

%%% planta real incerta
P = Ke/(1*s+1)*exp(-Lr*s);
Pz = c2d(P,Ts,'zoh');

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
Czps = c2d(Cps,Ts);

Gnz = c2d(G,Ts,'zoh')
Lnz = z^-(L1/Ts);

%%% critério de robustez
mf1 = feedback(Czps*Gnz,1);
mfpid = feedback(Cz*Gnz*Lnz,1);

%%% cenário de simulação
tsim = 30;
tref = 4;
ampref = 1;
tpert = 20;
amppert = 0.5;

% % 
% Czps = Cps;
% Gnz  = G;
% Lnz = exp(-L1*s);

out = sim('../Exemplo6.1/trocador_PS_PID_sim');

deltat = 5;
t1 = out.tout(1:deltat:end);
ref1 = out.simout(1:deltat:end,1);
y1 = out.simout(1:deltat:end,2:3);
u1 = out.simout(1:deltat:end,4:5);


%% caso 2 - ajuste robusto
Ke = 1;
tau = 1.5;
G = Ke/(tau*s+1)
L1 = 10*Ts;
Lr = 12*Ts;


%%% ajuste PID
alfa1=0.3;
T0=(sqrt(alfa1^2+alfa1)+alfa1)*L1/2;
kx=(2*tau)/((4*T0+L1)*Ke);

C = kx*(tau*s+1)*(L1/2*s+1)/(tau*s)/(alfa1*L1/2*s+1)
Cz = c2d(C,Ts,'tustin')


%%% ajuste PS - controlador primário
Ti = tau;
Kc = 2;

Cps = Kc*(1 + 1/Ti/s);
Czps = c2d(Cps,Ts);

Gnz = c2d(G,Ts,'zoh')
Lnz = z^-(L1/Ts);

%%% critério de robustez
mf2 = feedback(Czps*Gnz,1);
Pn = Gnz*Lnz;
dP = Pz/Pn - 1;

Czpseq = Czps/(1+Czps*Gnz*(1-Lnz));


%%% cenário de simulação
tsim = 30;
tref = 4;
ampref = 1;
tpert = 20;
amppert = 0.5;

out = sim('../Exemplo6.1/trocador_PS_PID_sim');


t2 = out.tout(1:deltat:end);
ref2 = out.simout(1:deltat:end,1);
y2 = out.simout(1:deltat:end,2:3);
u2 = out.simout(1:deltat:end,4:5);

%% dados dos critérios de robustez
w = logspace(-1,log10(2*pi*1/2/Ts),200);
[mag,pha] = bode([mf1,mf2,1/dP,mfpid],w);

mag1 = reshape(mag(1,1,:),200,1);
mag2 = reshape(mag(1,2,:),200,1);
mag3 = reshape(mag(1,3,:),200,1);
mag4 = reshape(mag(1,4,:),200,1);

w2 = logspace(-3,log10(2*pi*1/2/Ts),200);
[mag,pha] = bode([Cz,Czpseq],w2);
magcpid = reshape(mag(1,1,:),200,1);
magcps = reshape(mag(1,2,:),200,1);

phacpid = reshape(pha(1,1,:),200,1);
phacps = reshape(pha(1,2,:),200,1);


%%
cores = gray(4);
cores = cores(1:end-1,:);

%%% figura caso 1
hf = figure
h=subplot(2,1,1)
plot(t1,y1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,y1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,ref1,'-.','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 2])
hl = legend('PID','PS','Referência','Location','NorthWest')
hl.Position = [0.7351 0.5255 0.2375 0.1619];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t1,u1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,u1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)

hf.Position = tamfigura;
hl.Position = [0.6720 0.4928 0.2054 0.1806]; 
% print('trocador_PS_PID_erro_1','-depsc')

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

xlabel('Tempo (s)','FontSize',tamletra)


hl.Position = [0.6774 0.5831 0.2054 0.1806]; 
% print('trocador_PS_PID_erro_2','-depsc')


%%
cores = gray(5);
cores = cores(1:end-1,:);

hf = figure
h= subplot(1,1,1)
semilogx(w,20*log10(mag1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w,20*log10(mag2),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx(w,20*log10(mag4),':','LineWidth',tamlinha,'Color',cores(3,:))
semilogx(w,20*log10(mag3),'-.','LineWidth',tamlinha,'Color',cores(4,:))
semilogx(2*pi/2/Ts*[1 1],[-500 500],'LineWidth',tamlinha,'Color',cores(1,:))
grid on

hl = legend('PS(K_c=5)','PS(K_c=2)','PID','1/dP','Location','NorthEast')
xlabel('Frequência (rad/s)','FontSize', tamletra)
ylabel('Magnitude (dB)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


ylim([-20 20])
xlim([10^-1 3*10^1])

hf.Position = tamfigura;
hl.Position = [0.6935 0.6156 0.1929 0.2758]; 
% print('trocador_PS_PID_erro_robustez','-depsc')


%%
hf = figure
h= subplot(2,1,1)
semilogx(w2,20*log10(magcps),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w2,20*log10(magcpid),'--','LineWidth',tamlinha,'Color',cores(3,:))
semilogx(2*pi/2/Ts*[1 1],[-500 500],'LineWidth',tamlinha,'Color',cores(1,:))
grid on

hl = legend('PS(K_c=2)','PID','Location','NorthEast')
ylabel(['Magnitude',newline,'(dB)'],'FontSize', tamletra)

xlim([10^-2 3*10^1])
ylim([-5 40])

set(h, 'FontSize', tamletra);

h= subplot(2,1,2)
semilogx(w2,(phacps),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w2,(phacpid),'--','LineWidth',tamlinha,'Color',cores(3,:))
semilogx(2*pi/2/Ts*[1 1],[-500 500],'LineWidth',tamlinha,'Color',cores(1,:))
grid on

xlim([10^-2 3*10^1])
ylim([-90 45])
ylabel('Fase (graus)','FontSize', tamletra)
xlabel('Frequência (rad/s)','FontSize', tamletra)


set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.6935 0.7605 0.1929 0.1435]; 
% print('trocador_PS_PID_erro_ceq','-depsc')
