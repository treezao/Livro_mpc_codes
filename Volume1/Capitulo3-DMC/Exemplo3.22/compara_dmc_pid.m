%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% caso DMC
caso_dmc

saidasDMC = saidas;
entradasDMC = entradas;
duDMC = du;
tDMC = t;
vxDMC = vx;

Cr
Cff


%% caso GDMC
caso_pid

saidasPID = saidas;
entradasPID = entradas;
duPID = du;
tPID = t;
vxPID = vx;

Cpid
Fpid
Cffpid

%% dados frequenciais dos controladores
[mag1,pha1,w1] = bode(Cr);
[mag2,pha2] = bode(Cpid,w1);
[mag3,pha3,w2] = bode(Cff);
[mag4,pha4] = bode(Cffpid,w2);


%% figuras
cores = gray(4);
cores = cores(1:end-1,:);


%%% figura das respostas
hf = figure
h=subplot(2,1,1)
plot(tDMC,saidasDMC(vxDMC),'-.','LineWidth',tamlinha,'Color',cores(2,:))
hold on
plot(tPID,saidasPID(vxPID),'LineWidth',tamlinha,'Color',cores(1,:))
plot(tPID,refs(vxPID),'--','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 1.6])
% h.YTick = [0 0.5 1 1.5];
hl = legend('DMC','PID','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(tDMC,entradasDMC(vxDMC),'-.','LineWidth',tamlinha,'Color',cores(2,:))
hold on
plot(tPID,entradasPID(vxPID),'LineWidth',tamlinha,'Color',cores(1,:))

% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7470 0.4659 0.2054 0.1806];




%%% figuras das respostas em frequencia
hf = figure
h = subplot(2,1,1)
semilogx(w1,20*log10(mag1(:)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w1,20*log10(mag2(:)),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx([2*pi/Ts/2, 2*pi/Ts/2],[-100 100],'k') 
set(h,'Fontsize',tamletra)
xlim([min(w1) max(w1)])
hl = legend('DMC','PID','Location','NorthEast')
grid on

title('Comparação entre os controladores DMC e PID')
ylabel('Magnitude (dB)','FontSize', tamletra)
grid on

h = subplot(2,1,2)
semilogx(w1,pha1(:),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w1,pha2(:),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx([2*pi/Ts/2, 2*pi/Ts/2],[-100 100],'k') 
grid on
xlim([min(w1) max(w1)])
ylabel('Fase (graus)', 'FontSize',tamletra)

xlabel('Frequência (rad/s)','FontSize', tamletra)
set(h,'Fontsize',tamletra)

hf.Position = tamfigura;

ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))






hf = figure
h = subplot(2,1,1)
semilogx(w2,20*log10(mag3(:)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w2,20*log10(mag4(:)),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx([2*pi/Ts/2, 2*pi/Ts/2],[-100 100],'k') 

xlim([10^-1 max(w1)])
hl = legend('DMC','PID','Location','SouthWest')
grid on

title('Comparação entre os controladores antecipativos DMC e PID')
ylabel('Magnitude (dB)','FontSize', tamletra)
grid on
set(h,'Fontsize',tamletra)

h = subplot(2,1,2)
semilogx(w2,pha3(:),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w2,pha4(:),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx([2*pi/Ts/2, 2*pi/Ts/2],[-100 100],'k') 
grid on
xlim([10^-1 max(w1)])
ylabel('Fase (graus)', 'FontSize',tamletra)

xlabel('Frequência (rad/s)','FontSize', tamletra)

set(h,'Fontsize',tamletra)

hf.Position = tamfigura;

