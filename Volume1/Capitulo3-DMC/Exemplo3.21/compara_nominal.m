%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% caso DMC
gdmc1_nominal

saidasGDMC1 = saidas;
entradasGDMC1 = entradas;
duGDMC1 = du;
tGDMC1 = t;
vxGDMC1 = vx;


%% caso GDMC
gdmc2_nominal

saidasGDMC2 = saidas;
entradasGDMC2 = entradas;
duGDMC2 = du;
tGDMC2 = t;
vxGDMC2 = vx;

%% figuras


cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(tGDMC1,saidasGDMC1(vxGDMC1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tGDMC2,saidasGDMC2(vxGDMC2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(tGDMC2,refs(vxGDMC2),'--','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 1.6])
% h.YTick = [0 0.5 1 1.5];
hl = legend('GDMC1','GDMC2','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(tGDMC1,entradasGDMC1(vxGDMC1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tGDMC2,entradasGDMC2(vxGDMC2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7470 0.4659 0.2054 0.1806];

