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


%% caso GDMC1
caso_gdmc1

saidasGDMC = saidas;
entradasGDMC = entradas;
duGDMC = du;
tGDMC = t;
vxGDMC = vx;

%% caso GDMC2
caso_gdmc2

saidasGDMC2 = saidas;
entradasGDMC2 = entradas;
duGDMC2 = du;
tGDMC2 = t;
vxGDMC2 = vx;

%% figuras


cores = gray(5);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(tDMC,saidasDMC(vxDMC),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tGDMC,saidasGDMC(vxGDMC),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(tGDMC2,saidasGDMC2(vxGDMC2),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(tGDMC,refs(vxGDMC),'--','LineWidth',tamlinha,'Color',cores(1,:))
ylim([0 1.5])
h.YTick = [0 0.5 1 1.5];
hl = legend('DMC','GDMC1','GDMC2','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(tDMC,entradasDMC(vxDMC),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tGDMC,entradasGDMC(vxGDMC),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(tGDMC2,entradasGDMC2(vxGDMC2),':','LineWidth',tamlinha,'Color',cores(3,:))

ylim([-1 2])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7167 0.4605 0.2054 0.2371];

