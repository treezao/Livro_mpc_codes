%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% caso DMC
exemplo_dmc

saidasDMC = saidas;
entradasDMC = entradas;
duDMC = du;
tDMC = t;
vxDMC = vx;


%% caso GDMC
exemplo_gdmc

saidasGDMC = saidas;
entradasGDMC = entradas;
duGDMC = du;
tGDMC = t;
vxGDMC = vx;

%% figuras


cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(tGDMC,saidasGDMC(vxGDMC),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tDMC,saidasDMC(vxDMC),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(tGDMC,refs(vxGDMC),'--','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 1.6])
h.YTick = [0 0.5 1 1.5];
hl = legend('GDMC','DMC','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(tGDMC,entradasGDMC(vxGDMC),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tDMC,entradasDMC(vxDMC),'-.','LineWidth',tamlinha,'Color',cores(2,:))
h.YTick = [-2 -1 0 1 2]
ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7470 0.4659 0.2054 0.1806];

ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))