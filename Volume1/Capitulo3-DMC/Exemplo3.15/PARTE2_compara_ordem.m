%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% caso GDMC com filtro de ordem 1
PARTE2_gdmc_ordem_1


%% caso GDMC com filtro de ordem 2
PARTE2_gdmc_ordem_2 

%% figuras
cores = gray(4);
cores = cores(1:end-1,:);

hf = figure
h=subplot(2,1,1)
plot(tout2,yGDMC2,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tout1,yGDMC1,'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(tout2,ref,'--','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 1.6])
h.YTick = [0 0.5 1 1.5];
hl = legend('Ordem 2','Ordem 1','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(tin2,uGDMC2,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tin1,uGDMC1,'-.','LineWidth',tamlinha,'Color',cores(2,:))
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
