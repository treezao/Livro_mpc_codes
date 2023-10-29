%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% caso DMC
estudo1_dmc

saidasDMC = saidas;
entradasDMC = entradas;
duDMC = du;
tDMC = t;
vxDMC = vx;


%% caso GDMC
estudo1_gdmc

saidasGDMC = saidas;
entradasGDMC = entradas;
duGDMC = du;
tGDMC = t;
vxGDMC = vx;

%% figuras


cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(3,1,1)
plot(tGDMC,saidasGDMC(vxGDMC),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tDMC,saidasDMC(vxDMC),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(tGDMC,refs(vxGDMC),'--','LineWidth',tamlinha,'Color',cores(3,:))
ylim([1.3 1.45])
hl = legend('GDMC','DMC','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('C_o (mol/L)','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(tGDMC,entradasGDMC(vxGDMC),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tDMC,entradasDMC(vxDMC),'-.','LineWidth',tamlinha,'Color',cores(2,:))
h.YTick = [3.2400 3.2800 3.3200]

ylabel('C_f (mol/L)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

h = subplot(3,1,3)
plot(tGDMC,duGDMC(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tDMC,duDMC(vx),'-.','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (s)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);



hf.Position = tamfigura;
hl.Position = [0.7785 0.5637 0.2054 0.1789];

ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))