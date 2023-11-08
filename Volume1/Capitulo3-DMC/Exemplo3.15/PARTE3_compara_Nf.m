%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% caso GDMC com Nf=40
PARTE3_gdmc_Nf_40

saidas1 = saidas;
entradas1 = entradas;
du1 = du;
t1 = t;
vx1 = vx;


%% caso GDMC Nf = 20
PARTE3_gdmc_Nf_20

saidas2 = saidas;
entradas2 = entradas;
du2 = du;
t2 = t;
vx2 = vx;

%% figuras


cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(t1,saidas1(vx1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,saidas2(vx2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,refs(vx1),'--','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 1.6])
h.YTick = [0 0.5 1 1.5];
hl = legend('N_f=40','N_f=20','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t1,entradas1(vx1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,entradas2(vx2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
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