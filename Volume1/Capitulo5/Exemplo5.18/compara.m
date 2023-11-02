%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% caso GPC
caso0_gpc

saidas0 = saidas;
entradas0 = entradas;
du0 = du;
t0 = t;
vx0 = vx;


%% caso GPC 1 - usar solução = 0
caso1_gpc

saidas1 = saidas;
entradas1 = entradas;
du1 = du;
t1 = t;
vx1 = vx;

%% caso GPC 2 - usar solução sub-ótima
caso2_gpc

saidas2 = saidas;
entradas2 = entradas;
du2 = du;
t2 = t;
vx2 = vx;

%% caso GPC 3 - usar solução sub-ótima
caso3_gpc

saidas3 = saidas;
entradas3 = entradas;
du3 = du;
t3 = t;
vx3 = vx;

%% caso GPC 4 - usar melhor solução das anteriores
caso4_gpc

saidas4 = saidas;
entradas4 = entradas;
du4 = du;
t4 = t;
vx4 = vx;

%% figuras
cores = gray(7);
cores = cores(1:end-1,:);

hf = figure
h=subplot(2,1,1)
plot(t0,saidas0(vx0),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,saidas1(vx1),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,saidas2(vx2),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t3,saidas3(vx3),'-','LineWidth',tamlinha,'Color',cores(4,:))
plot(t4,saidas4(vx4),'--','LineWidth',tamlinha,'Color',cores(5,:))

plot(t4,refs(vx4),'-.','LineWidth',tamlinha,'Color',cores(6,:))
ylim([0 1.2])
% h.YTick = [0 0.5 1 1.5];
hl = legend('GPC','GPC 1','GPC 2','GPC 3','GPC 4','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t0,entradas0(vx0),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,entradas1(vx1),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,entradas2(vx2),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t3,entradas3(vx3),'-','LineWidth',tamlinha,'Color',cores(4,:))
plot(t4,entradas4(vx4),'--','LineWidth',tamlinha,'Color',cores(5,:))


% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)

grid on

xlabel('Tempo (segundos)','FontSize', tamletra)

set(h, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.7470 0.3651 0.2054 0.3500];

% hf.Position = tamfigura3;
% hl.Position = [0.7238 0.6671 0.2054 0.0917];

ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))


print('tempo_insuficiente2_v2','-depsc')
