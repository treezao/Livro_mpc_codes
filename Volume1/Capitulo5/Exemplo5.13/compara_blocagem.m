%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% caso DMC sem blocagem
caso1_sem_blocagem

saidas1 = saidas;
entradas1 = entradas;
du1 = du;
t1 = t;
vx1 = vx;


%% caso DMC com blocagem uniforme
caso2_bloc_unif

saidas2 = saidas;
entradas2 = entradas;
du2 = du;
t2 = t;
vx2 = vx;

%% caso DMC com blocagem exponencial
caso3_bloc_exp

saidas3 = saidas;
entradas3 = entradas;
du3 = du;
t3 = t;
vx3 = vx;

%% figuras


cores = gray(5);
cores = cores(1:end-1,:);


hf = figure
h=subplot(3,1,1)
plot(t1,saidas1(vx1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,saidas2(vx2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(t3,saidas3(vx3),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t1,refs(vx1),'--','LineWidth',tamlinha,'Color',cores(4,:))

hl = legend('S/ Bloc.','Bloc. Unif.','Bloc. Exp.','Referência','Location','SouthEast')
ylabel(['Temperatura',newline,'(^oC)'], 'FontSize', tamletra);
% axis([0 60 0 350]);
grid on
set(h, 'FontSize', tamletra);

h = subplot(3,1,2)
plot(t1,entradas1(vx1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,entradas2(vx2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(t3,entradas3(vx3),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel(['Potência',newline,'(kW)'], 'FontSize', tamletra);
% axis([0 60 0 450]);
grid on
set(h, 'FontSize', tamletra);

h=subplot(3,1,3)
plot(t1,du1(vx1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,du2(vx2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(t3,du3(vx3),':','LineWidth',tamlinha,'Color',cores(3,:))


ylabel(['\Delta u'], 'FontSize', tamletra);
ylim([-2.5 2.5])
xlabel('Tempo (minutos)', 'FontSize', tamletra);
% axis([0 60 0 450]);
grid on
set(h, 'FontSize', tamletra);


hf.Position = [680 558 560 420];
hl.Position = [0.7148 0.6514 0.2054 0.1750];

ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))