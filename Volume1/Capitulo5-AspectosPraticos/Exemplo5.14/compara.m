%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% caso 1
caso1_sem_u_target

saidas1 = saidas;
entradas1 = entradas;
du1 = du;
t1 = t;
vx1 = vx;


%% caso 2
caso2_u_target_04

saidas2 = saidas;
entradas2 = entradas;
du2 = du;
t2 = t;
vx2 = vx;

%% caso 3
caso3_u_target_5

saidas3 = saidas;
entradas3 = entradas;
du3 = du;
t3 = t;
vx3 = vx;

%% figuras
cores = gray(5);
cores = cores(1:end-1,:);

hf = figure
h=subplot(2,1,1)
plot(t1,saidas1(vx1)+52,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,saidas2(vx2)+52,'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,saidas3(vx2)+52,':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t2,refs(vx2)+52,'--','LineWidth',tamlinha,'Color',cores(4,:))




ylabel('Nível (%)', 'FontSize', tamletra);
ylim([45 56]);
hl = legend('Caso 1','Caso 2','Caso 3','Referência')

grid on


set(h, 'FontSize', tamletra);

h=subplot(2,1,2)
plot(t1,entradas1(vx1)+30,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,entradas2(vx2)+30,'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(t3,entradas3(vx3)+30,':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t2,refsu(vx2)+30,'--','LineWidth',tamlinha,'Color',cores(4,:))

ylabel('Manipulada (%)', 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
ylim([29 35]);
set(h, 'FontSize', tamletra);

grid on


hf.Position = tamfigura;
hl.Position = [0.1448 0.4987 0.2054 0.2348];