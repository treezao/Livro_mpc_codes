%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% caso 1
caso1_s_warmstart

saidas1 = saidas;
entradas1 = entradas;
du1 = du;
t1 = t;
vx1 = vx;
tempo1 = tempo;

%% caso 2
caso2_c_warmstart

saidas2 = saidas;
entradas2 = entradas;
du2 = du;
t2 = t;
vx2 = vx;
tempo2 = tempo


%% tempos de cômputo
'mediana'
[median(tempo1(vx1)) median(tempo2(vx2))]*1000
1-ans(2)/ans(1)
'media'
[mean(tempo1(vx1)) mean(tempo2(vx2))]*1000
1-ans(2)/ans(1)


%% figuras
cores = gray(4);
cores = cores(1:end-1,:);

hf = figure;
h=subplot(2,1,1);
plot(t1,saidas1(vx1)+52,'LineWidth',tamlinha,'Color',cores(1,:));
hold on;
plot(t2,saidas2(vx2)+52,'-.','LineWidth',tamlinha,'Color',cores(2,:));

plot(t2,bandasup(vx2)+52,'--','LineWidth',tamlinha,'Color',cores(3,:));
plot(t2,bandainf(vx2)+52,'--','LineWidth',tamlinha,'Color',cores(3,:));



ylabel('Nível (%)', 'FontSize', tamletra);
% ylim([50 60]);
hl = legend('Caso 1','Caso 2','Bandas');

grid on;


set(h, 'FontSize', tamletra);

h=subplot(2,1,2);
plot(t1,entradas1(vx1)+30,'LineWidth',tamlinha,'Color',cores(1,:));
hold on
plot(t2,entradas2(vx2)+30,'-.','LineWidth',tamlinha,'Color',cores(2,:));

ylabel('Manipulada (%)', 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
% ylim([20 50]);
set(h, 'FontSize', tamletra);

grid on;


hf.Position = tamfigura;
hl.Position = [0.7488 0.4589 0.1714 0.1806];

%%% figura dos tempos

hf = figure;
plot(t1,tempo1(vx1)*1000,'LineWidth',tamlinha,'Color',cores(1,:))
hold on;
plot(t2,tempo2(vx2)*1000,'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 15]);
legend('Caso 1','Caso 2');

ylabel(['Tempo',' ','(milisegundos)'], 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
% ylim([20 50]);
set(h, 'FontSize', tamletra);

grid on

hf.Position = tamfigura;

