%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% caso 1
deltaM = 1000; % ponderação das relaxações
PARTE2_caso1_restr_saida

saidas1 = saidas;
entradas1 = entradas;
du1 = du;
t1 = t;
vx1 = vx;


%% caso 2
deltaM = 100; % ponderação das relaxações
PARTE2_caso1_restr_saida

saidas2 = saidas;
entradas2 = entradas;
du2 = du;
t2 = t;
vx2 = vx;

%% caso 3
deltaM = 1; % ponderação das relaxações
PARTE2_caso1_restr_saida

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
plot(t3,saidas3(vx2)+52,':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t3,bandasup(vx3)+52,'--','LineWidth',tamlinha,'Color',cores(4,:))
plot(t3,bandainf(vx3)+52,'--','LineWidth',tamlinha,'Color',cores(4,:))



ylabel('Nível (%)', 'FontSize', tamletra);
% ylim([50 60]);
hl = legend('\psi_i = 1000','\psi_i = 100','\psi_i = 1','Bandas')

grid on


set(h, 'FontSize', tamletra);

h=subplot(2,1,2)
plot(t1,entradas1(vx1)+30,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,entradas2(vx2)+30,'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(t3,entradas3(vx3)+30,':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('Manipulada (%)', 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
% ylim([20 50]);
set(h, 'FontSize', tamletra);

grid on


hf.Position = tamfigura;
hl.Position = [0.7399 0.4895 0.1714 0.2371]
