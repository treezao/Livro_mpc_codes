%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% caso 1
caso1_pond_variavel

saidas1 = saidas;
entradas1 = entradas;
du1 = du;
t1 = t;
vx1 = vx;


%% caso 2
caso2_pond_fixa

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
plot(t1,saidas1(vx1)+52,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,saidas2(vx2)+52,'-.','LineWidth',tamlinha,'Color',cores(2,:))

plot(t2,refs(vx2)+52,':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t2,bandasup(vx2)+52,'--','LineWidth',tamlinha,'Color',cores(3,:))
plot(t2,bandainf(vx2)+52,'--','LineWidth',tamlinha,'Color',cores(3,:))



ylabel('Nível (%)', 'FontSize', tamletra);
xlim([0 26]);
ylim([48 61])
hl = legend('Caso 1','Caso 2','Referência','Limites')

grid on


set(h, 'FontSize', tamletra);

h=subplot(2,1,2)
plot(t1,entradas1(vx1)+30,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,entradas2(vx2)+30,'-.','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada (%)', 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
xlim([0 26]);
ylim([14 34])
set(h, 'FontSize', tamletra);

grid on


hf.Position = tamfigura;
hl.Position = [0.7077 0.4588 0.2054 0.2371];