%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% caso 1
caso_com_rest_inversa

saidas1 = saidas;
entradas1 = entradas;
du1 = du;
t1 = t;
vx1 = vx;
xi1 = xi;


%% caso 2
caso_com_rest_inversa_mais_relax

saidas2 = saidas;
entradas2 = entradas;
du2 = du;
t2 = t;
vx2 = vx;
xi2 = xi;


%% figuras
cores = gray(4);
cores = cores(1:end-1,:);

hf = figure
h=subplot(2,1,1)
plot(t1,saidas1(vx1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,saidas2(vx2),'-.','LineWidth',tamlinha,'Color',cores(2,:))

plot(t2,refs(vx2),'--','LineWidth',tamlinha,'Color',cores(3,:))


ylabel('Controlada', 'FontSize', tamletra);
% ylim([50 60]);
hl = legend('Caso 1','Caso 2','Referência')
% hl.Position = [0.6685 0.6029 0.2375 0.1107];

grid on


set(h, 'FontSize', tamletra);

h=subplot(2,1,2)
plot(t1,entradas1(vx1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,entradas2(vx2),'-.','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada', 'FontSize', tamletra);
xlabel('Tempo (segundos)', 'FontSize', tamletra);
% ylim([20 50]);
set(h, 'FontSize', tamletra);

grid on

% h=subplot(3,1,3)
% plot(t1,erro(1:nit-N2-1),'LineWidth',tamlinha)
% title('Erro', 'FontSize', tamletra);
% xlabel('tempo', 'FontSize', tamletra);
% set(h, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.7077 0.5411 0.2054 0.1242];

%% figura das relaxações
hf = figure
plot(t,xi1(vx1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,max(xi2(vx2,:)'),'-.','LineWidth',tamlinha,'Color',cores(2,:))


ylabel('Valor Relaxação', 'FontSize', tamletra);
xlabel('Tempo (segundos)', 'FontSize', tamletra);
% set(hf, 'FontSize', tamletra);
hl = legend('Caso 1','Caso 2')

grid on

hf.Position = tamfigura;
