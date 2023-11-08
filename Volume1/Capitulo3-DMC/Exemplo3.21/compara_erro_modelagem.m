%%% Exemplo comparativo entre GDMC e DMC

clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% caso DMC
gdmc1_erro_modelagem

saidasGDMC1 = saidas;
entradasGDMC1 = entradas;
duGDMC1 = du;
tGDMC1 = t;
vxGDMC1 = vx;

mf1 = mf;


%% caso GDMC
gdmc2_erro_modelagem

saidasGDMC2 = saidas;
entradasGDMC2 = entradas;
duGDMC2 = du;
tGDMC2 = t;
vxGDMC2 = vx;

mf2 = mf;

%% dados para a análise da robustez
dP = Gzr/Gz - 1; %% erro de modelagem
w = logspace(-4,log10(2*pi/2/Ts),200); % faixa de frequencias

[magdP,phadP] = bode(dP,w);
[mag1,pha1] = bode(1/mf1,w);
[mag2,pha2] = bode(1/mf2,w);

%% figuras


cores = gray(4);
cores = cores(1:end-1,:);

%%% figura respostas
hf = figure
h=subplot(2,1,1)
plot(tGDMC1,saidasGDMC1(vxGDMC1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tGDMC2,saidasGDMC2(vxGDMC2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
plot(tGDMC2,refs(vxGDMC2),'--','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 1])
% h.YTick = [0 0.5 1 1.5];
hl = legend('GDMC1','GDMC1','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(tGDMC1,entradasGDMC1(vxGDMC1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(tGDMC2,entradasGDMC2(vxGDMC2),'-.','LineWidth',tamlinha,'Color',cores(2,:))
% h.YTick = [-2 -1 0 1 2]
ylim([-1.1 0.7])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7470 0.4659 0.2054 0.1806];


%%% figura robustez
hf = figure
semilogx(w,20*log10(magdP(:)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w,20*log10(mag1(:)),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx(w,20*log10(mag2(:)),'-.','LineWidth',tamlinha,'Color',cores(3,:))
semilogx([2*pi/Ts/2, 2*pi/Ts/2],[-100 100],'k') 

ylim([-20 15])
% h.YTick = [0 0.5 1 1.5];
hl = legend('\delta G(z)','GDMC1','GDMC2','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
% ylabel('Controlada','FontSize', tamletra)
grid on


ylabel('Magnitude (dB)','FontSize', tamletra)
grid on
xlabel('Frequência (rad/s)','FontSize', tamletra)

% set(gcf, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.1527 0.1769 0.2054 0.1806];

ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))

