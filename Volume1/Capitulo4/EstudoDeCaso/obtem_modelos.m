%%% Estudo de Caso GPC - obtenção dos modelos
clear all
close all
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%% carrega o ponto de operação
pontoOperacao


%% obtendo a resposta ao degrau para a entrada
%%% aplicação de uma variação de 10% em ae
nit = 500
y = h0+zeros(1,nit);
u = ae0+zeros(1,nit);
q = as0+zeros(1,nit);

u(10:end) = ae0*1.1;

for i=2:nit
    y(i) = simModelo(u(i-1),q(i-1),y(i-1));
end

%% obtendo a resposta ao degrau para a perturbação
%%% aplicação de uma variação de 10% em as
nit = 500
y2 = h0+zeros(1,nit);
u2 = ae0+zeros(1,nit);
q2 = as0+zeros(1,nit);

q2(10:end) = as0*1.1;

for i=2:nit
    y2(i) = simModelo(u2(i-1),q2(i-1),y2(i-1));
end

%% obtenção da resposta ao degrau
%%% a partir da resposta do sistema, dividir pela variação da entrada.
Gstep = (y(11:end)-h0)'/(1.1*ae0-ae0); % modelo para ae
Gstepq = (y2(11:end)-h0)'/(1.1*as0-as0); % modelo para as
save('modeloDegrau','Gstep','Gstepq')


%% obtenção da função de transferência
%%% aproximação por FT de primeira ordem

%%% modelo para ae
Ke = (y(end)-h0)/(u(end)-ae0)

y63 = (y(end)-h0)*0.63 + h0 % valor de 63%
t63 = 85-10 % tempo de 63%

pd = (1-0.63)^(1/t63) % polo discreto

Gz = tf(Ke*(1-pd),[1 -pd],-1)
figure
step(Gz)

%%% modelo para as
%%% analisando o modelo, só muda o ganho
Kes = (y2(end)-h0)/(q2(end)-as0)

y63_2 = (y2(end)-h0)*0.63 + h0 % valor de 63%
t63_2 = 71-10 % tempo de 63%

pd2 = (1-0.63)^(1/t63_2) % polo discreto



Gzq = tf(Kes*(1-pd2),[1 -pd2],-1)
figure
step(Gzq)

save('modeloFT','Gz','Gzq')


%% Geração das figuras
ampdegrau = 0.05;

[s1,t1] = step(ampdegrau*Gz,490);
[s2,t2] = step(ampdegrau*Gzq,490);

cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(t1,s1,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(ampdegrau*Gstep,'-.','LineWidth',tamlinha,'Color',cores(2,:))
ylim([0 1])
hl = legend('Mod. FT','Mod. Degrau','Location','SouthEast')
% hl.Position = [0.6921    0.6355    0.2439    0.1662];
title(['\fontsize{' num2str(tamletra) '}' 'Resposta ao degrau para a entrada']) 
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h=subplot(2,1,2)
plot(t2,s2,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(ampdegrau*Gstepq,'-.','LineWidth',tamlinha,'Color',cores(2,:))
title(['\fontsize{' num2str(tamletra) '}' 'Resposta ao degrau para a perturbação']) 
ylabel('Controlada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)
set(h, 'FontSize', tamletra);

hf.Position = tamfigura;


