% para resposta ao degrau trocador de calor
close all
clear all

run('../../../Bibliotecas/parametrosFiguras.m')

%% Resposta do sistema sem atraso (Exemplo 3.1)
G=tf(2,[2 3 1]);

[g,Tg]=step(G*10); % aplicação do degrau com amplitude 10

y=g+70; %% conversão dos dados para unidades adequadas



hf = figure
h=subplot(1,1,1);
plot(Tg',y','LineWidth',tamlinha)
grid on
ylabel('Temperatura (^oC)', 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
set(h, 'FontSize', tamletra);

hf.Position = tamfigura;

%% resposta dos sistema com atraso (Exemplo 3.2)

G=tf(2,[2 3 1]);

set(G,'InputDelay',1)

[g,Tg]=step(G*10); % aplicação do degrau com amplitude 10

y=g+70; %% conversão dos dados para unidades adequadas

hf2 = figure

h=subplot(1,1,1);
plot(Tg',y','LineWidth',tamlinha)
grid on
ylabel('Temperatura (^oC)', 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
set(h, 'FontSize', tamletra);

hf2.Position = tamfigura;
