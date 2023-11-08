clear all
close all
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%

N_entradas = 2; % número de entradas da planta - manipuladas
N_saidas = 2; %número de saidas da planta - controladas
N_perturbacoes = 1;

s = tf('s');
% Definir os modelos que relacionam {saida x,entrada y};
G{1,1} = 1.43/(7.5*s+1)*exp(-s);
G{1,2} = -4.5/(13.4*s+1)*exp(-s);
G{2,1} = 0*s;
G{2,2} = 5.72/(14.7*s+1);

% Definir os modelos que relacionam {saida x,perturbação z};
Gp{1,1} = 2.1/(13*s+1)*exp(-s);
Gp{2,1} = -3.12/(14.7*s+1);

Ts = 1.0; %período de amostragem
Nrp = [30 50]; %horizonte de regime permanente - dimensão 1 x N_saidas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%respostas ao degrau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nrp = [100 100];
for i = 1:N_saidas
    for j = 1:N_entradas
        g{i,j} = step(G{i,j},0:Ts:(Nrp(i)-1)*Ts);
    end
    for j = 1:N_perturbacoes
        gp{i,j} = step(Gp{i,j},0:Ts:(Nrp(i)-1)*Ts);
    end
end

t = Ts*(0:max(Nrp)-1);


%%
cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h = subplot(2,1,1)
plot(t,g{1,1}(1:Nrp(1)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,g{1,2}(1:Nrp(1)),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,g{2,2}(1:Nrp(2)),'-.','LineWidth',tamlinha,'Color',cores(3,:))

hl1 = legend('y_1 x u_1','y_1 x u_2','y_2 x u_2','Location','SouthEast')

axis([0 max(Nrp)*Ts -5.0 8.0]);
% xlabel('Tempo (s)', 'FontSize', tamletra);
ylabel('Controladas', 'FontSize', tamletra);
set(h, 'FontSize', tamletra);
grid on
% print -depsc results_degrau_56.eps

h = subplot(2,1,2)
plot(t,gp{1,1}(1:Nrp(1)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,gp{2,1}(1:Nrp(1)),'--','LineWidth',tamlinha,'Color',cores(2,:))

hl2 = legend('y_1 x q','y_2 x q','Location','SouthEast')
axis([0 max(Nrp)*Ts -4.0 3.0]);
 %   ylabel('Controladas e Referências do Sistemas', 'FontSize', tamletra);
xlabel('Tempo (s)', 'FontSize', tamletra);
ylabel('Controladas', 'FontSize', tamletra);
set(h,'FontSize',tamletra);
grid on

hf.Position = tamfigura;
hl1.Position = [0.8113 0.5052 0.1654 0.2332];
hl2.Position = [0.8070 0.2181 0.1539 0.1592];


% print('tanque1_modelo','-depsc')

