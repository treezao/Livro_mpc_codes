%%% Exemplo duplo integrador do capítulo GPC SISO
%%% cálculo do ganho irrestrito
%%% cálculo das funções de transferência equivalentes.

clear all
close all
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% Modelo da planta
b = 0.4;
a = 0.8;

z = tf('z',1);
Gz = b/(z-a)

dr=1; % atraso real

nb = 0; % ordem do polinômio B
na = 1; % ordem do polinômio A

% f = 0.95;
f = 1;
Gzr = (1-a*f)/(1-a)*b/(z-a*f)


%% Sintonia do GPC
Nu = 2; %% horizonte de controle

lambda = .2; % ponderação do incremento de controle
delta = 1;  % ponderação dos erros futuros



%% Ajuste do GPC 1 -> d = 1
d = 1;
N1 = d+1; %% horizonte de predição inicial
N2 = d+10; %% horizonte predição
Ny = N2-N1+1;

Btil = conv(Gz.num{1}(2:end),[zeros(1,d) 1]); % incorporação do atraso no numerador B


G = zeros(Ny,Ny); % matriz dinâmica G
H = zeros(Ny,nb+d);

[E,F] = diofantina(conv(Gz.den{1},[1 -1]),N1,N2);

for i=N1:N2
    EjB = conv(E(i-N1+1,1:i),Btil);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

    
end
G = G(:,1:Nu);
G,F,H

%% obtenção do ganho irrestrito
Qu = lambda*eye(Nu); % matriz de ponderação dos incrementos de controle
Qe = delta*eye(Ny);  % matriz de ponderaão dos erros futuros

Kmpc = (G'*Qe*G+Qu)\G'*Qe;

Kmpc1 = Kmpc(1,:);

%% vetores de simulação
nin = 50; % iteração inicial (para inicialização correta)
nit = nin+30; % número de iterações da simulação

refs = zeros(1,nit); % vetor das referências
saidas = zeros(1,nit); % vetor da saída
perts = zeros(1,nit); % vetor da perturbação
entradas = zeros(1,nit); % vetor do sinal de controle
du = zeros(1,nit); % vetor dos incrementos de controle

refs(nin+5:end) = 1;
% perts(1,nin+50:end) = 0.5;

for i=nin+1:nit
    %%% simulação do processo
    saidas(1,i) = -Gzr.den{1}(2:end)*saidas(1,i-1:-1:i-na)' ...
                  +Gzr.num{1}(2:end)*entradas(1,i-dr-1:-1:i-1-nb-dr)';
              
    %%% vetor de referências (caso futuro desconhecido)
    R = refs(i)*ones(Ny,1);
    
    %%% cálculo da resposta livre;
    f = F*saidas(1,i:-1:i-na)';
    
    
    if(~isempty(H))
        f = f+ H*du(1,i-1:-1:i-nb-d)'; % parcela dos incrementos de controle
    end
    
    %%% cálculo do incremento de controle ótimo
    du(1,i) = Kmpc1*(R-f);
    
    %%% cálculo do sinal de controle ótimo
    entradas(1,i) = du(1,i)+entradas(1,i-1);
    
              
end

%%
saidas1 = saidas;
entradas1 = entradas;
du1 = du;

% figure
% subplot(2,1,1)
% plot(saidas1)
% subplot(2,1,2)
% plot(entradas1)



%% Ajuste do GPC 2 -> d = 0
d = 0;
N1 = d+1; %% horizonte de predição inicial
N2 = d+50; %% horizonte predição
Ny = N2-N1+1;

Btil = conv(Gz.num{1}(2:end),[zeros(1,d) 1]); % incorporação do atraso no numerador B


G = zeros(Ny,Ny); % matriz dinâmica G
H = zeros(Ny,nb+d);

[E,F] = diofantina(conv(Gz.den{1},[1 -1]),N1,N2);

for i=N1:N2
    EjB = conv(E(i-N1+1,1:i),Btil);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

    
end
G = G(:,1:Nu);
G,F,H

%% obtenção do ganho irrestrito
Qu = lambda*eye(Nu); % matriz de ponderação dos incrementos de controle
Qe = delta*eye(Ny);  % matriz de ponderaão dos erros futuros

Kmpc = (G'*Qe*G+Qu)\G'*Qe;

Kmpc1 = Kmpc(1,:);

%% vetores de simulação
saidas = zeros(1,nit); % vetor da saída
entradas = zeros(1,nit); % vetor do sinal de controle
du = zeros(1,nit); % vetor dos incrementos de controle


for i=nin+1:nit
    %%% simulação do processo
    saidas(1,i) = -Gzr.den{1}(2:end)*saidas(1,i-1:-1:i-na)' ...
                  +Gzr.num{1}(2:end)*entradas(1,i-dr-1:-1:i-1-nb-dr)';
              
    %%% vetor de referências (caso futuro desconhecido)
    R = refs(i)*ones(Ny,1);
    
    %%% cálculo da resposta livre;
    f = F*saidas(1,i:-1:i-na)';
    
    
    if(~isempty(H))
        f = f+ H*du(1,i-1:-1:i-nb-d)'; % parcela dos incrementos de controle
    end
    
    %%% cálculo do incremento de controle ótimo
    du(1,i) = Kmpc1*(R-f);
    
    %%% cálculo do sinal de controle ótimo
    entradas(1,i) = du(1,i)+entradas(1,i-1);
    
              
end

%%
saidas2 = saidas;
entradas2 = entradas;
du2 = du;
% 
% figure
% subplot(2,1,1)
% plot(saidas2)
% subplot(2,1,2)
% plot(entradas2)


%% Geração das figuras
cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(saidas1(1,nin+1:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(saidas2(1,nin+1:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(refs(1,nin+1:nit),'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 1.5])
hl = legend('GPC1','GPC2','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(entradas1(1,nin+1:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(entradas2(1,nin+1:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
ylim([-1.5 1.5])
ylabel('Manipulada','FontSize', tamletra)
xlabel('Tempo (amostras)','FontSize',tamletra)
grid on
set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7081 0.5377 0.2054 0.1789];
