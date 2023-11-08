% clear all, 
% close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%% modelo não-linear e seus parâmetros
global k1 k2 Vc Qf
k1 = 10;
k2 = 10;
Vc = 1;
Qf = 1/30;


Cfbar = 3.288
Cobar = 1.316

%% modelo linearizado
s = tf('s');

Ts = 5; % Periodo de amostragem do processo em segundos
z = tf('z',Ts);

G = 3.433/(103.1*s-1); % Modelo por função de transferência do processo
Gz = c2d(G,Ts,'zoh'); % modelo discretizado

%% parâmetros de ajuste
N1 = 1; %horizonte de predição inicial
N2 = 30; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

umax = [0.05]+Cfbar; % valor máximo do sinal de controle
umin = [-0.05]+Cfbar; % valor mínimo do sinal de controle
dumax = [0.025]; % valor máximo do incremento de controle
dumin = -dumax; % valor mínimo do incremento de controle
ymax = 0.12+Cobar; % valor máximo para a saída
ymin = -0.05+Cobar; % valor mínimo para a saída


Nss=40; % horizonte de modelo
Nf = 25; % horizonte de modelo filtrado
betaf = 0.8; % polo do filtro do gdmc

Gcoef = step(G,Ts:Ts:2*Nss*Ts); % coeficientes da resposta ao degrau


%% montando as matrizes do DMC recursivo

G = zeros(N2,Nu);
G(:,1) = Gcoef(1:N2,1);


for i=2:Nu
    G(i:end,i) = G(1:end-(i-1),1);    
end

G = G(N1:end,:);

Qy = delta*eye(N2-N1+1);
Qu = lambda*eye(Nu);

Kdmc = inv(G'*Qy*G+Qu)*G'*Qy

Kdmc1 = Kdmc(1,:);

%%% cálculo dos filtros dos erros de predição (SISO) para o GDMC

F = tf(0,1,Ts);
nf = 2; %ordem do filtro

pz = pole(Gz); % obtem o polo indesejado (instável) de malha aberta

for i=N1(1):N2(1)
    %%% monta sistema para obtenção dos parâmetros do filtro
    %%% primeira equação (z^i - F_i(z) = 0 -> para z=pz
    %%% segunda equação F_i(1) = 0 -> força ganho unitário
    indf = i-N1(1)+1;
    Af = [pz^2 pz;
          1 1];
    bf = [pz^i*(pz-betaf)^2;
          (1-betaf)^2];
    X = Af\bf;
    F(indf,1) = (X(1)*z^2+X(2)*z)/(z-betaf)^2;
    %%% teste da condição
%     pz^i-(X(1)*pz+X(2))/(pz-betaf)^2

    %%% armazena coeficientes gtil
    modDegrauUF{i} = filter(F(indf,1).num{1},F(indf,1).den{1},Gcoef);

end


%%% calcula a matriz H para o cálculo da resposta livre no caso GDMC
H1 = [];
H2 = [];

for i=N1(1):N2(1)
    H1 = [H1;Gcoef(i+1:i+Nf)'];
    H2 = [H2;modDegrauUF{i}(1:Nf)'];
    
end
H = H1-H2



%% matrizes para o caso com restrições
Hqp = 2*(G'*Qy*G+Qu);
fqp1 = -2*G'*Qy; 

LB = repelem(dumin,Nu')';
UB = repelem(dumax,Nu')';
Rbar = tril(ones(Nu));
Rbar = [Rbar;-Rbar;G;-G];


%% inicialização vetores
nin = Nss+1;
nit = round(400/Ts) + nin; % número de iterações da simulação

entradas = Cfbar*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = Cobar*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,2); % vetor com as perturbações do sistema
perts(nin+round(200/Ts):end) = -.03;

refs = Cobar*ones(nit,1); % vetor de referências
refs(nin+round(50/Ts):end) = 0.1+Cobar;

erro = zeros(nit,1); % vetor de erros

yfilt = Cobar*ones(nit,N); % saídas filtradas

%% simulação com referencia futura
for i = nin:nit
    %% simulação do modelo processo
    saidas(i) = simModelo(saidas(i-1),entradas(i-1)+perts(i-1),Ts);
    
    erro(i) = refs(i)-saidas(i);
    %% Controlador
    %%% referencias
    R = refs(i)*ones(N,1);
    
    %%% calculo da resposta livre
    for k=1:N(1)
        yfilt(i,k) = -F(k,1).den{1}(2:end)*yfilt(i-1:-1:i-nf,k) + F(k,1).num{1}*saidas(i:-1:i-nf);
    end
    
    f = H*du(i-1:-1:i-Nf) + yfilt(i,:)';
    
    %%% calculo do incremento de controle ótimo    
    % com restrições
    fqp = fqp1*(R-f);
    
    rbar = [repelem(umax-entradas(i-1),Nu)';
             repelem(entradas(i-1)-umin,Nu)';
             ymax-f;
             -ymin+f];
    [X,FVAL,EXITFLAG] = quadprog(Hqp,fqp,Rbar,rbar,[],[],LB,UB);
    
    %%% caso dê infactível, remover a restrição na saída
    if(EXITFLAG==-2)
        rbar = [repelem(umax-entradas(i-1),Nu)';
                 repelem(entradas(i-1)-umin,Nu)'];
        X = quadprog(Hqp,fqp,Rbar(1:2*Nu,:),rbar,[],[],LB,UB);
    end

    du(i) = X(1);
    
    
    % sem restrições
%     du(i) = Kdmc1*(R-f);
    
    %%% aplicar no sistema
    entradas(i) = entradas(i-1)+du(i);    
    
end

%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;

cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(3,1,1)
plot(t,saidas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx),'--','LineWidth',tamlinha,'Color',cores(2,:))
ylim([1.3 1.45])
hl = legend('GDMC','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('C_o (mol/L)','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
h.YTick = [3.2400 3.2800 3.3200]
ylabel('C_f (mol/L)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

h = subplot(3,1,3)
plot(t,du(vx),'LineWidth',tamlinha,'Color',cores(1,:))
ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (s)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.7370 0.6532 0.2054 0.1231];

