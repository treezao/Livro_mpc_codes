clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');
Ts = 1; % perído de amostragem em minutos

G = 1/s;
z = tf('z',Ts);
Gz = c2d(G,Ts,'zoh'); % modelo discretizado
num = Gz.num{1}; % numerador do modelo discreto
den = Gz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Gz.inputdelay; % atraso discreto


%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 20; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

Nss=80; % horizonte de modelo
Nf = 20; % horizonte de modelo filtrado
betaf = 0.6; % polo do filtro do GDMC

Gcoef = step(Gz,Ts:Ts:Nss*Ts); % coeficientes da resposta ao degrau

%% montando as matrizes do DMC

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

pz = 1; % obtem o polo indesejado (instável) de malha aberta

for i=N1(1):N2(1)
    %%% monta sistema para obtenção dos parâmetros do filtro
    %%% primeira equação (z^i - F_i(z) = 0 -> para z=pz
    %%% segunda equação F_i(1) = 0 -> força ganho unitário
    indf = i-N1(1)+1;
    Af = [-2*betaf,(-1-betaf);
          1 1];
    bf = [i*(1-betaf)^3;
          (1-betaf)^2];
    X = Af\bf;
    F(indf,1) = (X(1)*z^2+X(2)*z)/(z-betaf)^2;
    %%% teste da condição
%     zpk(z^i - (X(1)*z^2+X(2)*z)/(z-betaf)^2)

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

%% inicialização vetores
nit =round(60/Ts); %tempo de simulação
nin = max(Nss,Nf)+10; %número da iteração inicial
nit = nit+nin;


entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+round(30/Ts):end) = 5;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+round(5/Ts):end) = 10;


erro = zeros(nit,1); % vetor de erros
yfilt = zeros(nit,N(1)); % vetor com as saidas filtras


%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = -den(2:end)*saidas(k-1:-1:k-na) + num*(entradas(k-dd:-1:k-nb-dd-1) + perts(k-dd:-1:k-dd-nb-1));
    
    erro(k) = refs(k)-saidas(k);
    
    %% -- Controlador GDMC 
    %%% referencias
    R = refs(k)*ones(N,1);
    
    %%% calculo da resposta livre
    for i=1:N(1)
        yfilt(k,i) = -F(i,1).den{1}(2:end)*yfilt(k-1:-1:k-nf,i) + F(i,1).num{1}*saidas(k:-1:k-nf);
    end
    
    f = H*du(k-1:-1:k-Nf) + yfilt(k,:)';
    
    %% Resolve o problema de otimização
    du(k) = Kdmc1*(R-f);
    entradas(k) = entradas(k-1)+du(k);
    
end

%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;

cores = gray(3);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(t,saidas(vx)+30,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx)+30,'--','LineWidth',tamlinha,'Color',cores(2,:))
hl = legend('Controlada','Referência','Location','NorthEast')
ylabel('Controlada (%)','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t,entradas(vx)+47,'LineWidth',tamlinha,'Color',cores(1,:))

ylabel('Manipulada (%)','FontSize', tamletra)
grid on
xlabel('Tempo (minutos)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7220 0.6314 0.2071 0.1242];





