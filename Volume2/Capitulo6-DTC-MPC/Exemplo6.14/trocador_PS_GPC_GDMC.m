clear all, 
close all, 
clc
close_system

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%
Ts = 0.15;
z =tf('z',Ts);
s = tf('s');

%% caso 1
Ke = 1;
tau = 1.5;
G = Ke/(tau*s+1)
L1 = 10*Ts;

P = G*exp(-L1*s);

%%% ajuste PS - controlador primário
Czps = 3.70*(z-0.77)/(z-1);
Fzps = 0.23/(z-0.77);

Gnz = c2d(G,Ts,'zoh')
d = (L1/Ts);
Lnz = z^-(d);



%% ajuste GPC e GDMC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

m = 1; % número de entradas da planta - manipuladas
n = 1; % número de saidas da planta - controladas

atrasoEntrada = L1/Ts;

Nu = 5; %horizontes de controle - dimensão 1 x m
N1 = 11; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = 25; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar
Nf = 40;% para gdmc somente;

delta = 1; %ponderação nos erros - dimensão 1 x n
lambda = 1; %ponderação nas ações de controle - dimensão 1 x m      

betaf = 0.8; % para gdmc somente, polo do filtro

%% Configuração GPC
[Amfd,Bmfd] = MFDredux([Gnz.den],[Gnz.num]);
%%% montagem das matrizes de predição

F = {};
H = {};
G = {};

nH = zeros(1,m); % qt de dados passados necessários das entradas
nF = zeros(1,n); % qt de dados passados necessários das saídas.

for i=1:n
    
    [Ei,Fi] = diofantina(conv(Amfd{i},[1 -1]),N1(i),N2(i));
    
    F{i} = Fi;
    nF(i) = size(Fi,2);
    
    for j=1:m
        % incorporação do atraso no numerador B
        Btil = conv(Bmfd{i,j}(2:end),[zeros(1,atrasoEntrada(i,j)) 1]); 
        
        Gtemp = [];
        Htemp = [];
        for k=N1(i):N2(i)
            EjB = conv(Ei(k-N1(i)+1,1:k),Btil);
            
            Gtemp(k-N1(i)+1,1:k) = EjB(k:-1:1);    
            Htemp(k-N1(i)+1,:) = EjB(k+1:end);
        end
        G{i,j} = Gtemp(:,1:Nu(j));
        H{i,j} = Htemp;
        
        nhi = size(Htemp,2);
        if(nhi>nH(j))
            nH(j)=nhi;
        end
    end
end

%%% ajuste dos tamanhos das matrizes Hi e Hqi para concatenação

for i=1:n
    for j=1:m
        nHi = size(H{i,j},2);
        if(nHi<nH(j))
            H{i,j} = [H{i,j} zeros(Ny(i),nH(j)-nHi)];
        end
    end    
end

Hgpc = cell2mat(H);
G = cell2mat(G);
Ggpc = G;

Ftemp = [];
for i=1:n
    Ftemp = blkdiag(Ftemp,F{i});
end
Fgpc = Ftemp;
        
%%% matrizes de ponderação
Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

Kgpc = (G'*Qy*G+Qu)\G'*Qy;
Kgpc1 = [];
for i=1:m
    Kgpc1 = [Kgpc1; Kgpc(sum(Nu(1:i-1))+1,:)];
end

%%% preparação das maztrizes por FTs.
zy = [];
for i=1:n
    ztemp = []
    for k=1:nF(i)
        ztemp = [ztemp;z^(-(k-1))];
    end
    zy = blkdiag(zy,ztemp);
end

zu = [];
for i=1:m
    ztemp = []
    for k=1:nH(i)
        ztemp = [ztemp;z^(-k)];
    end
    zu = blkdiag(zu,ztemp);
end

zr = [];
for i=1:n
    zr = blkdiag(zr,ones(Ny(i),1));
end

%% configuração GDMC
modDegrauU = step(Gnz*Lnz,Ts:Ts:3*Nf*Ts);

%%% Monta as matrizes offline que definem o problema DMC
G = [];
for i=1:Nu
    G(i:Ny,i) = modDegrauU(N1:N2-i+1);
end

Kgdmc = (G'*Qy*G+Qu)\G'*Qy; % utilizado apenas no caso irrestrito
Kgdmc1 = Kgdmc(1,:);


%%% cálculo dos filtros dos erros de predição (SISO) para o GDMC
Fegdmc = tf(0,1,Ts);
nf = 1;

pz1 = pole(Gnz); % polo da FT da entrada

for i=N1(1):N2(1)
    %%% monta sistema para obtenção dos parâmetros do filtro
    %%% primeira equação (z^i - F_i(z) = 0 -> para z=pz1
    %%% segunda equação (z^i - F_i(z) = 0 -> para z=pz2
    %%% terceira equação F_i(1) = 0 -> força ganho unitário
    indf = i-N1(1)+1;
    Af = [pz1 1;
          1 1];
    bf = [pz1^i*(pz1-betaf);
          (1-betaf)];
    X = Af\bf;
    Fegdmc(indf,1) = (X(1)*z+X(2))/(z-betaf);
    %%% teste da condição
%     pz1^i-(X(1)*pz1^2+X(2)*pz1)/(pz1-betaf)^2

    %%% armazena coeficientes gtil
    modDegrauUF{i} = filter(Fegdmc(indf,1).num{1},Fegdmc(indf,1).den{1},modDegrauU);

end


%%% calcula a matriz H para o cálculo da resposta livre no caso GDMC
H1 = [];
H2 = [];



for i=N1(1):N2(1)
    H1 = [H1;modDegrauU(i+1:i+Nf)'];
    H2 = [H2;modDegrauUF{i}(1:Nf)'];
        
end
Hgdmc = H1-H2;

zudmc = [];
for k=1:size(Hgdmc,2)
    zudmc = [zudmc;z^(-k)];
end

%%
%%% cenário de simulação
tsim = 35;
tref = 4;
ampref = 1;
tpert = 15;
amppert = 0.5;
tpert2 = 25;
amppert2 = -0.1;



out = sim('trocador_PS_GPC_GDMC_sim');

deltat = 5;
t1 = out.tout(1:deltat:end);
ref1 = out.simout(1:deltat:end,1);
y1 = out.simout(1:deltat:end,2:4);
u1 = out.simout(1:deltat:end,5:7);


%% caso 2
Lr = 12*Ts;
%%% planta real incerta
P = Ke/(1*s+1)*exp(-Lr*s);

out = sim('trocador_PS_GPC_GDMC_sim');


t2 = out.tout(1:deltat:end);
ref2 = out.simout(1:deltat:end,1);
y2 = out.simout(1:deltat:end,2:4);
u2 = out.simout(1:deltat:end,5:7);

%%
cores = gray(5);
cores = cores(1:end-1,:);

%%% figura caso 1
hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t1,y1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,y1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,y1(:,3),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t1,ref1,'-.','LineWidth',tamlinha,'Color',cores(4,:))
% ylim([0 0.7])
xlim([0 35])
hl = legend('GPC','PS','GDMC','Referência','Location','SouthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h = subplot(2,1,2)
plot(t1,u1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,u1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,u1(:,3),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([0 35])
% ylim([0 1.3])

set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)


hl.Position = [0.6827 0.5226 0.2054 0.2371]; 
% print('trocador_PS_GPC_GDMC_1','-depsc')

%%
%%% figura caso 1- perturbacao
hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t1,y1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,y1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,y1(:,3),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t1,ref1,'-.','LineWidth',tamlinha,'Color',cores(4,:))
% ylim([0.38 0.6])
xlim([15 35])
hl = legend('GPC','PS','GDMC','Referência','Location','SouthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h1 = h;

h = subplot(2,1,2)
plot(t1,u1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,u1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,u1(:,3),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([15 35])
% ylim([0.3 1.2])

set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)


hl.Position = [0.6827 0.5226 0.2054 0.2371]; 
% print('trocador_PS_GPC_GDMC_1_pert','-depsc')

%%
%%% figura caso 2

hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t2,y2(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,y2(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,y2(:,3),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t2,ref2,'-.','LineWidth',tamlinha,'Color',cores(4,:))
ylim([0 1.5])
xlim([0 35])

hl = legend('GPC','PS','GDMC','Referência','Location','SouthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h = subplot(2,1,2)
plot(t2,u2(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,u2(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,u2(:,3),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([0 35])
ylim([0 2])


set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)


hl.Position = [0.6827 0.5226 0.2054 0.2371]; 
% print('trocador_PS_GPC_GDMC_2','-depsc')