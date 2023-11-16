clear all, 
close all, 
clc
close_system

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
Lnz = z^-(L1/Ts);


%% ajuste GPC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

m = 1; % número de entradas da planta - manipuladas
n = 1; % número de saidas da planta - controladas

atrasoEntrada = L1/Ts;

Nu = 5; %horizontes de controle - dimensão 1 x m
N1 = 11; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = 25; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar

delta = 1; %ponderação nos erros - dimensão 1 x n
lambda = 1; %ponderação nas ações de controle - dimensão 1 x m       

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
            H{i,j} = [H{i,j} zeros(N(i),nH(j)-nHi)];
        end
    end    
end

H = cell2mat(H);
G = cell2mat(G);

Ftemp = [];
for i=1:n
    Ftemp = blkdiag(Ftemp,F{i});
end
F = Ftemp;
        
%%% matrizes de ponderação
Qe = diag(repelem(delta,1,N)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

% H = 2*(G'*Qe*G+Qu);
% invH = inv(H); 

Kgpc = (G'*Qe*G+Qu)\G'*Qe;
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
    zr = blkdiag(zr,ones(N(i),1));
end

%%
%%% cenário de simulação
tsim = 35;
tref = 4;
ampref = 0.5;
tpert = 15;
amppert = 0.3;
tpert2 = 25;
amppert2 = -0.1;



out = sim('trocador_PS_GPC_sim');

deltat = 5;
t1 = out.tout(1:deltat:end);
ref1 = out.simout(1:deltat:end,1);
y1 = out.simout(1:deltat:end,2:3);
u1 = out.simout(1:deltat:end,4:5);


%% caso 2
Lr = 12*Ts;
%%% planta real incerta
P = Ke/(1*s+1)*exp(-Lr*s);

out = sim('trocador_PS_GPC_sim');


t2 = out.tout(1:deltat:end);
ref2 = out.simout(1:deltat:end,1);
y2 = out.simout(1:deltat:end,2:3);
u2 = out.simout(1:deltat:end,4:5);

%%
cores = gray(4);
cores = cores(1:end-1,:);

%%% figura caso 1
hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t1,y1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,y1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,ref1,'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 0.7])
xlim([0 35])
hl = legend('GPC','PS','Referência','Location','SouthEast')
% hl.Position = [0.6889 0.5943 0.2375 0.1619];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h1 = h;

h = subplot(2,1,2)
plot(t1,u1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,u1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([0 35])
ylim([0 1.3])

set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)

% h.YTickLabel = trocaponto(h.YTickLabel)
% h1.YTickLabel{2} = '0,5';


hl.Position = [0.6845 0.5863 0.2054 0.1806]; 
% print('trocador_PS_GPC_1','-depsc')

%%
%%% figura caso 2

hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t2,y2(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,y2(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,ref2,'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 0.7])
xlim([0 35])

hl = legend('GPC','PS','Referência','Location','SouthEast')
% hl.Position = [0.6889 0.5943 0.2375 0.1619];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

% h.YTickLabel = trocaponto(h.YTickLabel)
h1 = h;

h = subplot(2,1,2)
plot(t2,u2(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,u2(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([0 35])
ylim([0 1.3])


set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)
% h.YTickLabel = trocaponto(h.YTickLabel)
% h1.YTickLabel{2} = '0,5';

hl.Position = [0.6845 0.5863 0.2054 0.1806]; 
% print('trocador_PS_GPC_2','-depsc')