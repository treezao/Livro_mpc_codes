clear all, 
close all, 
clc
close_system

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%% Modelo nominal

Ts = 0.2;
z =tf('z',Ts);
s = tf('s');


Gn = 3/(2*s+1)/(s+1);
L = 5;
Ls = exp(-5*s);

P = Gn*Ls;

d = L/Ts;
Gnz = c2d(Gn,Ts,'zoh')
Lnz = z^-d

%% ajuste GPC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

m = 1; % número de entradas da planta - manipuladas
n = 1; % número de saidas da planta - controladas

atrasoEntrada = L/Ts;

Nu = 2; %horizontes de controle - dimensão 1 x m
N1 = d+1; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = d+20; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar

delta = 1; %ponderação nos erros - dimensão 1 x n
lambda = 2; %ponderação nas ações de controle - dimensão 1 x m       

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

Ftemp = [];
for i=1:n
    Ftemp = blkdiag(Ftemp,F{i});
end
Fgpc = Ftemp;
        
%%% matrizes de ponderação
Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

% H = 2*(G'*Qy*G+Qu);
% invH = inv(H); 

Kgpc = (G'*Qy*G+Qu)\G'*Qy;
Kgpc1 = [];
for i=1:m
    Kgpc1 = [Kgpc1; Kgpc(sum(Nu(1:i-1))+1,:)];
end

%%% preparação das maztrizes por FTs.
zy_gpc = [];
for i=1:n
    ztemp = [];
    for k=1:nF(i)
        ztemp = [ztemp;z^(-(k-1))];
    end
    zy_gpc = blkdiag(zy_gpc,ztemp);
end

zu_gpc = [];
for i=1:m
    ztemp = [];
    for k=1:nH(i)
        ztemp = [ztemp;z^(-k)];
    end
    zu_gpc = blkdiag(zu_gpc,ztemp);
end

zr_gpc = [];
for i=1:n
    zr_gpc = blkdiag(zr_gpc,ones(Ny(i),1));
end

%% ajuste DTC-GPC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

N1 = N1-d; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = N2-d; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar

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
        Btil = Bmfd{i,j}(2:end); 
        
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

Hdtc = cell2mat(H);
G = cell2mat(G);

Ftemp = [];
for i=1:n
    Ftemp = blkdiag(Ftemp,F{i});
end
Fdtc = Ftemp;
        
%%% matrizes de ponderação
Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

% H = 2*(G'*Qy*G+Qu);
% invH = inv(H); 

Kdtc = (G'*Qy*G+Qu)\G'*Qy;
Kdtc1 = [];
for i=1:m
    Kdtc1 = [Kdtc1; Kdtc(sum(Nu(1:i-1))+1,:)];
end

%%% preparação das maztrizes por FTs.
zy_dtc = [];
for i=1:n
    ztemp = []
    for k=1:nF(i)
        ztemp = [ztemp;z^(-(k-1))];
    end
    zy_dtc = blkdiag(zy_dtc,ztemp);
end

zu_dtc = [];
for i=1:m
    ztemp = []
    for k=1:nH(i)
        ztemp = [ztemp;z^(-k)];
    end
    zu_dtc = blkdiag(zu_dtc,ztemp);
end

zr_dtc = [];
for i=1:n
    zr_dtc = blkdiag(zr_dtc,ones(Ny(i),1));
end


%%% ajuste do filtro do erro de predição
a = 0.915;
b = 0.6;
Fe = (1-b)/(1-a)*(z-a)/(z-b)
S = Gnz*(1-z^-d*Fe)





%%
%%% cenário de simulação
tsim = 40;
tref = 2;
ampref = 0.5;
tpert = 20;
amppert = 0.05;
tpert2 = 0;
amppert2 = 0;



out = sim('dtc_gpc_temperatura_sim');

deltat = 5;
t1 = out.tout(1:deltat:end);
ref1 = out.simout(1:deltat:end,1);
y1 = out.simout(1:deltat:end,2:3);
u1 = out.simout(1:deltat:end,4:5);


%% caso 2

%%% reajuste do preditor do DTC-GPC
%%% ajuste do filtro do erro de predição
a = 0.7;
b = 0.9;
Fe = (1-b)/(1-a)*(z-a)/(z-b)
S = Gnz*(1-z^-d*Fe)

out = sim('dtc_gpc_temperatura_sim');


t2 = out.tout(1:deltat:end);
ref2 = out.simout(1:deltat:end,1);
y2 = out.simout(1:deltat:end,2:3);
u2 = out.simout(1:deltat:end,4:5);

%% caso 3 - com, erro de modelagem


Lr = 5.5;
%%% planta real incerta
P = Gn*exp(-Lr*s);

out = sim('dtc_gpc_temperatura_sim');


t3 = out.tout(1:deltat:end);
ref3 = out.simout(1:deltat:end,1);
y3 = out.simout(1:deltat:end,2:3);
u3 = out.simout(1:deltat:end,4:5);

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
% ylim([0 0.7])
xlim([0 40])
hl = legend('GPC','DTC-GPC','Referência','Location','SouthEast')
% hl.Position = [0.6889 0.5943 0.2375 0.1619];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on



h = subplot(2,1,2)
plot(t1,u1(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t1,u1(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([0 40])
% ylim([0 1.3])

set(h, 'FontSize', tamletra);

xlabel('Tempo (minutos)','FontSize',tamletra)



hl.Position = [0.6810 0.6411 0.2054 0.1806]; 
% print('dtc_gpc_temperatura_1','-depsc')






%%
%%% figura caso 2

hf = figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t2,y2(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,y2(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,ref2,'-.','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 0.7])
xlim([0 40])

hl = legend('GPC','DTC-GPC','Referência','Location','SouthEast')
% hl.Position = [0.6889 0.5943 0.2375 0.1619];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on



h = subplot(2,1,2)
plot(t2,u2(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t2,u2(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([0 40])
% ylim([0 1.3])


set(h, 'FontSize', tamletra);

xlabel('Tempo (minutos)','FontSize',tamletra)


hl.Position = [0.6810 0.6411 0.2054 0.1806]; 
% print('dtc_gpc_temperatura_2','-depsc')






%%
%%% figura caso 3

hf = figure
hf.Position = tamfigura;

h=subplot(2,1,1)
plot(t3,y3(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t3,y3(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t3,ref3,'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 0.65])
xlim([0 40])

hl = legend('GPC','DTC-GPC','Referência','Location','SouthEast')
% hl.Position = [0.6889 0.5943 0.2375 0.1619];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on



h = subplot(2,1,2)
plot(t3,u3(:,1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t3,u3(:,2),'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([0 40])
ylim([0 0.4])


set(h, 'FontSize', tamletra);

xlabel('Tempo (minutos)','FontSize',tamletra)


hl.Position = [0.6935 0.4863 0.2054 0.1806]; 
% print('dtc_gpc_temperatura_3','-depsc')