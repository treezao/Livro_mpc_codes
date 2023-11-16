clear all, 
close all, 
clc
close_system

run('../../../Bibliotecas/parametrosFiguras.m')
%% Modelo nominal
Ts = 1;
z =tf('z',Ts);


d = 1;
Gnz = (0.1*z^-1+0.2*z^-2)/(1-1.1*z^-1)
Lnz = z^-d

P = Gnz*Lnz


%% ajuste DTC-GPC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

m = 1; % número de entradas da planta - manipuladas
n = 1; % número de saidas da planta - controladas


Nu = 5; %horizontes de controle - dimensão 1 x m
N1 = 1; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = 20; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar

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

%%
%%% ajuste do filtro do erro de predição
d = 1
b = 0.9;

z0 = 1.1;
As = [z0 1;
      1 1];
bs= [z0^d*(z0-b);
     (1-b)];
 
X = As\bs;
a1 = X(1);
a2 = X(2);

Fe1 = (a1*z+a2)/(z-b)
Fe1 = zpk(Fe1)

% 1-z0^-d*(a1*z0+a2)/(z0-b)

% Fe = (1-b)/(1-a)*(z-a)/(z-b)
S1 = Gnz*(1-z^-d*Fe1)
S1 = minreal(S1)


%%% ajuste do filtro do erro de predição
d = 5
b = 0.9;

z0 = 1.1;
As = [z0 1;
      1 1];
bs= [z0^d*(z0-b);
     (1-b)];
 
X = As\bs;
a1 = X(1);
a2 = X(2);

Fe5 = (a1*z+a2)/(z-b)
Fe5 = zpk(Fe5)

% 1-z0^-d*(a1*z0+a2)/(z0-b)

% Fe = (1-b)/(1-a)*(z-a)/(z-b)
S5 = Gnz*(1-z^-d*Fe5)
S5 = minreal(S5)

%%% ajuste do filtro do erro de predição
d = 1
b = 0.95;

z0 = 1.1;
As = [z0 1;
      1 1];
bs= [z0^d*(z0-b);
     (1-b)];
 
X = As\bs;
a1 = X(1);
a2 = X(2);

Fe12 = (a1*z+a2)/(z-b)
Fe12 = zpk(Fe12)

% 1-z0^-d*(a1*z0+a2)/(z0-b)

% Fe = (1-b)/(1-a)*(z-a)/(z-b)
S12 = Gnz*(1-z^-d*Fe12)
S12 = minreal(S12)


%%% ajuste do filtro do erro de predição
d = 5
b = 0.95;

z0 = 1.1;
As = [z0 1;
      1 1];
bs= [z0^d*(z0-b);
     (1-b)];
 
X = As\bs;
a1 = X(1);
a2 = X(2);

Fe52 = (a1*z+a2)/(z-b)
Fe52 = zpk(Fe52)

% 1-z0^-d*(a1*z0+a2)/(z0-b)

% Fe = (1-b)/(1-a)*(z-a)/(z-b)
S52 = Gnz*(1-z^-d*Fe52)
S52 = minreal(S52)



w = logspace(-3,log10(2*pi*1/2/Ts),200);
[mag,pha] = bode([Fe1,Fe5,Fe12,Fe52],w);

mag1 = reshape(mag(1,1,:),200,1);
mag2 = reshape(mag(1,2,:),200,1);
mag3 = reshape(mag(1,3,:),200,1);
mag4 = reshape(mag(1,4,:),200,1);


%%
cores = gray(5);
cores = cores(1:end-1,:);

hf = figure
h= subplot(1,1,1)
semilogx(w,20*log10(mag1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w,20*log10(mag2),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx(w,20*log10(mag3),':','LineWidth',tamlinha,'Color',cores(3,:))
semilogx(w,20*log10(mag4),'-.','LineWidth',tamlinha,'Color',cores(4,:))
semilogx(2*pi/2/Ts*[1 1],[-500 500],'LineWidth',tamlinha,'Color',cores(1,:))

grid on

hl = legend('F_{e,1}(z), d=1, b=0,9','F_{e,5}(z), d=5, b=0,9','F_{e,1}(z), d=1, b=0,95','F_{e,5}(z), d=5, b=0,95','Location','NorthWest')
xlabel('Frequência (rad/amostra)','FontSize', tamletra)
ylabel('Magnitude (dB)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


ylim([0 8])
xlim([10^-3 10^1])

hf.Position = tamfigura;
% hl.Position = [0.6810 0.6411 0.2054 0.1806]; 
% print('dtc_gpc_instavel_ex_2_fe1','-depsc')


%% Aumentando a ordem do filtro
%%% ajuste do filtro do erro de predição
d = 1
b = 0.8449; %% msm tempo de assentamento do caso b=0.9 p=1.

z0 = 1.1;
As = [z0^2 z0;
      1 1];
bs= [z0^d*(z0-b)^2;
     (1-b)^2];
 
X = As\bs;
a1 = X(1);
a2 = X(2);

Fe1p = (a1*z^2+a2*z)/(z-b)^2
Fe1p = zpk(Fe1p)

% 1-z0^-d*(a1*z0+a2)/(z0-b)

% Fe = (1-b)/(1-a)*(z-a)/(z-b)
S1p = Gnz*(1-z^-d*Fe1p)
S1p = minreal(S1p)


%%% ajuste do filtro do erro de predição
d = 5
b = 0.8449;

z0 = 1.1;
As = [z0^2 z0;
      1 1];
bs= [z0^d*(z0-b)^2;
     (1-b)^2];
 
X = As\bs;
a1 = X(1);
a2 = X(2);

Fe5p = (a1*z^2+a2*z)/(z-b)^2
Fe5p = zpk(Fe5p)

% 1-z0^-d*(a1*z0+a2)/(z0-b)

% Fe = (1-b)/(1-a)*(z-a)/(z-b)
S5p = Gnz*(1-z^-d*Fe5p)
S5p = minreal(S5p)


w = logspace(-3,log10(2*pi*1/2/Ts),200);
[mag,pha] = bode([Fe1,Fe5,Fe1p,Fe5p],w);

mag1 = reshape(mag(1,1,:),200,1);
mag2 = reshape(mag(1,2,:),200,1);
mag3 = reshape(mag(1,3,:),200,1);
mag4 = reshape(mag(1,4,:),200,1);

%%


hf = figure
h= subplot(1,1,1)
semilogx(w,20*log10(mag1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w,20*log10(mag2),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx(w,20*log10(mag3),':','LineWidth',tamlinha,'Color',cores(3,:))
semilogx(w,20*log10(mag4),'-.','LineWidth',tamlinha,'Color',cores(4,:))
semilogx(2*pi/2/Ts*[1 1],[-500 500],'LineWidth',tamlinha,'Color',cores(1,:))

grid on

hl = legend('F_{e,1}(z), d=1, b=0,9, p=1','F_{e,5}(z), d=5, b=0,9, p=1','F_{e,1}(z), d=1, b=0,8449, p=2','F_{e,5}(z), d=5, b=0,8449, p=2','Location','SouthWest')
xlabel('Frequência (rad/amostra)','FontSize', tamletra)
ylabel('Magnitude (dB)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


ylim([-15 10])
xlim([10^-3 10^1])

hf.Position = tamfigura;
% hl.Position = [0.6810 0.6411 0.2054 0.1806]; 
% print('dtc_gpc_instavel_ex_2_fe2','-depsc')


%%
%%% cenário de simulação
tsim = 300;
tref = 5;
ampref = 0.6;
tpert = 30;
amppert = -0.2;

%%% caso nominal com d=1 e b=0.9
S = S1
Fe = Fe1
P = Gnz*z^-1

out = sim('dtc_gpc_instavel_sim');

deltat = 1;
t1 = out.tout(1:deltat:end);
ref1 = out.simout(1:deltat:end,1);
y1 = out.simout(1:deltat:end,2);
u1 = out.simout(1:deltat:end,3);


%%% caso nominal com d=5 e b=0.9
S = S5
Fe = Fe5
P = Gnz*z^-5

t2 = out.tout(1:deltat:end);
ref2 = out.simout(1:deltat:end,1);
y2 = out.simout(1:deltat:end,2);
u2 = out.simout(1:deltat:end,3);


%%% caso nominal com d=1 e b=0.95
S = S12
Fe = Fe12
P = Gnz*z^-1

out = sim('dtc_gpc_instavel_sim');

deltat = 1;
t3 = out.tout(1:deltat:end);
ref3 = out.simout(1:deltat:end,1);
y3 = out.simout(1:deltat:end,2);
u3 = out.simout(1:deltat:end,3);


%%% caso nominal com d=5 e b=0.95
S = S52
Fe = Fe52
P = Gnz*z^-5

t4 = out.tout(1:deltat:end);
ref4 = out.simout(1:deltat:end,1);
y4 = out.simout(1:deltat:end,2);
u4 = out.simout(1:deltat:end,3);


amppert = 0;
%%% caso erro com d=1 e b=0.9
S = S1
Fe = Fe1
P = Gnz*0.68*z^-1

out = sim('dtc_gpc_instavel_sim');

t5 = out.tout(1:deltat:end);
ref5 = out.simout(1:deltat:end,1);
y5 = out.simout(1:deltat:end,2);
u5 = out.simout(1:deltat:end,3);


%%% caso erro com d=5 e b=0.9
S = S5
Fe = Fe5
P = Gnz*0.68*z^-5

out = sim('dtc_gpc_instavel_sim');

t6 = out.tout(1:deltat:end);
ref6 = out.simout(1:deltat:end,1);
y6 = out.simout(1:deltat:end,2);
u6 = out.simout(1:deltat:end,3);

%%% caso erro com d=1 e b=0.9 p=2
S = S1p
Fe = Fe1p
P = Gnz*0.68*z^-1

out = sim('dtc_gpc_instavel_sim');

t7 = out.tout(1:deltat:end);
ref7 = out.simout(1:deltat:end,1);
y7 = out.simout(1:deltat:end,2);
u7 = out.simout(1:deltat:end,3);


%%% caso erro com d=5 e b=0.9 p=2
S = S5p
Fe = Fe5p
P = Gnz*0.68*z^-5

out = sim('dtc_gpc_instavel_sim');

t8 = out.tout(1:deltat:end);
ref8 = out.simout(1:deltat:end,1);
y8 = out.simout(1:deltat:end,2);
u8 = out.simout(1:deltat:end,3);



%%
cores = gray(4);
cores = cores(1:end-1,:);

%%% figura caso nominal
hf = figure
hf.Position = tamfigura;

h=subplot(2,1,1)
plot(t1,y1,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t3,y3,'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t1,ref1,'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 0.7])
xlim([0 50])
% hl = legend('b=0,9','b=0,95','Referência','Location','SouthEast')
% hl.Position = [0.6889 0.5943 0.2375 0.1619];
ylabel(['Controlada',newline,'Caso d=1'],'FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h1 = h;

h = subplot(2,1,2)
plot(t2,y2,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t4,y4,'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t2,ref2,'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 0.7])
xlim([0 50])
hl = legend('b=0,9','b=0,95','Referência','Location','SouthEast')
% hl.Position = [0.6889 0.5943 0.2375 0.1619];
ylabel(['Controlada',newline,'Caso d=5'],'FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


xlabel('Tempo (amostras)','FontSize',tamletra)

hl.Position = [0.6942 0.5072 0.2054 0.1789]; 
% print('dtc_gpc_instavel_ex_2_nominal','-depsc')


%%
%%% figura caso com erro
hf = figure
hf.Position = tamfigura;

h=subplot(2,1,1)
plot(t5,y5,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t7,y7,'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t5,ref5,'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 1])
xlim([0 150])
hl = legend('p=1','p=2','Referência','Location','SouthEast')

ylabel(['Controlada',newline,'Caso d=1'],'FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on



h = subplot(2,1,2)
plot(t6,y6,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t8,y8,'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t6,ref6,':','LineWidth',tamlinha,'Color',cores(3,:))
ylim([0 1])
xlim([0 150])
% hl = legend('p=1','p=2','Referência','Location','SouthEast')

ylabel(['Controlada',newline,'Caso d=5'],'FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


xlabel('Tempo (amostras)','FontSize',tamletra)

hl.Position = [0.6899 0.5103 0.2054 0.1789]; 
% print('dtc_gpc_instavel_ex_2_erro','-depsc')
