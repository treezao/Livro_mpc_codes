clear all
close all
clc

addpath('../../../Bibliotecas/')
%%
% MIMO DMC - nivel e temperatura 

m = 2; % número de entradas da planta - manipuladas
n = 2; % número de saidas da planta - controladas
mq = 1; % número de perturbações da planta

s = tf('s');
% Definir os modelos que relacionam {saida x,entrada y};
G(1,1) = 1.43/(7.5*s+1)*exp(-s);
G(1,2) = -4.5/(13.4*s+1)*exp(-s);
G(2,1) = 0;
G(2,2) = 5.72/(14.7*s+1);

% Definir os modelos que relacionam {saida x,perturbação z};
Gq(1,1) = 2.1/(13*s+1)*exp(-s);
Gq(2,1) = -3.12/(14.7*s+1);

Ts = 1.0; %período de amostragem, OBS: manter atrasos múltiplos inteiros do Ts

%%% discretização do processo
Gz = c2d(G,Ts,'zoh');
Gzq = c2d(Gq,Ts,'zoh'); % perts de entrada


atrasoEntrada = totaldelay(Gz);
	
if(mq>0)
    atrasoPert = totaldelay(Gzq);
end

%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [2 1]; %horizontes de controle - dimensão 1 x m
N1 = [1 1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [3 3]; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar
Nq = [2]; % horizonte de predição da perturbação - dimensão 1 x mq

delta = [1 9]./Ny; %ponderação nos erros - dimensão 1 x n
lambda = [1 20]./Nu; %ponderação nas ações de controle - dimensão 1 x m       

%% Controlador GPC - MIMO
%%% obtenção da representação MFD
[Amfd,Bmfd] = MFDredux([Gz.den Gzq.den],[Gz.num,Gzq.num]);

Bqmfd = Bmfd(:,m+1:end);
Bmfd = Bmfd(:,1:m);

%%% montagem das matrizes de predição

F = {};
H = {};
Hq = {};
G = {};
Gq ={};

nH = zeros(1,m); % qt de dados passados necessários das entradas
nHq = zeros(1,mq); %qt de dados passados necessários das perturbações
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
        
    for j=1:mq
        % incorporação do atraso no numerador Bq
        Bqtil = conv(Bqmfd{i,j},[zeros(1,atrasoPert(i,j)) 1]); 
    
        Gqtemp = [];
        Hqtemp = [];
        for k=N1(i):N2(i)
            EjBq = conv(Ei(k-N1(i)+1,1:k),Bqtil);

            Gqtemp(k-N1(i)+1,1:k) = EjBq(k:-1:1);    
            Hqtemp(k-N1(i)+1,:) = EjBq(k+1:end);
        end
        
        Gq{i,j} = Gqtemp(:,1:Nq(j));
        Hq{i,j} = Hqtemp;

        nhi = size(Hqtemp,2);
        if(nhi>nHq(j))
            nHq(j)=nhi;
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
    
    for j=1:mq
        nHi = size(Hq{i,j},2);
        if(nHi<nHq(j))
            Hq{i,j} = [Hq{i,j} zeros(Ny(i),nHq(j)-nHi)];
        end
    end        
end

H = cell2mat(H)
Hq = cell2mat(Hq)
G = cell2mat(G)
Gq = cell2mat(Gq)

Ftemp = [];
for i=1:n
    Ftemp = blkdiag(Ftemp,F{i});
end
F = Ftemp
   

        
%%% matrizes de ponderação
Qe = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle


Kmpc = (G'*Qe*G+Qu)\G'*Qe;
Kmpc1 = [];
for i=1:m
    Kmpc1 = [Kmpc1; Kmpc(sum(Nu(1:i-1))+1,:)];
end

