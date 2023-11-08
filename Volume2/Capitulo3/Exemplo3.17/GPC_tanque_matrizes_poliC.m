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

C{1} = [1 -0.8];%% polinômios C
C{2} = [1 -.9]; 


%% definição do cenário de simulação
nit =200; %tempo de simulação
nin = 10; %número da iteração inicial
nit = nit+nin;


refs = zeros(n,nit);
perts = zeros(mq,nit);

refs(1,nin+10:end) = 1;
refs(2,nin+80:end) = 0.5;

perts(1,nin+130:end) = 0.2;

%% Controlador GPC - MIMO

%%% obtenção da representação MFD

% [Amfd,Bmfd] = MFD([Gz.den Gzq.den],[Gz.num,Gzq.num]);
[Amfd,Bmfd] = MFDredux([Gz.den Gzq.den],[Gz.num,Gzq.num]);

%%%

Bqmfd = Bmfd(:,m+1:end);
Bmfd = Bmfd(:,1:m);

%%% montagem das matrizes de predição

F = {};
H = {};
G = {};

nH = zeros(1,m); % qt de dados passados necessários das entradas
nF = zeros(1,n); % qt de dados passados necessários das saídas.




for i=1:n
    
    [Ei,Fi] = diofantinaC(conv(Amfd{i},[1 -1]),C{i},N1(i),N2(i));
    
    F{i} = Fi;
    nF(i) = size(Fi,2);
    
    for j=1:m
        % incorporação do atraso no numerador B
        Btil = conv(Bmfd{i,j}(2:end),[zeros(1,atrasoEntrada(i,j)) 1]); 
        
        Gtemp = [];
        Htemp = [];
        for k=N1(i):N2(i)
            EjB = conv(Ei(k-N1(i)+1,1:k),Btil);
            
            [M,N] = diofantinaC(C{i},EjB,k,k);
            
            Gtemp(k-N1(i)+1,1:k) = M(end:-1:1);    
            Htemp(k-N1(i)+1,:) = N;
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

Hu = [];
for j=1:m
    Htemp = [];
    for i=1:n
        Htemp = blkdiag(Htemp,H{i,j});
    end
    Hu = [Hu, Htemp];
end
        
G = cell2mat(G);

Ftemp = [];
for i=1:n
    Ftemp = blkdiag(Ftemp,F{i});
end
F = Ftemp;
   

%%%
G
F
Htemp = H;
H = [blkdiag(Htemp{1,1},Htemp{2,1}),blkdiag(Htemp{1,2},Htemp{2,2})]
H_simplificado = [[Htemp{1,1}; zeros(size(Htemp{2,1},1),size(Htemp{1,1},2))],blkdiag(Htemp{1,2},Htemp{2,2})]

        
%% matrizes de ponderação
Qe = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle


Kmpc = (G'*Qe*G+Qu)\G'*Qe;
Kmpc1 = [];
for i=1:m
    Kmpc1 = [Kmpc1; Kmpc(sum(Nu(1:i-1))+1,:)];
end


