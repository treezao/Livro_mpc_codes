clear all, 
close all, 
clc


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

Gz = c2d(G,Ts,'zoh')
Lz = z^-(L1/Ts);

Pz = Gz*Lz;

Gzq = Pz;


%% ajuste GPC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

m = 1; % número de entradas da planta - manipuladas
n = 1; % número de saidas da planta - controladas
mq = 1; % número de perturabções

Nu = 5; %horizontes de controle - dimensão 1 x m
N1 = 11; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = 25; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar
Nss = 80;

delta = 1; %ponderação nos erros - dimensão 1 x n
lambda = 1; %ponderação nas ações de controle - dimensão 1 x m


%%% filtro do erro de predição
%%% (z^d-Fe) = 0 para z=0.9048
%%% Fe(z) = (az-b)/(z-c)
%%% ->(a-b)/(1-c) = 1
%%% ->z^d(z-c) -(az-b) = 0, para z = 0.9048
c = 0
Ab = [1, -1;
      0.9048,-1];
bb = [(1-c);
      0.9048^(L1/Ts)*(0.9048-c)];
X = Ab\bb

Fe{1} = tf([X(1) -X(2)],[1 -c],Ts)



%% definição do cenário de simulação
tsim = 25;
tref = 4;
tpert = 15;


nit =round(tsim/Ts); %tempo de simulação
nin = 20; % número da iteração inicial
nit = nit+nin;

refs = zeros(n,nit+max(N2));
perts = zeros(mq,nit+max(N2));

refs(1,nin+floor(tref/Ts):end) = 1;

perts(1,nin+floor(tpert/Ts):end) = .5;

%% montagem das matrizes de predição

%%% obtenção dos modelos
modDegrauU = {};
%%%%%%%%%
for i = 1:n
    for j = 1:m
        modDegrauU{i,j} = step(Pz(i,j),Ts:Ts:Nss(i)*Ts);
    end
end


%% Monta as matrizes offline que definem o problema MIMO DMC
Gm = {};
for k=1:n 
    for j=1:m
        for i=1:Nu(j)
            Gm{k,j}(i:Ny(k),i) = modDegrauU{k,j}(N1(k):N2(k)-i+1);
        end
    end
    
end
G = cell2mat(Gm);

Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

Kdmc = (G'*Qy*G+Qu)\G'*Qy;

Kdmc1 = [];
for i=1:m
    Kdmc1 = [Kdmc1; Kdmc(sum(Nu(1:i-1))+1,:)];
end


%% inicializacao dos estados e variaveis da simulaçãoo
% - obtencao dos numerados e denominadores
numR = cell(n,m);
denR = cell(n,m);
numFE = cell(n,1);
denFE = cell(n,1);

numDR = cell(n,mq);
denDR = cell(n,mq);


for i=1:n
	% obtencao das respostas ao degrau das entradas
	for j=1:m
		% obtencao dos numeradores e denominadores para o simulador
        [numR{i}{j}, denR{i}{j}] = tfdata(Pz(i,j),'v');
	end
	% obtencao das respostas ao degrau das perts
    for j=1:mq
        [numDR{i}{j},denDR{i}{j}]=tfdata(Gzq(i,j),'v');
    end    
    
    [numFE{i},denFE{i}] = tfdata(Fe{i},'v');
    
end


estadosEntr = cell(n,m);
for i=1:n
    for j=1:m
		estadosEntr{i}{j} = zeros(1,nit);
    end
end

estadosPert = cell(n,mq);
for i=1:n
	for j=1:mq
		estadosPert{i}{j}= zeros(1,nit);
	end    
end

saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle


fma = {};
for i=1:n
        fma{i} = zeros(Nss(i),1);
end

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numR{i}{j},2);
            sizeD = size(denR{i}{j},2);

            str1 = numR{i}{j}*entradas(j,k:-1:k+1-sizeN)';
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end

        for j=1:mq
            sizeN = size(numDR{i}{j}, 2);
            sizeD = size(denDR{i}{j}, 2);

            str1 = numDR{i}{j}*perts(j,k:-1:k+1-sizeN)';
            str2 = -denDR{i}{j}(1,2:end)*estadosPert{i}{j}(k-1:-1:k-sizeD+1)';
            estadosPert{i}{j}(k) = str1+str2;

            aux = aux + estadosPert{i}{j}(k);
        end
        saidas(i,k) = aux ;%;+ 0.01*randn(1,1);
    end
    %% -- Controlador DMC 
    %%% referencias
    R = refs(:,k)*ones(Ny,1);
    
    %%% calculo da resposta livre
    ep = zeros(n,1); % erros de prediçãoo
    
    f = []; %vetor de resposta livre corrigida
    
    %%%% atualizaçãoo da resposta livre
    for i = 1:n
        for j = 1:m
            fma{i} = fma{i} + modDegrauU{i,j}*du(j,k-1);
        end
        
        ep(i) = saidas(i,k)-fma{i}(1);
        
        fma{i} = [fma{i}(2:Nss(i)); fma{i}(Nss(i))];
        f = [f; fma{i}(N1(i):N2(i)) + ones(Ny(i),1)*ep(i)];
    end
        
    %% Resolve o problema de otimização
    
    du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

%%
t = (0:nit-nin)*Ts;
ref = refs(nin:nit);
y1 = saidas(nin:nit);
u1 = entradas(nin:nit);

%% Caso 2 - DMC filtrado
saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

epf = zeros(n,nit); % erro de predição filtrado
ep = zeros(n,nit); % erro de predição

fma = {};
for i=1:n
        fma{i} = zeros(Nss(i),1);
end

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numR{i}{j},2);
            sizeD = size(denR{i}{j},2);

            str1 = numR{i}{j}*entradas(j,k:-1:k+1-sizeN)';
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end

        for j=1:mq
            sizeN = size(numDR{i}{j}, 2);
            sizeD = size(denDR{i}{j}, 2);

            str1 = numDR{i}{j}*perts(j,k:-1:k+1-sizeN)';
            str2 = -denDR{i}{j}(1,2:end)*estadosPert{i}{j}(k-1:-1:k-sizeD+1)';
            estadosPert{i}{j}(k) = str1+str2;

            aux = aux + estadosPert{i}{j}(k);
        end
        saidas(i,k) = aux ;%;+ 0.01*randn(1,1);
    end
    %% -- Controlador DMC 
    %%% referencias
    R = refs(:,k)*ones(Ny,1);
    
    %%% calculo da resposta livre
    f = []; %vetor de resposta livre corrigida
    
    %%%% atualizaçãoo da resposta livre
    for i = 1:n
        for j = 1:m
            fma{i} = fma{i} + modDegrauU{i,j}*du(j,k-1);
        end
        
        ep(i,k) = saidas(i,k)-fma{i}(1);
        
        sizeFEN = size(numFE{i},2);
        sizeFED = size(denFE{i},2);
        epf(i,k) = -denFE{i}(1,2:end)*epf(i,k-1:-1:k-sizeFED+1)' + numFE{i}*ep(i,k:-1:k-sizeFEN+1)';
        
        fma{i} = [fma{i}(2:Nss(i)); fma{i}(Nss(i))];
        f = [f; fma{i}(N1(i):N2(i)) + ones(Ny(i),1)*epf(i,k)];
    end
        
    %% Resolve o problema de otimização
    
    du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

%%
y2 = saidas(nin:nit);
u2 = entradas(nin:nit);


%%
cores = gray(4);
cores = cores(1:end-1,:);


%%% figura caso 1
hf = figure
hf.Position = tamfigura;

h=subplot(2,1,1)
plot(t,y1,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,y2,'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,ref,'-.','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 0.7])
xlim([0 25])
hl = legend('DMC','DMC filtrado','Referência','Location','SouthEast')

ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

% h.YTickLabel = trocaponto(h.YTickLabel)

h = subplot(2,1,2)
plot(t,u1,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,u2,'--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
xlim([0 25])
% ylim([0 1.3])

set(h, 'FontSize', tamletra);

xlabel('Tempo (segundos)','FontSize',tamletra)
% h.YTickLabel = trocaponto(h.YTickLabel)


hl.Position = [0.6685 0.5444 0.2232 0.1806];
% print('trocador_PS_dmc_dtc','-depsc')
% print('trocador_PS_dmc_dtc_pb','-deps')

