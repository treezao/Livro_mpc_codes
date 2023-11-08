clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%

% MIMO GPC - fracionador

m = 3; % número de entradas da planta - manipuladas
n = 3; % número de saidas da planta - controladas
mq = 0; % número de perturbações da planta

s = tf('s');
% Definir os modelos que relacionam {saida x,entrada y};
G(1,1) = 4.05/(50*s+1)*exp(-27*s);
G(1,2) = 1.77/(60*s+1)*exp(-28*s);
G(1,3) = 5.88/(50*s+1)*exp(-27*s);
G(2,1) = 5.39/(50*s+1)*exp(-18*s);
G(2,2) = 5.72/(60*s+1)*exp(-14*s);
G(2,3) = 6.9/(40*s+1)*exp(-15*s);
G(3,1) = 4.38/(33*s+1)*exp(-20*s);
G(3,2) = 4.42/(44*s+1)*exp(-22*s);
G(3,3) = 7.2/(19*s+1);

Gq = []*tf(1,1);

% Definir os modelos que relacionam {saida x,perturbação z};

Ts = 1.0; %período de amostragem, OBS: manter atrasos múltiplos inteiros do Ts

%%% discretização do processo
Gz = c2d(G,Ts,'zoh');
Gzq = c2d(Gq,Ts,'zoh'); % perts de entrada


atrasoEntrada = totaldelay(Gz);
	
if(mq>0)
    atrasoPert = totaldelay(Gzq);
end

%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [5 5 5]; %horizontes de controle - dimensão 1 x m
N1 = [28 15 1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [88 75 45]; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar
Nq = [0]; % horizonte de predição da perturbação - dimensão 1 x mq

delta = [1 1 1]./N; %ponderação nos erros - dimensão 1 x n
lambda = [0.1 0.1 0.1]./Nu; %ponderação nas ações de controle - dimensão 1 x m       



%% definição do cenário de simulação
nit =floor(300/Ts); %tempo de simulação
nin = 50; %número da iteração inicial
nit = nit+nin;


refs = zeros(n,nit);
perts = zeros(mq,nit);

refs(1,nin+floor(10/Ts):end) = 1;
refs(2,nin+floor(100/Ts):end) = 0.5;
refs(3,nin+floor(170/Ts):end) = 0.2;

% perts(1,nin+floor(130/Ts):end) = 0.2;

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
            H{i,j} = [H{i,j} zeros(N(i),nH(j)-nHi)];
        end
    end
    
    for j=1:mq
        nHi = size(Hq{i,j},2);
        if(nHi<nHq(j))
            Hq{i,j} = [Hq{i,j} zeros(N(i),nHq(j)-nHi)];
        end
    end        
end

H = cell2mat(H);
Hq = cell2mat(Hq);
G = cell2mat(G);
Gq = cell2mat(Gq);

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

Kdmc = (G'*Qe*G+Qu)\G'*Qe;
Kdmc1 = [];
for i=1:m
    Kdmc1 = [Kdmc1; Kdmc(sum(Nu(1:i-1))+1,:)];
end

%%% inicialização dos vetores para o cálculo da resposta livre sem correção


%% inicializacao dos estados e variaveis da simulação
% - obtencao dos numerados e denominadores
numR = cell(n,m);
denR = cell(n,m);
numDR = cell(n,mq);
denDR = cell(n,mq);


for i=1:n
	% obtencao das respostas ao degrau das entradas
	for j=1:m
		% obtencao dos numeradores e denominadores para o simulador
        [numR{i}{j}, denR{i}{j}] = tfdata(Gz(i,j),'v');
	end
    
	% obtencao das respostas ao degrau das perts
    for j=1:mq
        [numDR{i}{j},denDR{i}{j}]=tfdata(Gzq(i,j),'v');
	end
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
dq = zeros(mq,nit); % vetor de incrementos de perturbações

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numR{i}{j},2);
            sizeD = size(denR{i}{j},2);

            str1 = numR{i}{j}*entradas(j,k-atrasoEntrada(i,j):-1:k+1-sizeN-atrasoEntrada(i,j))';
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end

        for j=1:mq
            sizeN = size(numDR{i}{j}, 2);
            sizeD = size(denDR{i}{j}, 2);

            str1 = numDR{i}{j}*perts(j,k-atrasoPert(i,j):-1:k+1-sizeN-atrasoPert(i,j))';
            str2 = -denDR{i}{j}(1,2:end)*estadosPert{i}{j}(k-1:-1:k-sizeD+1)';
            estadosPert{i}{j}(k) = str1+str2;

            aux = aux + estadosPert{i}{j}(k);
        end
        saidas(i,k) = aux ;%;+ 0.01*randn(1,1);
    end
    %% -- Controlador GPC 
    dq(:,k) = perts(:,k)-perts(:,k-1);
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),N(i),1)];
    end 
    
    %%% calculo da resposta livre
    du1 = [];
    for i=1:m
        du1 = [du1;du(i,k-1:-1:k-nH(i))'];
    end
    
    y1 = [];
    for i=1:n
        y1 = [y1;saidas(i,k:-1:k-nF(i)+1)'];
    end
    
    dq1 = [];
    for i=1:mq
        dq1 = [dq1;dq(i,k:-1:k-nHq(i)+1)'];
    end

    f = F*y1 + H*du1;
%     f = F*y1 + H*du1 + Hq*dq1;


    %% Resolve o problema de otimização
    du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end




%% Gera gráficos
cores = gray(5);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin)*Ts;
ind = nin:nit;


hf= figure
h=subplot(2,1,1)
plot(t,saidas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidas(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidas(3,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))
plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))
plot(t,refs(3,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))


hl = legend('y_1','y_2','y_3','r_1,r_2,r_3','Location','SouthEast')
ylim([-.3 1.2])
ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h=subplot(2,1,2)
plot(t,entradas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradas(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradas(3,ind),':','LineWidth',tamlinha,'Color',cores(3,:))


hl2 = legend('u_1','u_2','u_3','Location','SouthEast')
ylim([-1.3 1.4])
ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo (min)')

hf.Position = tamfigura;
hl.Position = [0.8024 0.6379 0.1643 0.3145];
hl2.Position = [0.8435 0.2250 0.1161 0.2387];

% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))

% print('gpc_fracionador_3','-depsc')

