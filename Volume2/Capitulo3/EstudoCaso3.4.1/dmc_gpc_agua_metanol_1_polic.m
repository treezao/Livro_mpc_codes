clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%
% MIMO DMC - nivel e temperatura 

m = 2; % número de entradas da planta - manipuladas
n = 2; %número de saidas da planta - controladas
mq = 1;

s = tf('s');
% Definir os modelos que relacionam {saida x,entrada y};
G(1,1) = 12.8/(16.7*s+1)*exp(-s);
G(1,2) = -18.9/(21*s+1)*exp(-3*s);
G(2,1) = 6.6/(10.9*s+1)*exp(-7*s);
G(2,2) = -19.4/(14.4*s+1)*exp(-3*s);

% Definir os modelos que relacionam {saida x,perturbação z};
Gq(1,1) = 3.8/(14.9*s+1)*exp(-8*s);
Gq(2,1) = 4.9/(13.2*s+1)*exp(-4*s);


Ts = 1.0; %período de amostragem, OBS: manter atrasos múltiplos inteiros do Ts

%%% discretização do processo
Gz = c2d(G,Ts,'zoh');
Gzq = c2d(Gq,Ts,'zoh'); % perts de entrada


atrasoEntrada = totaldelay(Gz);
	
if(mq>0)
    atrasoPert = totaldelay(Gzq);
end

%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [5 5]; %horizontes de controle - dimensão 1 x m
N1 = [2 3]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [30 40]; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar
Nq = 0;

Nss = [100 100]; %horizonte de regime permanente - dimensão 1 x n - apenas par ao DMC

delta = [1 1]./N; %ponderação nos erros - dimensão 1 x n
lambda = [1 1]./Nu; %ponderação nas ações de controle - dimensão 1 x m       

Qe = diag(repelem(delta,1,N)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

C{1} = [1 -0.9];% polinomios C para o caso GPC-C
C{2} = [1 -0.9]; 


%% definição do cenário de simulação
nit =160; %tempo de simulação
nin = 20; %número da iteração inicial
nit = nit+nin;


refs = zeros(n,nit);
perts = zeros(mq,nit);

refs(1,nin+floor(10/Ts):end) = 1;
refs(2,nin+floor(60/Ts):end) = 0.5;

perts(1,nin+floor(110/Ts):end) = -0.2;






%% Controlador DMC - MIMO

%%% obtenção dos modelos
modDegrauU = {};
modDegrauQ = {}; 
%%%%%%%%%
for i = 1:n
    for j = 1:m
        modDegrauU{i,j} = step(Gz(i,j),Ts:Ts:Nss(i)*Ts);
    end
end


%% Monta as matrizes offline que definem o problema MIMO DMC

for k=1:n 
    for j=1:m
        for i=1:Nu(j)
            Gm{k,j}(i:N(k),i) = modDegrauU{k,j}(N1(k):N2(k)-i+1);
        end
    end
end
G = cell2mat(Gm);

% H = 2*(G'*Qe*G+Qu);
% invH = inv(H); 

Kdmc = (G'*Qe*G+Qu)\G'*Qe;

Kdmc1 = [];
for i=1:m
    Kdmc1 = [Kdmc1; Kdmc(sum(Nu(1:i-1))+1,:)];
end


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
    %% -- Controlador DMC 
    
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),N(i),1)];
    end 
    
    %%% calculo da resposta livre
    ep = zeros(n,1); % erros de predição
    
    f = []; %vetor de resposta livre corrigida
    
    %%%% atualização da resposta livre
    for i = 1:n
        for j = 1:m
            fma{i} = fma{i} + modDegrauU{i,j}*du(j,k-1);
        end
        
        ep(i) = saidas(i,k)-fma{i}(1);
        
        fma{i} = [fma{i}(2:Nss(i)); fma{i}(Nss(i))];
        f = [f; fma{i}(N1(i):N2(i)) + ones(N(i),1)*ep(i)];
    end

    %% Resolve o problema de otimização
    du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasDMC = saidas;
entradasDMC = entradas;
duDMC = du;








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
   

        
Kgpc = (G'*Qe*G+Qu)\G'*Qe;
Kgpc1 = [];
for i=1:m
    Kgpc1 = [Kgpc1; Kgpc(sum(Nu(1:i-1))+1,:)];
end


%% inicializacao dos estados e variaveis da simulação
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
    
    f = F*y1 + H*du1;

    %% Resolve o problema de otimização
    du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end


saidasGPC = saidas;
entradasGPC = entradas;
duGPC = du;










%% Controlador GPC - MIMO - polinômio C
%%% obtenção da representação MFD

[Amfd,Bmfd] = MFDredux([Gz.den Gzq.den],[Gz.num,Gzq.num]);

Bqmfd = Bmfd(:,m+1:end);
Bmfd = Bmfd(:,1:m);

%%% montagem das matrizes de predição
F = {};
H = {};
G = {};

nH = zeros(1,m); % qt de dados passados necessários das entradas
nF = zeros(1,n); % qt de dados passados necessários das saídas.
nC = zeros(1,n); % qt de dados passados necessários do poli C

for i=1:n
    
    nC(i) = size(C{i},2);
    
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

            [Mi,Ni] = diofantinaC(C{i},EjB,k,k);
    
            Gtemp(k-N1(i)+1,1:k) = Mi(k:-1:1);    
            Htemp(k-N1(i)+1,:) = Ni;
        end
        G{i,j} = Gtemp(:,1:Nu(j));
        H{i,j} = Htemp;
        
        nH(i,j) = size(H{i,j},2);
    end
end

G = cell2mat(G);

Kgpc = (G'*Qe*G+Qu)\G'*Qe;
Kgpc1 = [];
for i=1:m
    Kgpc1 = [Kgpc1; Kgpc(sum(Nu(1:i-1))+1,:)];
end


%% inicializacao dos estados e variaveis da simulação
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

saidasC = zeros(n,nit);
for i=1:n
    for j=1:m
        duC{i,j} = zeros(1,nit);
    end
end

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
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),N(i),1)];
    end 
    
    %%% calculo da resposta livre
    yc = [];
    f = [];
    for i=1:n
        saidasC(i,k) = saidas(i,k) -C{i}(2:end)*saidasC(i,k-1:-1:k-nC(i)+1)'; 
        
        ftemp = F{i}*saidasC(i,k:-1:k-nF(i)+1)';
        
        for j=1:m
            duC{i,j}(1,k-1) = -C{i}(2:end)*duC{i,j}(k-2:-1:k-nC(i))' + du(j,k-1);
            
            ftemp = ftemp + H{i,j}*duC{i,j}(1,k-1:-1:k-nH(i,j))';
        end
        f = [f;ftemp];
    end

    
    %% Resolve o problema de otimização
    du(:,k) = Kgpc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end


saidasGPCc = saidas;
entradasGPCc = entradas;
duGPCc = du;

%% Gera gráficos
cores = gray(5);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin)*Ts;
ind = nin:nit;


hf= figure
h=subplot(2,1,1)
plot(t,saidasDMC(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGPC(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPCc(1,ind),':','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(1,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(4,:))

plot(t,saidasDMC(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,saidasGPC(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPCc(2,ind),':','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(2,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(4,:))


hl = legend('DMC','GPC','GPC-C','Referências','Location','SouthEast')

% hl.Position = [0.6916 0.4864 0.2518 0.1619];

ta1 = annotation('textarrow');
% ta1.Position = [0.3179 0.8071 -0.0732 0.0452];
ta1.String = 'y_1';

ta2 = annotation('textarrow');
% ta2.Position = [0.7535 0.7094 -0.0732 0.0452];
ta2.String = 'y_2';

xlim([0 110])
ylim([-.15 1.1])
ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h=subplot(2,1,2)
plot(t,entradasDMC(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGPC(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPCc(1,ind),':','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(3,:))


plot(t,entradasDMC(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,entradasGPC(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPCc(2,ind),':','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(3,:))

xlim([0 110])
ylim([-0.1 0.4])

ta3 = annotation('textarrow');
% ta3.Position = [0.3036 0.4072 -0.0715 -0.0215];
ta3.String = 'u_1';

ta4 = annotation('textarrow');
% ta4.Position = [0.4571 0.1714 -0.0589 0.0405];
ta4.String = 'u_2';


ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo (segundos)')

hf.Position = tamfigura;
hl.Position = [0.7365 0.4780 0.2179 0.2371];
ta1.Position = [0.3179 0.8071 -0.0732 0.0452];
ta2.Position = [0.6704 0.7026 -0.0629 0.0258];
ta3.Position = [0.3036 0.4072 -0.0715 -0.0215];
ta4.Position = [0.4032 0.1994 -0.0600 0.0258];

% print('dmc_gpc_agua_metanol_1_r','-depsc')

%%
hf= figure
h=subplot(2,1,1)
plot(t,saidasDMC(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGPC(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPCc(1,ind),':','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(1,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(4,:))


plot(t,saidasDMC(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,saidasGPC(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPCc(2,ind),':','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(2,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(4,:))


hl = legend('DMC','GPC','GPC-C','Referências','Location','SouthEast')

ta1 = annotation('textarrow');
ta1.String = 'y_1';

ta2 = annotation('textarrow');
ta2.String = 'y_2';

xlim([110 150])
ylim([-0.1 1.3])

ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h=subplot(2,1,2)
plot(t,entradasDMC(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGPC(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPCc(1,ind),':','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(3,:))


plot(t,entradasDMC(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,entradasGPC(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPCc(2,ind),':','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(3,:))


xlim([110 150])
ylim([-0.3 0.3])
ta3 = annotation('textarrow');
ta3.String = 'u_1';

ta4 = annotation('textarrow');
ta4.String = 'u_2';


ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo (segundos)')

hf.Position = tamfigura;
hl.Position = [0.7365 0.4780 0.2179 0.2371];
ta1.Position = [0.6143 0.9032 -0.0714 -0.0194];
ta2.Position = [0.6500 0.7935 -0.0633 -0.0271];
ta3.Position = [0.4561 0.4161 -0.0671 -0.0465];
ta4.Position = [0.5840 0.2258 -0.0686 0.0335];

% print('dmc_gpc_agua_metanol_1_q','-depsc')

