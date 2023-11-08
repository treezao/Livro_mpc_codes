clear all
close all
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
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

Nu = [5 10]; %horizontes de controle - dimensão 1 x m
N1 = [2 1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [17 31]; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar
Nq = [0]; % horizonte de predição da perturbação - dimensão 1 x mq

delta = [1 9]./N; %ponderação nos erros - dimensão 1 x n
lambda = [1 30]./Nu; %ponderação nas ações de controle - dimensão 1 x m       



%% definição do cenário de simulação
nit =200; %tempo de simulação
nin = 10; %número da iteração inicial
nit = nit+nin;


refs = zeros(n,nit);
perts = zeros(mq,nit+max(Nq));

refs(1,nin+10:end) = 1;
refs(2,nin+80:end) = 0.5;

perts(1,nin+130:end) = 0.2;

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

Kmpc = (G'*Qe*G+Qu)\G'*Qe;
Kmpc1 = [];
for i=1:m
    Kmpc1 = [Kmpc1; Kmpc(sum(Nu(1:i-1))+1,:)];
end

%%% inicialização dos vetores para o cálculo da resposta livre sem correção


%% inicializacao dos estados e variaveis da simulação - GPC sem FF
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
    
    f = F*y1 + H*du1;


    %% Resolve o problema de otimização
    du(:,k) = Kmpc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasSFF = saidas;
entradasSFF = entradas;
duSFF = du;









%% simulação GPC com FF

%%% "zerando" estados do simulador
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
    dqfut = [];
    for i=1:mq
        dq1 = [dq1;dq(i,k:-1:k-nHq(i)+1)'];
        dqfut = [dqfut;(perts(i,k+1:k+Nq(i))-perts(i,k:k+Nq(i)-1))'];
    end

    f = F*y1 + H*du1 + Hq*dq1;
    
    if(~isempty(Gq))
        f = f+Gq*dqfut;
    end

    %% Resolve o problema de otimização
    du(:,k) = Kmpc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasCFF = saidas;
entradasCFF = entradas;
duCFF = du;


%% Gera gráficos
t = ((nin:nit)-nin)*Ts;
ind = nin:nit;

cores = gray(4);
cores = cores(1:end-1,:);


hf= figure
h=subplot(2,1,1)
plot(t,saidasSFF(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasCFF(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,refs(1,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(3,:))

plot(t,saidasSFF(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,saidasCFF(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,refs(2,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(3,:))


hl = legend('S/ FF','C/ FF','Referências','Location','SouthEast')


ta1 = annotation('textarrow');
ta1.String = 'y_1';

ta2 = annotation('textarrow');
ta2.String = 'y_2';

xlim([0 120])

ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h=subplot(2,1,2)
plot(t,entradasSFF(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasCFF(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

plot(t,entradasSFF(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,entradasCFF(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

xlim([0 120])

ta3 = annotation('textarrow');
ta3.String = 'u_1';

ta4 = annotation('textarrow');
ta4.String = 'u_2';


ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo (segundos)')


hf.Position = tamfigura;
hl.Position = [0.7685 0.4665 0.2168 0.1789];
ta1.Position = [0.3075 0.8187 -0.0571 0.0439];
ta2.Position = [0.7889 0.8497 -0.0629 -0.0465];
ta3.Position = [0.3418 0.3929 -0.0800 -0.0155];
ta4.Position = [0.4886 0.2086 -0.0668 -0.0273];

% print('tanque1_gpc_controle','-depsc')


%%

hf= figure
h=subplot(2,1,1)
plot(t,saidasSFF(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasCFF(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,refs(1,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(3,:))

plot(t,saidasSFF(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,saidasCFF(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))
plot(t,refs(2,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(3,:))


hl = legend('S/ FF','C/ FF','Referências','Location','SouthEast')


ta1 = annotation('textarrow');
ta1.String = 'y_1';

ta2 = annotation('textarrow');
ta2.String = 'y_2';

xlim([120 160])

ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h=subplot(2,1,2)
plot(t,entradasSFF(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasCFF(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

plot(t,entradasSFF(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,entradasCFF(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

xlim([120 160])

ta3 = annotation('textarrow');
ta3.String = 'u_1';

ta4 = annotation('textarrow');
ta4.String = 'u_2';


ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo (segundos)')

hf.Position = tamfigura;
hl.Position = [0.7170 0.4665 0.2168 0.1789];
ta1.Position = [0.3661 0.8084 -0.0629 0.0697];
ta2.Position = [0.2618 0.7232 0.0471 -0.0516];
ta3.Position = [0.4989 0.3542 -0.0529 0.0568];
ta4.Position = [0.6432 0.2587 -0.0957 -0.0335];


% print('tanque1_gpc_controle_q','-depsc')
