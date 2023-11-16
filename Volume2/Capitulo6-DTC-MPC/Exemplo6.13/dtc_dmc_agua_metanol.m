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
z = tf('z',Ts);

atrasoEntrada = totaldelay(Gz);
	
if(mq>0)
    atrasoPert = totaldelay(Gzq);
end

%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [5 5]; %horizontes de controle - dimensão 1 x m
N1 = [2 3]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [40 30]; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar
Nq = 0; % % horizonte de predição da perturbação - dimensão 1 x mq

Nss = [100 100]; %horizonte de regime permanente - dimensão 1 x n - apenas par ao DMC

delta = [1 1]./Ny; %ponderação nos erros - dimensão 1 x n
lambda = [1 1]./Nu; %ponderação nas ações de controle - dimensão 1 x m       

Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle


%% ajuste do filtro do erro de predição
pi1 = exp(-Ts/14.9) % polo a ser cancelado para a primeira saída
pi2 = exp(-Ts/13.2) % polo a ser cancelado para a segunda saída

A = [1 1;
     pi1 1];
B = [1;pi1^1*pi1]
X1 = A\B
Fe{1} = (X1(1)*z+X1(2))/z

1-pi1^-1*(X1(1)*pi1+X1(2))/pi1

A = [1 1;
     pi2 1];
B = [1;pi2^3*pi2]
X2 = A\B
Fe{2} = (X2(1)*z+X2(2))/z

1-pi2^-3*(X2(1)*pi2+X2(2))/pi2

%% definição do cenário de simulação
nit =80; %tempo de simulação
nin = 20; %número da iteração inicial
nit = nit+nin;


refs = zeros(n,nit);
perts = zeros(mq,nit);

perts(1,nin+floor(10/Ts):end) = -0.2;


%% Controlador DMC - MIMO

%%% obtenção dos modelos
modDegrauU = {};
modDegrauQ = {}; 
%%%%%%%%%
for i = 1:n
    for j = 1:m
        modDegrauU{i,j} = step(Gz(i,j),Ts:Ts:Nss(i)*Ts);
    end
    for j = 1:mq
        modDegrauQ{i,j} = step(Gzq(i,j),0:Ts:(Nss(i)-1)*Ts);
    end
end


%% Monta as matrizes offline que definem o problema MIMO DMC

for k=1:n 
    for j=1:m
        for i=1:Nu(j)
            Gm{k,j}(i:Ny(k),i) = modDegrauU{k,j}(N1(k):N2(k)-i+1);
        end
    end
end
G = cell2mat(Gm);

% H = 2*(G'*Qy*G+Qu);
% invH = inv(H); 

Kdmc = (G'*Qy*G+Qu)\G'*Qy;

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
             repelem(refs(i,k),Ny(i),1)];
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
        f = [f; fma{i}(N1(i):N2(i)) + ones(Ny(i),1)*ep(i)];
    end

    %% Resolve o problema de otimização
    du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end


figure
subplot(2,1,1)
plot(saidas')
subplot(2,1,2)
plot(entradas')

saidasDMC = saidas;
entradasDMC = entradas;
duDMC = du;

%% Controlador DMC com Fe
%%% inicializacao dos estados e variaveis da simulação
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

ep = zeros(n,nit); % vetor de erros de predição
epf = zeros(n,nit); % vetor de erros de predição filtrados

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
             repelem(refs(i,k),Ny(i),1)];
    end 
    
    %%% calculo da resposta livre
    f = []; %vetor de resposta livre corrigida
    
    %%%% atualização da resposta livre
    for i = 1:n
        for j = 1:m
            fma{i} = fma{i} + modDegrauU{i,j}*du(j,k-1);
        end
        
        ep(i,k) = saidas(i,k)-fma{i}(1);
        epf(i,k) = Fe{i}.num{1}*ep(i,k:-1:k-1)';
        
        fma{i} = [fma{i}(2:Nss(i)); fma{i}(Nss(i))];
        f = [f; fma{i}(N1(i):N2(i)) + ones(Ny(i),1)*epf(i,k)];
    end

    %% Resolve o problema de otimização
    du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end


figure
subplot(2,1,1)
plot(saidas')
subplot(2,1,2)
plot(entradas')

saidasDMC2 = saidas;
entradasDMC2 = entradas;
duDMC2 = du;

%% Gera gráficos
cores = gray(4);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin)*Ts;
ind = nin:nit;


hf= figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t,saidasDMC(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasDMC2(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

plot(t,refs(1,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(3,:))


plot(t,saidasDMC(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,saidasDMC2(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))

plot(t,refs(2,ind),'-.','Color',[0 0 0],'LineWidth',tamlinha,'Color',cores(3,:))


hl = legend('DMC','DMC com F_e(z)','Referências','Location','SouthEast')

ta1 = annotation('textarrow');
ta1.String = 'y_1';

ta2 = annotation('textarrow');
ta2.String = 'y_2';

% xlim([110 150])
% ylim([-0.1 1.3])

ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h=subplot(2,1,2)
plot(t,entradasDMC(1,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasDMC2(1,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))


plot(t,entradasDMC(2,ind),'Color',[0 0.4470 0.7410],'LineWidth',tamlinha,'Color',cores(1,:))
plot(t,entradasDMC2(2,ind),'--','Color',[1 0 0],'LineWidth',tamlinha,'Color',cores(2,:))


% xlim([110 150])
% ylim([-0.3 0.3])
ta3 = annotation('textarrow');
ta3.String = 'u_1';

ta4 = annotation('textarrow');
ta4.String = 'u_2';


ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo (segundos)')

% hl.Position = [0.7365 0.4780 0.2179 0.2371];



ta1.Position = [0.3839 0.9065 -0.0624 -0.0292];
ta2.Position = [0.2214 0.7581 0.0518 0.0258];
ta3.Position = [0.3679 0.4290 -0.0681 -0.0401];
ta4.Position = [0.5000 0.1968 -0.0829 0.0076];


% print('dtc_dmc_agua_metanol','-depsc')

