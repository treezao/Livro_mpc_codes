clear all
close all
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
% MIMO DMC - nivel e temperatura 

m = 2; % número de entradas da planta - manipuladas
n = 2; %número de saidas da planta - controladas
mq = 1; % número de entradas da planta - perturbações

s = tf('s');
% Definir os modelos que relacionam {saida x,entrada y};
G(1,1) = 1.43/(7.5*s+1)*exp(-s);
G(1,2) = -4.5/(13.4*s+1)*exp(-s);
G(2,1) = 0*s;
G(2,2) = 5.72/(14.7*s+1);

% Definir os modelos que relacionam {saida x,perturbação z};
Gq(1,1) = 2.1/(13*s+1)*exp(-s);
Gq(2,1) = -3.12/(14.7*s+1);

Ts = 1.0; %período de amostragem, OBS: manter atrasos múltiplos inteiros do Ts

%%% discretização do processo
Gz = c2d(G,Ts,'zoh');
Gzq = c2d(Gq,Ts,'zoh'); % perts de entrada


atrasoPlantaInput = totaldelay(Gz);
	
if(mq>0)
    atrasoPlantaDist = totaldelay(Gzq);
end

%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nss = [50 50]; %horizonte de regime permanente - dimensão 1 x n
Nu = [5 10]; %horizontes de controle - dimensão 1 x m
N1 = [2 1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [17 31]; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar

delta = [1 9]./N; %ponderação nos erros - dimensão 1 x n
lambda = [1 20]./Nu; %ponderação nas ações de controle - dimensão 1 x m       



%% definição do cenário de simulação
nit =200; %tempo de simulação
nin = 10; %número da iteração inicial
nit = nit+nin;


refs = zeros(n,nit); % vetor de referências
perts = zeros(mq,nit); % vetor de perturbações

refs(1,nin+10:end) = 1;
refs(2,nin+80:end) = 0.5;

perts(1,nin+130:end) = 0.2;

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

Qy = diag(repelem(delta,1,N)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

 

%% Monta as matrizes offline que definem o problema MIMO DMC

for k=1:n 
    for j=1:m
        for i=1:Nu(j)
            Gm{k,j}(i:N(k),i) = modDegrauU{k,j}(N1(k):N2(k)-i+1);
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

            str1 = numR{i}{j}*entradas(j,k-atrasoPlantaInput(i,j):-1:k+1-sizeN-atrasoPlantaInput(i,j))';
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end

        for j=1:mq
            sizeN = size(numDR{i}{j}, 2);
            sizeD = size(denDR{i}{j}, 2);

            str1 = numDR{i}{j}*perts(j,k-atrasoPlantaDist(i,j):-1:k+1-sizeN-atrasoPlantaDist(i,j))';
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
    eta = zeros(n,1); % erros de predição
    
    f = []; %vetor de resposta livre corrigida
    
    %%%% atualização da resposta livre
    for i = 1:n
        for j = 1:m
            fma{i} = fma{i} + modDegrauU{i,j}*du(j,k-1);
        end
        
        eta(i) = saidas(i,k)-fma{i}(1);
        
        fma{i} = [fma{i}(2:Nss(i)); fma{i}(Nss(i))];
        f = [f; fma{i}(N1(i):N2(i)) + ones(N(i),1)*eta(i)];
    end

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

plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))

hl = legend('y_1','y_2','r_1','r_2','Location','SouthEast')

ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h=subplot(2,1,2)
plot(t,entradas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradas(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))

hl2 = legend('u_1','u_2','Location','SouthEast')

ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo [s]')

hf.Position = tamfigura;
hl.Position = [0.8570 0.5548 0.1154 0.3072];
hl2.Position = [0.8556 0.2155 0.1168 0.1592];

% print('tanque1_controle','-depsc')

