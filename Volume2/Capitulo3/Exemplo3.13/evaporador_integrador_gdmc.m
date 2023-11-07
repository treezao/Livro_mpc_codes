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
G(1,1) = 3.5/s*exp(-0.1*s);
G(1,2) = -1/(2*s+1)*exp(-0.5*s);
G(2,1) = -2/(1.5*s+1)*exp(-.7*s);
G(2,2) = 1/(3.2*s+1)*exp(-.5*s);


% Definir os modelos que relacionam {saida x,perturbação z};
Gq(1,1) = 3.5/s*exp(-.1*s);
Gq(2,1) = -4.5/(2*s+1)*exp(-0.2*s);

Ts = 0.1; %período de amostragem, OBS: manter atrasos múltiplos inteiros do Ts

%%% discretização do processo
Gz = c2d(G,Ts,'zoh');
Gzq = c2d(Gq,Ts,'zoh'); % perts de entrada


atrasoPlantaInput = totaldelay(Gz);
	
if(mq>0)
    atrasoPlantaDist = totaldelay(Gzq);
end

%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nss = [90 220]; %horizonte de regime permanente - dimensão 1 x n
Nf = [100 220]; % horizonte de modelo do filtro do GDMC - dimensão 1 x n
Nu = [5 5]; %horizontes de controle - dimensão 1 x m
N1 = [1+1 5+1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [1+20 5+20]; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar

delta = [5 1]./N; %ponderação nos erros - dimensão 1 x n
lambda = [5 1]./Nu; %ponderação nas ações de controle - dimensão 1 x m       

af = [0.93;0]; % polos do filtro de referência discreto

%% definição do cenário de simulação
nit =floor(40/Ts); %tempo de simulação
nin =max(max(Nf),max(Nss))+10; %número da iteração inicial
nit = nit+nin;


refs = zeros(n,nit); % vetor de referências
perts = zeros(mq,nit); % vetor de perturbações

refs(1,nin+floor(5/Ts):end) = 1;
refs(2,nin+floor(15/Ts):end) = 0.5;

perts(1,nin+floor(25/Ts):end) = 0.2;

%% Controlador DMC - MIMO

%%% obtenção dos modelos
modDegrauU = {};
%%%%%%%%%
for i = 1:n
    for j = 1:m
        modDegrauU{i,j} = step(Gz(i,j),Ts:Ts:5*Nss(i)*Ts);
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


%% Cálculo do filtro do erro de predição
% Nf = 35;
betaf = 0.8; % formato do filtro (az+b)/(z-betaf)

z = tf('z',Ts);
nf = 1;

pz = 1; % obtem o polo indesejado (instável) de malha aberta


%%% para a saída 1
F{1} = tf(0,1,Ts);
for i=N1(1):N2(1)
    %%% monta sistema para obtenção dos parâmetros do filtro
    %%% primeira equação (z^i - F_i(z) = 0 -> para z=pz
    %%% segunda equação F_i(1) = 0 -> força ganho unitário
    indf = i-N1(1)+1;
    Af = [1 1;
          -betaf, -1];
    bf = [(1-betaf);
          i*(1-betaf)^2];
    X = Af\bf;
    F{1}(indf,1) = (X(1)*z+X(2))/(z-betaf);
    %%% teste da condição
%     zpk(z^i - (X(1)*z+X(2))/(z-betaf))

    %%% armazena coeficientes gtil
    for k=1:m
        modDegrauUF{1,k,i} = filter(F{1}(indf,1).num{1},F{1}(indf,1).den{1},modDegrauU{1,k});    
    end
end

%%% armazena coeficientes gtil
for k=1:m
    %%% calcula a matriz H para o cálculo da resposta livre no caso GDMC
    H1 = [];
    H2 = [];

    for j=N1(1):N2(1)
        H1 = [H1;modDegrauU{1,k}(j+1:j+Nf(1))'];
        H2 = [H2;modDegrauUF{1,k,j}(1:Nf(1))'];

    end
    H{1,k} = H1-H2;

end



%%% calculo das matrizes H para a saída 2 que não utilizará filtros de erro
%%% de predição

%%% calcula a matriz H para o cálculo da resposta livre no caso GDMC


for k=1:m
    H1 = [];
    H2 = [];
    for i=N1(2):N2(2)
        H1 = [H1;modDegrauU{2,k}(i+1:i+Nss(2))'];
        H2 = [H2;modDegrauU{2,k}(1:Nss(2))'];

    end
    H{2,k} = H1-H2;
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
yfilt = zeros(N(1),nit);

fma = {};
for i=1:n
        fma{i} = zeros(Nss(i),1);
end

rf = zeros(n,1);

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
    rf = af.*rf + (1-af).*refs(:,k);
    
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(rf(i),N(i),1)];
    end 
    
    %%% calculo da resposta livre
    f = []; %vetor de resposta livre corrigida
    
    %%%% cálculo da resposta livre filtrada
    %%% saída 1 filtrada
    for i=1:N(1)
        yfilt(i,k) = -F{1}(i,1).den{1}(2:end)*yfilt(i,k-1:-1:k-nf)' + F{1}(i,1).num{1}*saidas(1,k:-1:k-nf)';
    end
    f = [(H{1,1}*du(1,k-1:-1:k-Nf(1))' +H{1,2}*du(2,k-1:-1:k-Nf(1))'  +yfilt(:,k));
         (H{2,1}*du(1,k-1:-1:k-Nss(2))'+H{2,2}*du(2,k-1:-1:k-Nss(2))'+saidas(2,k))];

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

xlabel('Tempo (s)')

hf.Position = tamfigura;
hl.Position = [0.8599 0.5987 0.1154 0.3072];
hl2.Position = [0.8570 0.1872 0.1168 0.1592];

% print('evaporador_integrador_gdmc','-depsc')

