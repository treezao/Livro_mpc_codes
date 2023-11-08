clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%% parâmetros do sistema 4 tanques
global A1 A2 A3 A4 a1 a2 a3 a4 g gamma1 gamma2 k1 k2 Ts
    
A1 = 28;
A2 = 32;
A3 = 28;
A4 = 32;
a1 = 0.071;
a2 = 0.057;
a3 = 0.071;
a4 = 0.057;
g = 981;
gamma1 = 0.7;
gamma2 = 0.6;
k1 = 3.33;
k2 = 3.35;

%%% ponto de operação
% h0 = [12.4;12.7;1.8;1.4];
h0 = [12.2630;12.7831;1.6339;1.4090];
v0 = [3;3];

%%% modelo em espaço de estados linearizado em torno do ponto de operação
T1 = A1/a1*sqrt(2*h0(1)/g);
T2 = A2/a2*sqrt(2*h0(2)/g);
T3 = A3/a3*sqrt(2*h0(3)/g);
T4 = A4/a4*sqrt(2*h0(4)/g);

Aee = [-1/T1, 0, A3/A1/T3, 0;
       0, -1/T2, 0, A4/A2/T4;
       0, 0, -1/T3, 0;
       0,0,0,-1/T4];
Bee = [gamma1*k1/A1,0;
       0,gamma2*k2/A2;
       0,(1-gamma2)*k2/A3;
       (1-gamma1)*k1/A4 0];
Cee = eye(4);

Gee = ss(Aee,Bee,Cee,[])

G = tf(Gee); % obtenção da matriz de transferência do sistema

Ts = 2.0; %período de amostragem, OBS: manter atrasos múltiplos inteiros do Ts

%%

m = 2; % número de entradas da planta - manipuladas
n = 4; % número de saidas da planta - controladas
mq = 0; % número de perturbações da planta

%%% discretização do processo
Gz = c2d(G,Ts,'zoh');

atrasoEntrada = totaldelay(Gz);


%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [5 5]; %horizontes de controle - dimensão 1 x m
N1 = [1 1 1 1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [50 50 50 50]; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar

delta = [1 1 1 1]./N; %ponderação nos erros - dimensão 1 x n
lambda = [1 1]./Nu; %ponderação nas ações de controle - dimensão 1 x m       

umax = [12 12];
umin = [0 0];

%% definição do cenário de simulação
nit =floor(190/Ts); %tempo de simulação
nin = 20; %número da iteração inicial
nit = nit+nin;


refs = repmat(h0,[1,nit]);
perts = zeros(mq,nit);

refs(1,nin+floor(5/Ts):end) = h0(1)*1.2;
refs(2,nin+floor(60/Ts):end) = h0(2)*0.8;

refs(1,nin+floor(120/Ts):end) = [8.7664];
refs(2,nin+floor(120/Ts):end) = [9.3579];
refs(3,nin+floor(120/Ts):end) = [1.2273];
refs(4,nin+floor(120/Ts):end) = [0.9785];

% perts(1,nin+130:end) = 0.2;

%% Controlador GPC - MIMO

%%% obtenção da representação MFD
[Amfd,Bmfd] = MFDredux(Gz.den,Gz.num);

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
        
end

%%% ajuste dos tamanhos das matrizes Hi e Hqi para concatenação

for i=1:n
    for j=1:m
        nHi = size(H{i,j},2);
        if(nHi<nH(j))
            H{i,j} = [H{i,j} zeros(N(i),nH(j)-nHi)];
        end
    end
    

end

H = cell2mat(H);
G = cell2mat(G);

Ftemp = [];
for i=1:n
    Ftemp = blkdiag(Ftemp,F{i});
end
F = Ftemp;
   

        
%%% matrizes de ponderação
Qe = diag(repelem(delta,1,N)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

Kgpc = (G'*Qe*G+Qu)\G'*Qe;
Kgpc1 = [];
for i=1:m
    Kgpc1 = [Kgpc1; Kgpc(sum(Nu(1:i-1))+1,:)];
end

%%% matrizes do problema de otimização
Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;

Rbar = [];
for i=1:m
    Rbar = blkdiag(Rbar,tril(ones(Nu(i))));
end
Rbar = [Rbar;-Rbar];

%% inicializacao dos estados e variaveis da simulação
saidas = repmat(h0,[1,nit]); % vetor da saída
entradas = repmat(v0,[1,nit]); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    saidas(:,k) = modelo4tanques(saidas(:,k-1),entradas(:,k-1));
    
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
    
    %%% montagem das matrizes de restrição
    rbar = repelem(umax'-entradas(:,k-1),Nu');
    rbar = [rbar;
            repelem(-umin'+entradas(:,k-1),Nu')];
    
    fqp = fqp1*(R-f);
    X = quadprog(Hqp,fqp,Rbar,rbar);
    for i=1:m
        du(i,k) = X(sum(Nu(1:i-1))+1,1);
    end
    %% Resolve o problema de otimizacao sem restricoes
%     du(:,k) = Kgpc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end



%% Gera gráficos
cores = gray(4);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin)*Ts;
ind = nin:nit;


hf= figure
h=subplot(3,1,1)
plot(t,saidas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidas(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))

plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))

hl = legend('y_1','y_2','r_1,r_2','Location','SouthEast')
hl.Position = [0.8070 0.7694 0.1382 0.1719];

ylabel(['Controladas',newline,'(cm)'],'FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on
xlim([0 190])
ylim([8 15])

h=subplot(3,1,2)
plot(t,saidas(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidas(4,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))

plot(t,refs(3,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,refs(4,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))

h2 = legend('y_3','y_4','r_3','r_4','Location','SouthEast')
h2.Position = [0.8742 0.4309 0.1154 0.2264];

ylabel(['Controladas',newline,'(cm)'],'FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on
xlim([0 190])
ylim([0 2.6])

h=subplot(3,1,3)
plot(t,entradas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradas(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))

h3 = legend('u_1','u_2','Location','SouthEast')
h3.Position = [0.8642 0.2379 0.1168 0.1174];

ylabel(['Manipuladas',newline,'(V)'],'FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlim([0 190])
ylim([-.5 5])
xlabel('Tempo (segundos)')

hf.Position = [488 342 560 420];
% print('gpc_4tanques_1','-depsc')


%% modelo do sistema 4 tanques
function h = modelo4tanques(h1,u1)
    
    global A1 A2 A3 A4 a1 a2 a3 a4 g gamma1 gamma2 k1 k2 Ts
    
    opts = odeset('NonNegative',[1 2 3 4]);
    
    f4tanques= @(t,y) [-a1/A1*sqrt(2*g*y(1))+a3/A1*sqrt(2*g*y(3))+gamma1*k1/A1*u1(1);
                       -a2/A2*sqrt(2*g*y(2))+a4/A2*sqrt(2*g*y(4))+gamma2*k2/A2*u1(2);
                       -a3/A3*sqrt(2*g*y(3))+(1-gamma2)*k2/A3*u1(2);
                       -a4/A4*sqrt(2*g*y(4))+(1-gamma1)*k1/A4*u1(1)];
                       
    [t,x] = ode45(f4tanques,[0 Ts],h1,opts);
    h = x(end,:)';
end

