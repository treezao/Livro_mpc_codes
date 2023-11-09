%%% Exemplo duplo integrador do capítulo GPC SISO
%%% cálculo do ganho irrestrito
%%% cálculo das funções de transferência equivalentes.

clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% Modelo da planta em função de transferência
A = [1 -1.905, 0.905];
num = [0.09675 0.09358]; 
B = 4*num;
Bq = -[0 num];

n = 1; % número de saídas
m = 1; % número de manipuladas
mq = 1; % número de perturbações

%% obtenção do modelo em espaço de estados
Gz = [tf([0 B],A,1), tf(Bq,A,1)] % modelo em função de transferência

sys = ss(Gz)

%% Sintonia do SSMPC

N1 = 1; %% horizonte de predição inicial
N2 = 30; %% horizonte predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; %% horizonte de controle
Nq = 3; %% horizonte de perturbação

lambda = 1; % ponderação do incremento de controle para o GPC1
delta = 1;  % ponderação dos erros futuros

af = 0; % polo do filtro de referência


Umax = 100; % limite superior do sinal de controle
Umin = 0; % limite inferior do sinal de controle
dumax = 3; % incremento máximo do sinal de controle
dumin = -dumax; % incremento mínimo do sinal de controle

%% obtenção das matrizes do SSMPC
A = sys.A;
B = sys.B(:,1);
Bq = sys.B(:,2);
C = sys.C;
Cq = sys.D(:,2);

temp1 = 0;

F = [];
Ii = [];
for i=1:N
    temp1 = temp1+A^i;
    
    Ftemp = [C*temp1 -C*temp1];
    F = [F;Ftemp];
    
    Ii = [Ii; eye(n)];
end

Fq = [C*Bq+Cq, -C*Bq];
temp1 = eye(size(A));
for i=1:N-1
    temp1 = temp1+A^i;
    
    Fqtemp = [C*temp1*Bq+Cq, -C*temp1*Bq];
    Fq = [Fq;Fqtemp];
    
end



G0 = [];
temp1 = eye(size(A));
for i=1:N
    G0 = [G0;C*temp1*B];
    
    temp1 = temp1 + A^i;
end


G = G0;
for i= 1:Nu-1
    temp = [zeros(i*n,m);G0(1:end-i*n,:)];
    G = [G,temp];
end


G0q = [Cq];
temp1 = eye(size(A));
for i=1:N-1
    
    G0q = [G0q;C*temp1*Bq+Cq];
    
    temp1 = temp1 + A^i;
end

Gq = [];
for i=0:Nq-1
    temp = [zeros(i*n,mq);G0q(1:end-i*n,:)];
    Gq = [Gq,temp];
end


%%% matrizes de ponderações
Qei = diag(delta);

Qe = [];
for i=1:N
    Qe = blkdiag(Qe,Qei);
end

Qui = diag(lambda);

Qu = [];
for i=1:Nu
    Qu = blkdiag(Qu,Qui);
end

Kmpc = (G'*Qe*G+Qu)\(G'*Qe);
Kmpc1 = Kmpc(1:m,:);

Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;

%%% matrizes de restrições
T = tril(ones(Nu));

Rbar = [eye(Nu);
         -eye(Nu);
         T;
        -T];

rbar1 = [dumax*ones(Nu,1);
         -dumin*ones(Nu,1)];

%% vetores de simulação
nin = 10; % iteração inicial (para inicialização correta)
nit = nin+100; % número de iterações da simulação

refs = zeros(1,nit); % vetor das referências
perts = zeros(1,nit+Nq); % vetor da perturbação

refs(nin+40:end) = 30;
perts(1,nin+5:end) = 20;
perts(1,nin+70:end) = 30;

%% simulação do SSMPC - caso 1 - sem FF
estados = zeros(size(A,1),nit);
saidas = zeros(1,nit); % vetor da saída
entradas = zeros(1,nit); % vetor do sinal de controle
du = zeros(1,nit); % vetor dos incrementos de controle

rfant = 0;

for k=nin+1:nit
    %%% simulação do processo
    estados(:,k) = A*estados(:,k-1)+B*entradas(:,k-1)+Bq*perts(:,k-1);
    saidas(1,k) = C*estados(:,k);
              
    %%% vetor de referências
    rf = af*rfant + (1-af)*refs(k);
    rfant = rf;
    R = rf*ones(N,1);
    
    %%% cálculo da resposta livre;
    f = Ii*saidas(:,k)+ F*[estados(:,k);
                           estados(:,k-1)];
       
    %%% cálculo das matrizes de restrição
    rbar = [rbar1;
             (Umax-entradas(1,k-1))*ones(Nu,1);
             (-Umin+entradas(1,k-1))*ones(Nu,1)];
    
    fqp = fqp1*(R-f);
    %%% cálculo do incremento de controle ótimo
%     duOti= Kmpc1*(R-f);
    duOti = quadprog(Hqp,fqp,Rbar,rbar);
    du(1:m,k) = duOti(1:m);
    
    %%% cálculo do sinal de controle ótimo
    entradas(:,k) = du(:,k)+entradas(:,k-1);
    
    
              
end

%%% 
saidas1 = saidas;
estados1 = estados;
entradas1 = entradas;
du1 = du;










%% simulação do SSMPC - caso 2 - com FF mas Nq = 0
estados = zeros(size(A,1),nit);
saidas = zeros(1,nit); % vetor da saída
entradas = zeros(1,nit); % vetor do sinal de controle
du = zeros(1,nit); % vetor dos incrementos de controle

rfant = 0;

for k=nin+1:nit
    %%% simulação do processo
    estados(:,k) = A*estados(:,k-1)+B*entradas(:,k-1)+Bq*perts(:,k-1);
    saidas(1,k) = C*estados(:,k);
              
    %%% vetor de referências
    rf = af*rfant + (1-af)*refs(k);
    rfant = rf;
    R = rf*ones(N,1);
    
    %%% cálculo da resposta livre;
    f = Ii*saidas(:,k)...
         +F*[estados(:,k);
           estados(:,k-1)]...
        +Fq*[perts(:,k);
             perts(:,k-1)];
       
    %%% cálculo das matrizes de restrição
    rbar = [rbar1;
             (Umax-entradas(1,k-1))*ones(Nu,1);
             (-Umin+entradas(1,k-1))*ones(Nu,1)];
    
    fqp = fqp1*(R-f);
    %%% cálculo do incremento de controle ótimo
%     duOti= Kmpc1*(R-f);
    duOti = quadprog(Hqp,fqp,Rbar,rbar);
    du(1:m,k) = duOti(1:m);
    
    %%% cálculo do sinal de controle ótimo
    entradas(:,k) = du(:,k)+entradas(:,k-1);
    
    
              
end

saidas2 = saidas;
estados2 = estados;
entradas2 = entradas;
du2 = du;









%% simulação do SSMPC - caso 3 - com FF mas Nq = 5
estados = zeros(size(A,1),nit);
saidas = zeros(1,nit); % vetor da saída
entradas = zeros(1,nit); % vetor do sinal de controle
du = zeros(1,nit); % vetor dos incrementos de controle

rfant = 0;

for k=nin+1:nit
    %%% simulação do processo
    estados(:,k) = A*estados(:,k-1)+B*entradas(:,k-1)+Bq*perts(:,k-1);
    saidas(1,k) = C*estados(:,k);
              
    %%% vetor de referências
    rf = af*rfant + (1-af)*refs(k);
    rfant = rf;
    R = rf*ones(N,1);

    
    %%% cálculo da resposta livre;
    f = Ii*saidas(:,k)...
        +F*[estados(:,k);
           estados(:,k-1)]...
        +Fq*[perts(:,k);
             perts(:,k-1)]...
        +Gq*(perts(:,k+1:k+Nq)-perts(:,k:k+Nq-1))';
       
    %%% cálculo das matrizes de restrição
    rbar = [rbar1;
             (Umax-entradas(1,k-1))*ones(Nu,1);
             (-Umin+entradas(1,k-1))*ones(Nu,1)];
    
    fqp = fqp1*(R-f);
    %%% cálculo do incremento de controle ótimo
%     duOti= Kmpc1*(R-f);
    duOti = quadprog(Hqp,fqp,Rbar,rbar);
    du(1:m,k) = duOti(1:m);
    
    %%% cálculo do sinal de controle ótimo
    entradas(:,k) = du(:,k)+entradas(:,k-1);
    
    
              
end

saidas3 = saidas;
estados3 = estados;
entradas3 = entradas;
du3 = du;

% norm(saidas3(:,nin+60:nin+80)-refs(:,nin+60:nin+80))
% norm(saidas2(:,nin+60:nin+80)-refs(:,nin+60:nin+80))
% norm(saidas1(:,nin+60:nin+80)-refs(:,nin+60:nin+80))

%% Geração das figuras
cores = gray(5);
cores = cores(1:end-1,:);

hf = figure
h=subplot(3,1,1)
plot(saidas1(1,nin+1:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(refs(1,nin+1:nit),'-.','LineWidth',tamlinha,'Color',cores(4,:))
ylim([-15 35])
hl = legend('Saída','Referência','Location','SouthEast')
% hl.Position = [0.6810 0.7246 0.2375 0.1107];
h.YTick = [0 15 30]
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(entradas1(1,nin+1:nit),'LineWidth',tamlinha,'Color',cores(1,:))

ylabel('Manipulada','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);
% ylim([-5 5])

h = subplot(3,1,3)
plot(du1(1,nin+1:nit),'LineWidth',tamlinha,'Color',cores(1,:))

h.YTick = [-3 0 3]
ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (amostras)','FontSize', tamletra)
grid on
ylim([-4 4])
set(h, 'FontSize', tamletra);

% h=subplot(3,1,3)
% plot(t,erro(1:nit-N2-1),'LineWidth',tamlinha,'Color',cores(1,:))
% title('Erro', 'FontSize', tamletra);
% xlabel('tempo', 'FontSize', tamletra);
% set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.6845 0.7108 0.2071 0.1242];
% print('exemploMotorDC_ssmpc','-depsc')


%%
limx = [60 80]

hf = figure
h=subplot(3,1,1)
plot(saidas1(1,nin+1:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(saidas2(1,nin+1:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(saidas3(1,nin+1:nit),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(refs(1,nin+1:nit),'-.','LineWidth',tamlinha,'Color',cores(4,:))
% ylim([0 2])
xlim(limx)
hl = legend('s/ FF','c/ FF (N_q=0)','c/ FF (N_q=3)','Referência','Location','SouthWest')
% hl.Position = [0.1423 0.6451 0.2393 0.2298];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on
ylim([27 31])

h = subplot(3,1,2)
plot(entradas1(1,nin+1:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(entradas2(1,nin+1:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(entradas3(1,nin+1:nit),':','LineWidth',tamlinha,'Color',cores(3,:))
xlim(limx)

ylabel('Manipulada','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);


h = subplot(3,1,3)
plot(du1(1,nin+1:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(du2(1,nin+1:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(du3(1,nin+1:nit),':','LineWidth',tamlinha,'Color',cores(3,:))
xlim(limx)

ylim([-3 3])

ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (amostras)','FontSize', tamletra)
grid on
% ylim([-5 5])
set(h, 'FontSize', tamletra);

ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))


hf.Position = tamfigura;
hl.Position = [0.7530 0.5844 0.2089 0.2565];
% print('exemploMotorDC_ssmpc_comp','-depsc')






