%%% Exemplo duplo integrador do capítulo GPC SISO
%%% cálculo do ganho irrestrito
%%% cálculo das funções de transferência equivalentes.


clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%% Modelo da planta
s = tf('s');
G = [12.8/(16.7*s+1), -18.9/(21*s+1);
     6.6/(10.9*s+1), -19.4/(14.4*s+1)];
Gq = [3.8/(14.9*s+1);
      4.9/(13.2*s+1)];
  
Gz = c2d(G,1,'zoh');
Gqz = c2d(Gq,1,'zoh');

%%% adicionando os atrasos discretos manualmente
z  = tf('z',1)
Gz=  Gz.*[z^-1, z^-3;
          z^-7, z^-3;];
Gqz = Gqz.*[z^-8;
            z^-4];
        
sys = minreal(ss([Gz,Gqz]))

n = 2;
m = 2;
mq = 1;

A = sys.A;
B = sys.B(:,1:m);
Bq = sys.B(:,m+1:end);
C = sys.C;
Cq = sys.D(:,m+1:end);


        

%% Sintonia do SSMPC

N1 = [2,3]; %% horizonte de predição inicial
N2 = [40,30]; %% horizonte predição
N = N2-N1+1;
Nu = [5,5]; %% horizonte de controle
Nq = 0; %% horizonte de perturbação

lambda = [10,10]./Nu; % ponderação do incremento de controle
delta = [1,1]./N;  % ponderação dos erros futuros

af = 0; % polo do filtro de referência

Umax = [1,1];
Umin = [-1,-1];
dumax = [.1,.1];
dumin = -dumax;

%% obtenção das matrizes do SSMPC

%%% montando matrizes desconsiderando as diferenças de horizontes, assim,
%%% utiliza-se o valor máximo de N2 e Nu.
temp1 = 0;

F = [];
Ii = [];
for i=1:max(N2)
    temp1 = temp1+A^i;
    
    Ftemp = [C*temp1 -C*temp1];
    F = [F;Ftemp];
    Ii = [Ii; eye(n)];
end

Fq = [C*Bq+Cq, -C*Bq];
temp1 = eye(size(A));
for i=1:max(N2)-1
    temp1 = temp1+A^i;
    
    Fqtemp = [C*temp1*Bq+Cq, -C*temp1*Bq];
    Fq = [Fq;Fqtemp];
    
end



%%% G    
G0 = [];
temp1 = eye(size(A));
for i=1:max(N2)
    
    G0 = [G0;C*temp1*B];
    
    temp1 = temp1 + A^i;
end

G = G0;
for i= 1:max(Nu)-1
    temp = [zeros(i*n,m);G0(1:end-i*n,:)];
    G = [G,temp];
end


%%% Gq
G0q = [Cq];
temp1 = eye(size(A));
for i=1:max(N2)-1
    
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
for i=1:max(N2)
    Qe = blkdiag(Qe,Qei);
end

Qui = diag(lambda);

Qu = [];
for i=1:max(Nu)
    Qu = blkdiag(Qu,Qui);
end


%%% Selecionando as parcelas das matrizes de interesse de acordo com os
%%% horizontes estabelecidos na sintonia.

indc = [];
for i=1:m
%     indc = [indc, (1:Nu(i))*m-1+(i-1)];
    indc = [indc, i:m:(Nu(i)-1)*m+i];
end

indc = sort(indc);

indcq = [];
for i=1:mq
    indcq = [indcq, i:mq:(Nq(i)-1)*mq+i];
%     indcq = [indcq, (1:Nq(i))*mq-1+(i-1)];
end

indcq = sort(indcq);


indl = [];
for i=1:n
%     indl = [indl, (N1(i):N2(i))*n-1+(i-1)];
    indl = [indl, (N1(i)-1)*n+i:n:(N2(i)-1)*n+i];
end

indl = sort(indl);


F = F(indl,:);
Ii = Ii(indl,:);
Fq = Fq(indl,:);
G = G(indl,indc);

if(max(Nq)>0)
    Gq = Gq(indl,indcq);
end

Qu = Qu(indc,indc);
Qe = Qe(indl,indl);

Kmpc = (G'*Qe*G+Qu)\(G'*Qe);
Kmpc1 = Kmpc(1:m,:);

Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;

%%% matrizes de restrições
LB = repmat(dumin',max(Nu),1);
LB = LB(indc,:);
UB = repmat(dumax',max(Nu),1);
UB = UB(indc,:);

I = eye(m);
T = [];
for i=1:max(Nu)
    T = [T;
          repmat(I,1,i),zeros(m,max(Nu)*m-i*m)];
end

T = T(indc,indc);

Rbar = [T;
        -T];


%% vetores de simulação
nin = 10; % iteração inicial (para inicialização correta)
nit = nin+150; % número de iterações da simulação

refs = zeros(n,nit); % vetor das referências
perts = zeros(mq,nit); % vetor da perturbação

refs(1,nin+10:end) = 1;
refs(2,nin+50:end)=0.5;

perts(1,nin+100:end) = -0.4;

%% simulação do SSMPC
estados = zeros(size(A,1),nit);
saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

rfant = zeros(n,1);

for k=nin+1:nit
    %%% simulação do processo
    estados(:,k) = A*estados(:,k-1)+B*entradas(:,k-1)+Bq*perts(:,k-1);
    saidas(:,k) = C*estados(:,k);
              
    %%% vetor de referências
    rf = af*rfant + (1-af)*refs(:,k);
    rfant = rf;
    R = repmat(rf,max(N2),1);
    R = R(indl,:);
    
    %%% cálculo da resposta livre;
    f = Ii*saidas(:,k)...
        +F*[estados(:,k);
           estados(:,k-1)];
       
    %%% cálculo das matrizes de restrição
    rbarmax = repmat([Umax'-entradas(:,k-1)],max(Nu),1);
    rbarmin = repmat([-Umin'+entradas(:,k-1)],max(Nu),1);
    rbar = [rbarmax(indc,:);
             rbarmin(indc,:)];
    
    fqp = fqp1*(R-f);
         
    %%% cálculo do incremento de controle ótimo
%     duOti= Kmpc1*(R-f);
    duOti = quadprog(Hqp,fqp,Rbar,rbar,[],[],LB,UB);
    du(1:m,k) = duOti(1:m);
    
    %%% cálculo do sinal de controle ótimo
    entradas(:,k) = du(:,k)+entradas(:,k-1);
    
    
              
end


%% Geração das figuras
cores = gray(4);
cores = cores(1:end-1,:);

hf = figure
h=subplot(3,1,1)
plot(saidas(1,nin+1:nit)','LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(saidas(2,nin+1:nit)','--','LineWidth',tamlinha,'Color',cores(2,:))
plot(refs(:,nin+1:nit)','-.','LineWidth',tamlinha,'Color',cores(3,:))
% ylim([0 2])
ylabel('Controladas','FontSize', tamletra)
hl0 = legend('y_1','y_2','Ref.','Location','NorthEast')


set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(entradas(1,nin+1:nit)','LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(entradas(2,nin+1:nit)','--','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Manipuladas','FontSize', tamletra)
hl1 = legend('u_1','u_2','Location','NorthEast')

grid on
set(h, 'FontSize', tamletra);
% ylim([-5 5])

h = subplot(3,1,3)
plot(du(1,nin+1:nit)','LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(du(2,nin+1:nit)','--','LineWidth',tamlinha,'Color',cores(2,:))


ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (segundos)','FontSize', tamletra)
grid on
% ylim([-5 5])
set(h, 'FontSize', tamletra);

% h=subplot(3,1,3)
% plot(t,erro(1:nit-N2-1),'LineWidth',tamlinha,'Color',cores(1,:))
% title('Erro', 'FontSize', tamletra);
% xlabel('tempo', 'FontSize', tamletra);
% set(h, 'FontSize', tamletra);



hf.Position = tamfigura3;
hl0.Position = [0.8220 0.6325 0.1357 0.1619];
hl1.Position = [0.8274 0.4933 0.1161 0.1202];
% print('exemploColunaMentanol','-depsc')






