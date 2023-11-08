clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%% CSTR instável
global Ac dH p cp ER k0 V
%%% constantes
Ac = 0.05;
dH = -50e3;
p = 1e3;
cp = 1;
ER = 1000;
k0 = 4.872997584;
V = 0.05;

Ts = 3;

m = 3; % número de entradas da planta - manipuladas
n = 3; % número de saidas da planta - controladas
mq = 0; % número de perturbações da planta

%%% pontos de operação
%%%% entradas
%%%%% manipuladas
qo0 = 5e-3;
Caf0 = 5;
Qh0 = 0.75;

%%%%% perturbações
qi0 = 5e-3;
Ti0 = 350;

%%%% saídas %% PO do artigo, estão incorretos
h0 = 1;
Ca0 = 1;
T0 = 400;

%%% teste ponto de operação
x0 = [h0;Ca0;T0];
u0 = [qo0;Caf0;Qh0];
q0 = [qi0;Ti0];

fR1 = @(T) k0*exp(-ER/T);
fcstrpo= @(y) [(q0(1)-u0(1))/Ac;
                q0(1)*(u0(2)-y(2))/V-y(2)*fR1(y(3));
                q0(1)*(q0(2)-y(3))/V-(dH/p/cp)*fR1(y(3))*y(2)-u0(3)/V];

x0f = fsolve(fcstrpo,x0);
[x0,x0f]
            
            
            
%% obtenção dos modelos degrau
modDegrauU = {};
modDegrauQ = {};
%%%%%%%%%
for j = 1:m
    xs0 = x0;
    v = [];

    u00 = u0;
    u00(j) = 0.98*u0(j);
    for ll = 1:300
        xs0 = modeloCSTR(xs0,u00,q0,Ts);
        v = [v,xs0];
    end
    
    for i=1:n
        modDegrauU{i,j} = (v(i,:)'-x0(i))/(u00(j)-u0(j));
    end
end

for j = 1:2
    xs0 = x0;
    v = [];

    q00 = q0;
    q00(j) = 0.98*q0(j);
    for ll = 1:300
        xs0 = modeloCSTR(xs0,u0,q00,Ts);
        v = [v,xs0];
    end
    
    for i=1:n
        modDegrauQ{i,j} = (v(i,:)'-x0(i))/(q00(j)-q0(j));
    end
end


%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [5 5 5]; %horizontes de controle - dimensão 1 x m
N1 = [1 1 1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [30 30 30]; % fim dos horizontes de predição - dimensão 1 x n
N = N2-N1+1; %horizonte de predição - não editar

Nss = [30 30 30];
Nf = [30 30 30];

delta = [1 1 1]./N./(x0.^2)'; %ponderação nos erros - dimensão 1 x n
lambda = [1 1 1]./Nu./(u0.^2)'; %ponderação nas ações de controle - dimensão 1 x m       
betaf = 0.8*[1 1 1]; %% polos dos filtros dos erros de predição por saída

umin = [0, u0(2)-0.4,u0(3)-0.3];
umax = [2*u0(1), u0(2)+0.4,u0(3)+0.3];


%% definição do cenário de simulação
nit =floor(900/Ts); %tempo de simulação
nin = max(max(Nf),max(Nss))+10; %número da iteração inicial
nit = nit+nin;


refs = repmat(x0,[1,nit]);
perts = repmat(q0,[1,nit]);

refs(1,nin+floor(12/Ts):end) = x0(1)*1.05;
refs(2,nin+floor(100/Ts):end) = x0(2)*1.05;
refs(2,nin+floor(200/Ts):end) = x0(2)*1.1;
refs(2,nin+floor(300/Ts):end) = x0(2)*1;
refs(2,nin+floor(400/Ts):end) = x0(2)*0.95;
refs(2,nin+floor(500/Ts):end) = x0(2);


perts(1,nin+floor(700/Ts):end) = q0(1)*1.05;
perts(2,nin+floor(800/Ts):end) = q0(2)*0.95;
% 

saidas = repmat(x0,[1,nit]); % vetor da saída
entradas = repmat(u0,[1,nit]); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle



%% Controlador GDMC - MIMO


Qe = diag(repelem(delta,1,N)); %montagem da matriz de ponderação do erro
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

Kdmc = (G'*Qe*G+Qu)\G'*Qe; % utilizado apenas no caso irrestrito
Kdmc1 = [];
for i=1:m
    Kdmc1 = [Kdmc1;Kdmc(sum(Nu(1:i-1))+1,:)];
end

%%% matrizes do problema de otimização
Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;

Rbar = [];
for i=1:m
    Rbar = blkdiag(Rbar,tril(ones(Nu(i))));
end
Rbar = [Rbar;-Rbar];

%% Cálculo do filtro do erro de predição
% formato do filtro (az^2+bz)/(z-betaf)^2

z = tf('z',Ts);
nf = 2;



%%% para a saída 1 é especial por conta de ser integrador.
for l=1:1
    F{l} = tf(0,1,Ts);
    
    for i=N1(l):N2(l)
        %%% monta sistema para obtenção dos parâmetros do filtro
        %%% primeira equação (z^i - F_i(z)) = 0 -> para z=pz
        %%% segunda equação F_i(1) = 0 -> força ganho unitário
        indf = i-N1(1)+1;
        Af = [1 1;
              -2*betaf(l),(-1-betaf(l))];
        bf = [(1-betaf(l))^2;
              i*(1-betaf(l))^3];
        X = Af\bf;
        F{l}(indf,1) = (X(1)*z^2+X(2)*z)/(z-betaf(l))^2;
        %%% teste da condição
%         zero(z^i - F{l}(indf,1))

        %%% armazena coeficientes gtil
        for k=1:m
            modDegrauUF{l,k,i} = filter(F{l}(indf,1).num{1},F{l}(indf,1).den{1},modDegrauU{l,k});    
        end
    end
end


%%% sem filtro para as outras saídas
for l=2:n
    F{l} = tf(0,1,Ts);
    
    for i=N1(l):N2(l)
%         %%% monta sistema para obtenção dos parâmetros do filtro
%         %%% primeira equação (z^i - F_i(z)) = 0 -> para z=pz
%         %%% segunda equação F_i(1) = 0 -> força ganho unitário
        indf = i-N1(1)+1;

        F{l}(indf,1) = z^2/z^2;
        %%% armazena coeficientes gtil
        for k=1:m
            modDegrauUF{l,k,i} = filter(F{l}(indf,1).num{1},F{l}(indf,1).den{1},modDegrauU{l,k});    
        end
    end
end

%%% armazena coeficientes gtil
for l=1:n
    for k=1:m
        %%% calcula a matriz H para o cálculo da resposta livre no caso GDMC
        H1 = [];
        H2 = [];

        for j=N1(1):N2(1)
            H1 = [H1;modDegrauU{l,k}(j+1:j+Nf(l))'];
            H2 = [H2;modDegrauUF{l,k,j}(1:Nf(l))'];

        end
        H{l,k} = H1-H2;

    end
end


%% inicializacao dos estados e variaveis da simulação
for l=1:n
    yfilt{l} = x0(l)*ones(N(l),nit);
end

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    saidas(:,k) = modeloCSTR(saidas(:,k-1),entradas(:,k-1),perts(:,k-1),Ts);
    
    %% -- Controlador GDMC 
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),N(i),1)];
    end 
    
    %%% calculo da resposta livre
    f = [];
    for l=1:n
        for i=1:N(1)
            yfilt{l}(i,k) = -F{l}(i,1).den{1}(2:end)*yfilt{l}(i,k-1:-1:k-nf)' + F{l}(i,1).num{1}*saidas(l,k:-1:k-nf)';
        end
        
        ftemp = 0;
        for ll=1:m
            ftemp = ftemp + H{l,ll}*du(ll,k-1:-1:k-Nf(ll))';
        end
        f = [f; ftemp + yfilt{l}(:,k)];
    end
        
    
    %%% montagem das matrizes de restrição
    rbar = repelem(umax'-entradas(:,k-1),Nu');
    rbar = [rbar;
            repelem(-umin'+entradas(:,k-1),Nu')];
    
    fqp = fqp1*(R-f);
    X = quadprog(Hqp,fqp,Rbar,rbar);
    for i=1:m
        du(i,k) = X(sum(Nu(1:i-1))+1,1);
    end
    %% Resolve o problema de otimização
%     du(:,k) = Kdmc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end



%% figuras
cores = gray(5);
cores = cores(1:end-1,:);

nm = 30;

tm = (1:nm)*Ts;

hf = figure
h=subplot(3,3,1)
plot(tm,modDegrauU{1,1}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on        
ylabel('\Delta h (m)','FontSize', tamletra)
title('\Delta q_o','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,3,2)
plot(tm,modDegrauU{1,2}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
title('\Delta C_{af}','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,3,3)
plot(tm,modDegrauU{1,3}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
title('\Delta Q_h/(\rho c_p)','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,3,4)
plot(tm,modDegrauU{2,1}(1:nm)*0,'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('\Delta C_a (kmol/m^3)','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,3,5)
plot(tm,modDegrauU{2,2}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
set(h, 'FontSize', tamletra);

h=subplot(3,3,6)
plot(tm,modDegrauU{2,3}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
set(h, 'FontSize', tamletra);

h=subplot(3,3,7)
plot(tm,modDegrauU{3,1}(1:nm)*0,'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('\Delta T (K)','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,3,8)
plot(tm,modDegrauU{3,2}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
xlabel('Tempo (segundos)','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,3,9)
plot(tm,modDegrauU{3,3}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on

set(h, 'FontSize', tamletra);

hf.Position = [680 558 560 420];

% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))

% print('cstrTito_integrador_gdmc_novo_modelo','-depsc')

%%
nm = 30;

tm = (1:nm)*Ts;

hf = figure
h=subplot(3,2,1)
plot(tm,modDegrauQ{1,1}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on        
ylabel('\Delta h (m)','FontSize', tamletra)
title('\Delta q_i','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,2,2)
plot(tm,modDegrauQ{1,2}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
title('\Delta T_i','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,2,3)
plot(tm,modDegrauQ{2,1}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('\Delta C_a (kmol/m^3)','FontSize', tamletra)
set(h, 'FontSize', tamletra);

h=subplot(3,2,4)
plot(tm,modDegrauQ{2,2}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
set(h, 'FontSize', tamletra);


h=subplot(3,2,5)
plot(tm,modDegrauQ{3,1}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('\Delta T (K)','FontSize', tamletra)
set(h, 'FontSize', tamletra);
xlabel('Tempo (segundos)','FontSize', tamletra)

h=subplot(3,2,6)
plot(tm,modDegrauQ{3,2}(1:nm),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
xlabel('Tempo (segundos)','FontSize', tamletra)
set(h, 'FontSize', tamletra);

hf.Position = [680 558 560 420];

% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))

% print('cstrTito_integrador_gdmc_novo_modeloQ','-depsc')

%%
t = ((nin:nit)-nin)*Ts;
ind = nin:nit;


hf = figure
h=subplot(3,2,1)
plot(t,saidas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
grid on    
title('Controladas','FontSize', tamletra)
ylabel('h (m)','FontSize', tamletra)
xlim([0 600])
hl = legend('y_i','Referência','Location','SouthEast')
set(h, 'FontSize', tamletra);


h=subplot(3,2,3)
plot(t,saidas(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
hold on
plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylabel('C_{a} (kmol/m^3)','FontSize', tamletra)
xlim([0 600])
set(h, 'FontSize', tamletra);


h=subplot(3,2,5)
plot(t,saidas(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
hold on
plot(t,refs(3,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylabel('T (K)','FontSize', tamletra)
xlabel('Tempo (segundos)','FontSize', tamletra)
xlim([0 600])
set(h, 'FontSize', tamletra);



h=subplot(3,2,2)
plot(t,entradas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('q_o (m^3/s)','FontSize', tamletra)
xlim([0 600])
title('Manipuladas','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);


h2=subplot(3,2,4)
plot(t,entradas(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('C_{af} (kmol/m^3)','FontSize', tamletra)
xlim([0 600])
h2.YTick = [5 5.4];
set(h2,'YAxisLocation','right')
set(h2, 'FontSize', tamletra);


h=subplot(3,2,6)
plot(t,entradas(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('Q_h/(\rho c_p) (Km^3/s)','FontSize', tamletra)
xlim([0 600])
xlabel('Tempo (segundos)','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);




% print('cstrTito_integrador_gdmc_novo_ref','-depsc')

%%
hf = figure
h=subplot(3,2,1)
plot(t,saidas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
grid on    
title('Controladas','FontSize', tamletra)
ylabel('h (m)','FontSize', tamletra)
xlim([600 900])
hl = legend('y_i','Referência','Location','NorthEast')
set(h, 'FontSize', tamletra);
ylim([1.045 1.115])

h=subplot(3,2,3)
plot(t,saidas(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
hold on
plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylabel('C_{a} (kmol/m^3)','FontSize', tamletra)
xlim([600 900])
set(h, 'FontSize', tamletra);
ylim([0.98 1.07])

h=subplot(3,2,5)
plot(t,saidas(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
hold on
plot(t,refs(3,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
ylabel('T (K)','FontSize', tamletra)
xlabel('Tempo (segundos)','FontSize', tamletra)
xlim([600 900])
set(h, 'FontSize', tamletra);
ylim([385 405])


h=subplot(3,2,2)
plot(t,entradas(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('q_o (m^3/s)','FontSize', tamletra)
xlim([600 900])
title('Manipuladas','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);
ylim([0.00495 0.00535])

h2=subplot(3,2,4)
plot(t,entradas(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('C_{af} (kmol/m^3)','FontSize', tamletra)
xlim([600 900])
% h2.YTick = [5 5.4];
set(h2,'YAxisLocation','right')
set(h2, 'FontSize', tamletra);
ylim([4.7 5.02])

h=subplot(3,2,6)
plot(t,entradas(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
ylabel('Q_h/(\rho c_p) (Km^3/s)','FontSize', tamletra)
xlim([600 900])
xlabel('Tempo (segundos)','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);
ylim([0.63 0.76])



% hf.Position = [680 558 560 420];

% hf.Position = tamfigura;
hl.Position = [0.3350 0.8052 0.2054 0.1060];
% hl2.Position = [0.8342 0.1739 0.1225 0.2952];

% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))

% print('cstrTito_integrador_gdmc_novo_pert','-depsc')


%% modelo do sistema 4 tanques
%%% x = [h,Ca,T]
%%% u = [qo,Caf,Qh]
%%% q = [qi,Ti]
function h = modeloCSTR(x,u,q,Ts)
    
    global Ac dH p cp ER k0 V
    
    opts = odeset('NonNegative',[1 2 3]);
    fR1 = @(T) k0*exp(-ER/T);
    
    fcstr= @(t,y) [(q(1)-u(1))/Ac;
                    q(1)*(u(2)-y(2))/V-y(2)*fR1(y(3));
                    q(1)*(q(2)-y(3))/V-(dH/p/cp)*fR1(y(3))*y(2)-u(3)/V];
                       
    [t,xf] = ode45(fcstr,[0 Ts],x,opts);
    h = xf(end,:)';
end
