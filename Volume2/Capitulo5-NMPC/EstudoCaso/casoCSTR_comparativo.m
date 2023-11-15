%%% Simulação comparativa do PNMPC, NEPSAC e xxx
clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

import casadi.*
%% parametros CSTR completo
global ER V k0 k1 k2 k3 Ts optsode
%%% condicoes iniciais
Ca0 = 0.1 ; %[mol/l]
T0 = 438.5444; %[K]

f0 = 100; %[l/min] process flow rate
fc0 = 103.41; % [l/min] coolant flow rate

CAf0 = 1; %[mol/l] feed concentration
Tcf0 = 350.0030; %[K] inlet coolant temperature
Tf0 = 350; %[K] feed temperature

Ts = 0.1; % [min]

%%% parametros
ha = 7e5; %[cal/(min K)] heat transfer term
ER = 1e4; %[K] activation energy term
rho = 1e3; %[g/l] liquid densities rho=rho_c
rhoc = rho;
V = 100; %[l] CSTR volume
k0 = 7.2e10; %[1/min] reaction rate constant
dH = -2e5; %[cal/mol] heat of reaction
Cp = 1; %[cal/(g K)] specific heats
Cpc = 1;

k1 = -dH*k0/(rho*Cp);
k2 = rhoc*Cpc/(rho*Cp*V);
k3 = ha/(rhoc*Cpc);


%%% opções para o ODE45 que simulará o sistema contínuo do CSTR
% opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
% opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
optsode = odeset('NonNegative',[1 2]);

%% restrições
umax = [150, 130];
umin = [75, 60];

dumax = [2,2];
dumin = -dumax;


%% definicao das variaveis de simulacao
%%% Modelo CSTR listado em DU et al 
n = 2; % numero de saidas
m = 2; % numero de entradas
mq = 3; % numero de perturbacoes

na = 1; % numero de amostras de y atrasadas usadas no calculo da saida
nb = 1; % numero de amostras de u atrasadas usadas no calculo da saida
nq = 1; % numero de amostras de q atrasadas usadas no calculo da saida



%% Sintonia dos controladores

%%% horizontes
N1 = [1,1]
N2 = [50, 50]; 
N2max = max(N2);
N = N2-N1+1;
    
Nu = [30, 30]; % vetor contendo os horizontes de controle

delta = [1/Ca0^2,100/T0^2]./N;
lambda = [1/f0^2,1/fc0^2]./Nu;

alpha = 0.85; % polo do filtro de referência

C = 1; % polinômio C do modelo da perturbação
D = [1 -1]; % polinômio D do modelo da perturbação

nc = size(C,2);
nd = size(D,2);

Qe = [];
for i=1:n
    Qe = blkdiag(Qe,delta(i)*eye(N2(i)-N1(i)+1));
end

Qu = [];
for i=1:m
    Qu = blkdiag(Qu,lambda(i)*eye(Nu(i)));
end



%% configuracaoes de simulacao
%%% inicializacao de algumas variaveis
nin = 2; %%% iteracao inicial
nit = nin+430; %%% iteracao final

saida = zeros(n,nit);
saida(1,1:nin) = Ca0;
saida(2,1:nin) = T0;

entrada = zeros(m,nit);
entrada(1,:) = f0;
entrada(2,:) = fc0;

refs = zeros(n,nit+max(N2));
refs(1,:) = Ca0;
refs(2,:) = T0;
refs(1,nin+10:end) = 0.08;
refs(2,nin+70:end) = T0+5;
refs(1,nin+130:end) = 0.1;
refs(2,nin+190:end) = T0;


perts = zeros(mq,nit);
perts(1,:) = CAf0;
perts(2,:) = Tf0;
perts(3,:) = Tcf0;

perts(1,nin+250:nit) = CAf0*1.03;
perts(2,nin+310:nit) = Tf0+3;
perts(3,nin+370:nit) = Tcf0-3;


eta = zeros(n,nit+max(N2));
e = zeros(n,nit+max(N2));
saidaModelo = zeros(n,nit);
du = zeros(m,nit);

rfant = saida(:,1);

%% configuração do CASADI para a resoluação do problema de otimização para o PNMPC
Hs = MX.sym('H',sum(Nu),sum(Nu));
As = MX.sym('g',sum(Nu),sum(Nu));

%%% ipopt
prob = struct('h', Hs.sparsity(), 'a', As.sparsity())
opts = struct('print_time',false,'nlpsol','ipopt','verbose',false);
opts.nlpsol_options.ipopt.print_level = 0;
QPsolver = conic('QPsolver','nlpsol',prob,opts);


%%% qpoases
% prob = struct('h', Hs.sparsity(), 'a', As.sparsity())
% opts = struct('print_time',false,'printLevel','none','verbose',false);
% QPsolver = conic('QPsolver','qpoases',prob,opts);


%% matrizes para restrições com o PNMPC
LBdu = [];
UBdu = [];
rbarLB0 = [];
rbarUB0 = [];
Rbar = [];
for i=1:2
    LBdu = [LBdu;ones(Nu(i),1)*dumin(i)];
    UBdu = [UBdu;ones(Nu(i),1)*dumax(i)];
    
    T = tril(ones(Nu(i)));
    Rbar = blkdiag(Rbar,T);
    
    rbarLB0 = [rbarLB0;ones(Nu(i),1)*umin(i)];
    rbarUB0 = [rbarUB0;ones(Nu(i),1)*umax(i)];

end

X = zeros(sum(Nu),1);

%% simulação 

for k=nin:nit
    %% modelo do CSTR
    saida(:,k) = modeloCSTR(saida(:,k-1),entrada(:,k-1),perts(:,k-1));

    %% controlador PNMPC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(:,k) = modeloCSTRnominal(saida(:,k-1),entrada(:,k-1),[CAf0;Tf0;Tcf0]);

    %%% calculo dos etas futuros
    eta(:,k) = saida(:,k) - saidaModelo(:,k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:max(N2)
        eta(:,k+j) = -D(2:end)*eta(:,k+j-1:-1:k+j-nd+1)' + C*e(:,k+j:-1:k+j-nc+1)';
    end
    
    
    %%% cálculo da resposta livre corrigida
    yc = zeros(n,N2max);
    yc(:,1) = modeloCSTRnominal(saida(:,k),entrada(:,k-1),[CAf0;Tf0;Tcf0]) + eta(:,k+1);
    for i=2:N2max
        yc(:,i) = modeloCSTRnominal(yc(:,i-1),entrada(:,k-1),[CAf0;Tf0;Tcf0]) + eta(:,k+i);
    end
    
    %%% passando para o formato vetorial correto
    f = zeros(sum(N2-N1)+1,1);
    for i=1:n
        f(sum(N(1:i-1))+1:sum(N(1:i)),1) = yc(i,N1(i):N2(i))';
    end

    
    %%% calculo de G
    % calculo de y0
    y0 = zeros(n,N2max);
    y0(:,1) = modeloCSTRnominal(saida(:,k),entrada(:,k-1),[CAf0;Tf0;Tcf0]);
    for i=2:N2max
        y0(:,i) = modeloCSTRnominal(y0(:,i-1),entrada(:,k-1),[CAf0;Tf0;Tcf0]);
    end
    
    %%% passando para o formato vetorial correto
    y0v = zeros(sum(N2-N1)+1,1);
    for i=1:n
        y0v(sum(N(1:i-1))+1:sum(N(1:i)),1) = y0(i,N1(i):N2(i))';
    end

    %%% calculo de yi para obter G
    G = zeros(sum(N),sum(Nu));
    for i=1:m
        u0 = entrada(:,k-1);
        
        if(entrada(i,k-1) == 0)
            epsilon = 0.001;
        else
            epsilon = entrada(i,k-1)/1000;
        end

        for j=1:Nu
            yi = zeros(n,N2max);

            ui = repmat(u0,[1 N2max]);
            ui(i,j:end) = ui(i,j:end)+epsilon;

            yi(:,1) = modeloCSTRnominal(saida(:,k),ui(:,1),[CAf0;Tf0;Tcf0]);
            for jj=2:N2max
                yi(:,jj) = modeloCSTRnominal(yi(:,jj-1),ui(:,jj),[CAf0;Tf0;Tcf0]);
            end
            
            %%% passando para o formato vetorial correto
            yiv = zeros(sum(N2-N1+1),1);
            for jj=1:n
                yiv(sum(N(1:jj-1))+1:sum(N(1:jj)),1) = yi(jj,N1(jj):N2(jj))';
            end
            
            G(:,sum(Nu(1:i-1))+j) = (yiv-y0v)/epsilon;
            
        end
    end

    %%% referências futuras
    rf = alpha*rfant + (1-alpha)*refs(:,k);
    rfant = rf;
    rfut = [];
    for i=1:n
        rfut = [rfut;
                rf(i,1)*ones(N2(i)-N1(i)+1,1)];
    end
    
    %%% montagem do problema de otimização
    Hqp = 2*(G'*Qe*G+Qu);
    fqp = -2*G'*Qe*(rfut-f);
    
    temp1 = [];
    temp2 = [];
    for i=1:2
        temp1 = [temp1;-entrada(i,k-1)*ones(Nu(i),1)];
        temp2 = [temp2;-entrada(i,k-1)*ones(Nu(i),1)];
    end
    
    rbarLB = rbarLB0 + temp1;
    rbarUB = rbarUB0 + temp2;
    
    tic
    x_opt = QPsolver('h', Hqp, 'g', fqp,'x0',X,'lbx', LBdu, 'ubx', UBdu,'a', Rbar, 'lba',rbarLB,'uba',rbarUB);
    x2(k) = toc;
    
    X = full(x_opt.x);

    for i=1:m
        ind = sum(Nu(1:i-1));
        indf = sum(Nu(1:i));
        du(i,k) = X(ind+1);
    end
        
    entrada(:,k) = entrada(:,k-1)+du(:,k); 
    

end

saidaPNMPC = saida;
entradaPNMPC = entrada;
duPNMPC = du;








%% NEPSAC - inicialização e configuração

%%% matrizes de restrição
T = [];
for i=1:m
    T = blkdiag(T,tril(ones(Nu(i))));
end
T1 = inv(T);



%%% reinicialização das variáveis de simulação
saida = zeros(n,nit);
saida(1,1:nin) = Ca0;
saida(2,1:nin) = T0;

entrada = zeros(m,nit);
entrada(1,:) = f0;
entrada(2,:) = fc0;

eta = zeros(n,nit+max(N2));
e = zeros(n,nit+max(N2));
saidaModelo = zeros(n,nit);
du = zeros(m,nit);

rfant = saida(:,1);

%% matrizes para restrições com o NEPSAC
LBu0 = [];
UBu0 = [];
rbarLB0 = [];
rbarUB0 = [];
for i=1:2
    LBu0 = [LBu0;ones(Nu(i),1)*umin(i)];
    UBu0 = [UBu0;ones(Nu(i),1)*umax(i)];
        
    rbarLB0 = [rbarLB0;ones(Nu(i),1)*dumin(i)];
    rbarUB0 = [rbarUB0;ones(Nu(i),1)*dumax(i)];

end

ub = repmat(entrada(:,k-1),[1 N2max]); % u base inicial

%% simulação 

for k=nin:nit
    %% modelo do CSTR
    saida(:,k) = modeloCSTR(saida(:,k-1),entrada(:,k-1),perts(:,k-1));

    
    %% controlador NEPSAC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(:,k) = modeloCSTRnominal(saida(:,k-1),entrada(:,k-1),[CAf0;Tf0;Tcf0]);
    
    
    %%% calculo dos etas futuros
    eta(:,k) = saida(:,k) - saidaModelo(:,k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:max(N2)
        eta(:,k+j) = -D(2:end)*eta(:,k+j-1:-1:k+j-nd+1)' + C*e(:,k+j:-1:k+j-nc+1)';
    end
    
    %%% referências futuras
    rf = alpha*rfant + (1-alpha)*refs(:,k);
    rfant = rf;
    rfut = [];
    for i=1:n
        rfut = [rfut;
                rf(i,1)*ones(N2(i)-N1(i)+1,1)];
    end
    
    ub = repmat(entrada(:,k-1),[1 N2max]); % u base inicial
    dunorm = inf; % inicilizar norma de delta u como infinito
    
    cont = 0;
    contt = 0;
    u0 = entrada(:,k-1)/1000;
    while (dunorm > .1 && cont <= 10) % iterar até obter a resposta dentro do erro especificado
        %%% Obtenção da resposta base
        yb = zeros(n,N2max); % resposta base
        yb(:,1) =  modeloCSTRnominal(saida(:,k),entrada(:,k-1),[CAf0;Tf0;Tcf0]) + eta(:,k+1); % yb(k+1|k)
        for i=2:N2max
            yb(:,i) = modeloCSTRnominal(yb(:,i-1),entrada(:,k-1),[CAf0;Tf0;Tcf0]) + eta(:,k+i);
        end
        
        %%% passando para o formato vetorial correto
        ybv = zeros(sum(N2-N1)+1,1);
        for i=1:n
            ybv(sum(N(1:i-1))+1:sum(N(1:i)),1) = yb(i,N1(i):N2(i))';
        end
        
        
        %%% cálculo de Ge
        Ge = zeros(sum(N),sum(Nu));

        for jj=1:m
            for j=1:Nu(jj)
                yot = zeros(n,N2max);
                udelta = ub;
                if(j<Nu(jj))
                    udelta(jj,j) = udelta(jj,j)+u0(jj,1);
                else
                    udelta(jj,j:end) = udelta(jj,j:end)+u0(jj,1);
                end

                yot(:,1) = modeloCSTRnominal(saida(:,k),udelta(:,1),[CAf0;Tf0;Tcf0]) + eta(:,k+1);
                for i=2:N2max
                    yot(:,i) = modeloCSTRnominal(yot(:,i-1),udelta(:,i),[CAf0;Tf0;Tcf0]) + eta(:,k+i);
                end
                
                %%% passando para o formato vetorial correto
                yotv = zeros(sum(N2-N1)+1,1);
                for i=1:n
                    yotv(sum(N(1:i-1))+1:sum(N(1:i)),1) = yot(i,N1(i):N2(i))';
                end

                Ge(:,sum(Nu(1:jj-1))+j) = (yotv - ybv)/u0(jj);
                
            end
        end
        
        
        %%% passando para o formato vetorial correto
        ubv = zeros(sum(Nu),1);
        for i=1:m
            ubv(sum(Nu(1:i-1))+1:sum(Nu(1:i)),1) = ub(i,1:Nu(i))';
        end
        
        utemp = zeros(sum(Nu),1);
        for i=1:m
            utemp(sum(Nu(1:i-1))+1:sum(Nu(1:i)),1) = -entrada(i,k-1)*ones(Nu(i),1);
        end
        
        %%% montando as matrizes do problema de otimização
        temp = T1*(utemp+ubv);
        Hqp = 2*(Ge'*Qe*Ge+T1'*Qu*T1);
        fqp = 2*(Ge'*Qe*(ybv-rfut)+T1'*Qu*temp);
        
        rbarLB = rbarLB0 - temp;
        rbarUB = rbarUB0 - temp;
        LBu = LBu0-ubv;
        UBu = UBu0-ubv;
        
        tic
        x_opt = QPsolver('h', Hqp, 'g', fqp,'lbx', LBu, 'ubx', UBu,'a', T1, 'lba',rbarLB,'uba',rbarUB);
        contt = contt+ toc;
        X = full(x_opt.x);
        
        
        for i=1:m
            ind = sum(Nu(1:i-1))+1:sum(Nu(1:i)); 
            ub(i,1:Nu(i)) =  ub(i,1:Nu(i)) + X(ind)';
            ub(i,Nu(i)+1:end) = ub(i,Nu(i));
        end
        
        dunorm = norm(X);
        cont = cont+1;
    end
    nepsacCont(k) = cont; %% armazena número de iterações
    
    entrada(:,k) = ub(:,1);
    du(:,k) = entrada(:,k) - entrada(:,k-1);

    x2n(k) = contt; 
end

saidaNEPSAC = saida;
entradaNEPSAC = entrada;
duNEPSAC = du;












%% --- Simulação com NMPC
%%% configuração do CASADI
% função do modelo, para facilitar
f_modelo = @(x,u,q) [x(1) + Ts*(u(1)/V*(q(1)-x(1)) - k0*x(1)*exp(-ER/x(2)));
                     x(2) + Ts*(u(1)/V*(q(2)-x(2)) + k1*x(1)*exp(-ER/x(2)) + k2*u(2)*(1-exp(-k3/u(2)))*(q(3)-x(2)))];

                 
opti = casadi.Opti();   % Chamada da interface "Opti" de Casadi

f_custo = 0;    % Função objetivo

% Variáveis de otimização do NLP
X = opti.variable(n, max(N));    % Variável dos estados de k=1 até k=N
dU = opti.variable(m, max(Nu));   % Variável de control de k=0 até k=Nu

% Parâmetros do NLP
Xr = opti.parameter(n);     % Referência 
Xk = opti.parameter(n);     % Estado em k
Xk1 = opti.parameter(n);    % estado em k-1
Uk1 = opti.parameter(m);     % sinal de controle em k-1
Qk = opti.parameter(mq);     % Perturbação em k
Qk1 = opti.parameter(mq);    % perturbação em k-1


% montando a função custo e restrições
f_custo = 0;

% Montando as restrições e função custo
U = Uk1; % variável com o sinal de controle futuro
Xant2 = Xk1; % variável com o estado em k-1
Xant1 = Xk; % variável com o estado em k
Qant1 = Qk1; % variável com a perturbçaão em k-1
Qa = Qk; % perturbação em k 
for j=1:max(N)
    if(j<=max(Nu))
        Uant = U; % sinal de controle anterior
        U = U + dU(:,j);
        
        % Restrições de operação
        opti.subject_to(dumin' <= dU(:,j) <= dumax');
        opti.subject_to(umin' <= U <= umax');    
    else
        Uant = U;
    end

    %%% montando predições corrigidas
    X = f_modelo(Xant1,U,Qa) + Xant1 - f_modelo(Xant2,Uant,Qant1);

    %%% atuzliação das variáveis
    Xant2 = Xant1;
    Xant1 = X;
    
    Qant1 = Qa;
    Qa = Qk;

    % monta a função custo
    f_custo = f_custo + (Xr-X)'*diag(delta)*(Xr-X);
    
    if(j<=max(Nu))
        f_custo = f_custo + dU(:,j)'*diag(lambda)*dU(:,j);
    end

end                 

opti.minimize(f_custo);     % Declarando o objetivo: minimizar f_custo

%%% configurações do solver
solver_opts = struct('warm_start_init_point', 'yes');
casadi_opts = struct('print_time', false);
opti.solver('ipopt', casadi_opts, solver_opts);  % Construção do solver

%% config da simulação
saida = zeros(n,nit);
saida(1,1:nin) = Ca0;
saida(2,1:nin) = T0;

entrada = zeros(m,nit);
entrada(1,:) = f0;
entrada(2,:) = fc0;

eta = zeros(n,nit+max(N2));
e = zeros(n,nit+max(N2));
saidaModelo = zeros(n,nit);
du = zeros(m,nit);

rfant = saida(:,1);

du_oti_ant = zeros(m,max(Nu));
%% simulação 

for k=nin:nit
    %% modelo do CSTR
    saida(:,k) = modeloCSTR(saida(:,k-1),entrada(:,k-1),perts(:,k-1));

    %% controlador NMPC com CASADI
    %%% referências futuras
    rf = alpha*rfant + (1-alpha)*refs(:,k);
    rfant = rf;

    opti.set_value(Xr,rf);
    opti.set_value(Xk1,saida(:,k-1));
    opti.set_value(Xk,saida(:,k));
    opti.set_value(Uk1,entrada(:,k-1));
    opti.set_value(Qk,[CAf0;Tf0;Tcf0]);
    opti.set_value(Qk1,[CAf0;Tf0;Tcf0]);
    
    % Inicialização X e dU com os valores atuais (warm start)
    opti.set_initial(dU,[du_oti_ant(:,1:end-1), zeros(m,1)]);    

    % resolve problema de otimização
    tic
    sol = opti.solve();
    x2nm(k) = toc;

    du_oti_ant = sol.value(dU); 
    du(:,k) = sol.value(dU(:,1));
    entrada(:,k) = entrada(:,k-1) + du(:,k);
    
end

saidaNMPC = saida;
entradaNMPC = entrada;
duNMPC = du;









%% Geração das figuras
cores = gray(5);
cores = cores(1:end-1,:);

t = (0:(nit)-nin)*Ts;
ind = nin:nit;


hf= figure
hf.Position = [680-400 383-300 998 595];
h=subplot(4,1,1)
plot(t,saidaPNMPC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaNEPSAC(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaNMPC(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))

ylim([0.078 0.103])
hl = legend('PNMPC','NEPSAC','NMPC','Referência','Location','SouthEast')

ylabel('C_a (mol/L)','FontSize', tamletra)
xlim([0 24])
set(h, 'FontSize', tamletra);
grid on
h.YTickLabel = trocaponto(h.YTickLabel)

h=subplot(4,1,2)
plot(t,saidaPNMPC(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaNEPSAC(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaNMPC(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))

ylabel('T (K)','FontSize', tamletra)
xlim([0 24])
set(h, 'FontSize', tamletra);
grid on
ylim([438 444])

h=subplot(4,1,3)
plot(t,entradaPNMPC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradaNEPSAC(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradaNMPC(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('F (L/min)','FontSize', tamletra)
xlim([0 24])
set(h, 'FontSize', tamletra);
grid on
ylim([75 135])

h=subplot(4,1,4)
plot(t,entradaPNMPC(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradaNEPSAC(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradaNMPC(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('F_c (L/min)','FontSize', tamletra)
xlabel('Tempo (min)')
xlim([0 24])
set(h, 'FontSize', tamletra);
grid on
ylim([80 125])


hl.Position = [0.7819 0.7136 0.1152 0.1235];

% print('casoCSTR_comparativo_ref','-depsc')

%% 
cores = gray(5);
cores = cores(1:end-1,:);

t = (0:(nit)-nin)*Ts;
ind = nin:nit;


hf= figure
hf.Position = [680-400 383-300 998 595];
h=subplot(4,1,1)
plot(t,saidaPNMPC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaNEPSAC(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaNMPC(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))

ylim([0.098 0.103])
hl = legend('PNMPC','NEPSAC','NMPC','Referência','Location','SouthEast')

ylabel('C_a (mol/L)','FontSize', tamletra)
xlim([23 43])
set(h, 'FontSize', tamletra);
grid on
h.YTickLabel = trocaponto(h.YTickLabel)

h=subplot(4,1,2)
plot(t,saidaPNMPC(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaNEPSAC(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaNMPC(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))


ylim([437.5 440])
ylabel('T (K)','FontSize', tamletra)
xlim([23 43])
set(h, 'FontSize', tamletra);
grid on

h=subplot(4,1,3)
plot(t,entradaPNMPC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradaNEPSAC(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradaNMPC(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('F (L/min)','FontSize', tamletra)
xlim([23 43])
ylim([92.5 102.5])
set(h, 'FontSize', tamletra);
grid on

h=subplot(4,1,4)
plot(t,entradaPNMPC(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradaNEPSAC(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradaNMPC(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('F_c (L/min)','FontSize', tamletra)
xlabel('Tempo (min)')
xlim([23 43])
set(h, 'FontSize', tamletra);
grid on
ylim([103 113])


hl.Position = [0.8019 0.8447 0.1152 0.1235];

% print('casoCSTR_comparativo_pert','-depsc')

%%
fator = 4;

ind = nin:nit;
ind = ind(1:fator:end);
t = (0:(nit)-nin);
t = t(1:fator:end);
hf = figure
semilogy(t,x2(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogy(t,x2n(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogy(t,x2nm(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

xlim([0 nit])
hl = legend('PNMPC','NEPSAC','NMPC','Location','NorthWest')

ylabel('Tempo (segundos)','FontSize', tamletra)
xlabel('Tempos (amostras)')

set(gca, 'FontSize', tamletra);
grid on

% hl.Position = [0.7834 0.7304 0.1333 0.1143];

hf.Position = tamfigura;
% hl.Position = [0.7471 0.5088 0.2054 0.2371];
% print('casoCSTR_comparativo_tempos','-depsc')


%% indices
IAEpnmpc1 = sum(abs(saidaPNMPC(1,ind)-refs(1,ind))*Ts);
IAEnepsac1 = sum(abs(saidaNEPSAC(1,ind)-refs(1,ind))*Ts);
IAEnmpc1 = sum(abs(saidaNMPC(1,ind)-refs(1,ind))*Ts);

IAEpnmpc2 = sum(abs(saidaPNMPC(2,ind)-refs(2,ind))*Ts);
IAEnepsac2 = sum(abs(saidaNEPSAC(2,ind)-refs(2,ind))*Ts);
IAEnmpc2 = sum(abs(saidaNMPC(2,ind)-refs(2,ind))*Ts);

% 'IAEs'
% [IAEpnmpc1 IAEnepsac1 IAEnmpc1]
% [IAEpnmpc2 IAEnepsac2 IAEnmpc2]

'IAEs normalizados pelo NMPC'
[IAEpnmpc1 IAEnepsac1 IAEnmpc1]./IAEnmpc1
[IAEpnmpc2 IAEnepsac2 IAEnmpc2]./IAEnmpc2

%%%%
ind2 = nin+250:nit;
IAEpnmpcq1 = sum(abs(saidaPNMPC(1,ind2)-refs(1,ind2))*Ts);
IAEnepsacq1 = sum(abs(saidaNEPSAC(1,ind2)-refs(1,ind2))*Ts);
IAEnmpcq1 = sum(abs(saidaNMPC(1,ind2)-refs(1,ind2))*Ts);

IAEpnmpcq2 = sum(abs(saidaPNMPC(2,ind2)-refs(2,ind2))*Ts);
IAEnepsacq2 = sum(abs(saidaNEPSAC(2,ind2)-refs(2,ind2))*Ts);
IAEnmpcq2 = sum(abs(saidaNMPC(2,ind2)-refs(2,ind2))*Ts);


% 'IAEs para perts'
% [IAEpnmpcq1 IAEnepsacq1 IAEnmpcq1]
% [IAEpnmpcq2 IAEnepsacq2 IAEnmpcq2]

'IAEs para perts normalizados pelo NMPC'
[IAEpnmpcq1 IAEnepsacq1 IAEnmpcq1]./IAEnmpcq1
[IAEpnmpcq2 IAEnepsacq2 IAEnmpcq2]./IAEnmpcq2

%%%%
ind3 = nin:nin+250;
IAEpnmpcr1 = sum(abs(saidaPNMPC(1,ind3)-refs(1,ind3))*Ts);
IAEnepsacr1 = sum(abs(saidaNEPSAC(1,ind3)-refs(1,ind3))*Ts);
IAEnmpcr1 = sum(abs(saidaNMPC(1,ind3)-refs(1,ind3))*Ts);

IAEpnmpcr2 = sum(abs(saidaPNMPC(2,ind3)-refs(2,ind3))*Ts);
IAEnepsacr2 = sum(abs(saidaNEPSAC(2,ind3)-refs(2,ind3))*Ts);
IAEnmpcr2 = sum(abs(saidaNMPC(2,ind3)-refs(2,ind3))*Ts);


% 'IAEs para perts'
% [IAEpnmpcq1 IAEnepsacq1 IAEnmpcq1]
% [IAEpnmpcq2 IAEnepsacq2 IAEnmpcq2]

'IAEs para refs normalizados pelo NMPC'
[IAEpnmpcr1 IAEnepsacr1 IAEnmpcr1]./IAEnmpcr1
[IAEpnmpcr2 IAEnepsacr2 IAEnmpcr2]./IAEnmpcr2

%%
%%% Modelo CSTR MIMO completo livro Bequette Process Dynamics
function [x] = modeloCSTR(y,u,q)
    global ER V k0 k1 k2 k3 Ts optsode
    Ca = y(1);
    T = y(2);

    f = u(1);
    fc = u(2);

    caf = q(1);
    tf = q(2);
    tcf = q(3);
    %%% definição do modelo com as entradas e perturbações constantes. Isto
    %%% simula o sustentador de ordem zero.
    %%% y = [Ca;T];
    fcstr = @(t,y) [f/V*(caf-y(1))-k0*y(1)*exp(-ER/y(2)); % derivada de Ca
                    f/V*(tf-y(2)) + k1*y(1)*exp(-ER/y(2)) + k2*fc*(1-exp(-k3/fc))*(tcf-y(2))]; % derivada de T

    [t,x] = ode45(fcstr,[0 Ts],[Ca;T],optsode);
    x = x(end,:)';
                    
end

%% Modelo CSTR nominal utilizado pelos controladores
function [yo] = modeloCSTRnominal(y,u,q)
    global ER V k0 k1 k2 k3 Ts 
    Ca = y(1);
    T = y(2);

    f = u(1);
    fc = u(2);

    caf = q(1);
    tf = q(2);
    tcf = q(3);

    yo(1,1) = Ca+ Ts*(f/V*(caf-Ca) - k0*Ca*exp(-ER/T));
    yo(2,1) = T+ Ts*(f/V*(tf-T) + k1*Ca*exp(-ER/T) + k2*fc*(1-exp(-k3/fc))*(tcf-T));

    %%% saturação dos estados
    if(yo(1,1)<0)
        yo(1,1)= 0;
    end

    if(yo(2,1)<0)
        yo(2,1) = 0;
    end


end
