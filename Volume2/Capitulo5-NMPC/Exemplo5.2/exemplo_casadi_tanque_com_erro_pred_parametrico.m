%%% código que utiliza incrementos de controle, mas não tem estimação da
%%% perturbação não medida... daí não rejeita perts

clear all
close all
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%%
import casadi.*

% Define parâmetros do problema

K1=2;             
K2=.1;  
g=9.8;                    % Aceleração da gravidade  [m/s^2]
tsim=100;                 % Tempo de simulação  [s]
ts=1;                     % Período de amostragem [s]
nit=round(tsim/ts);       % número de iterações

ns=1;                 % Numero de estados
m=1;                  % Numero de entradas
n=1;                  % Numero de saídas
mq = 1;               % número de perturbações
N=20;                 % Horizonte de predição

Qx = eye(n);        % Ponderação estado
Qu = 10;             % Ponderação controle

qs0 = 50;               % valor inicial da perturbação
u0 = 31.9438;           % valor inicial do sinal de controle
xs = 50*ones(n,nit);    % Estados
us = u0*ones(m,nit);    % Controle
qs = qs0*ones(n,nit);   % Perturbação
xr = 50*ones(n,nit);    % Referência
du = zeros(m,nit);      % incrementos de controle
J = zeros(1,nit);       % Objetivo

Umax = 100;             % Valor máximo da entrada
Umin = 0;               % Valor mínimo da entrada

Xmax = 100;         % Valor máximo do estado
Xmin = 0;           % Valor mínimo do estado


% valor da referência ao longo da simulação
xr(1,round(tsim*0.1/ts):end) = 45;

% Valor da perturbação ao longo da simulação
qs(1,round(nit/2):end) = 45;

% definição da função do modelo para facilitar a notação
f_modelo = @(x,u,q) x + K1*q - K2*u*sqrt(2*g*max(x,1e-3));


%% Montagem do problema de otimização de forma paramétrica (acelera resolução)
opti = casadi.Opti();   % Chamada da interface "Opti" de Casadi

f_custo = 0;    % Função objetivo

% Variáveis de otimização do NLP
X = opti.variable(n, N);    % Variável dos estados de k=1 até k=N
dU = opti.variable(m, N);   % Variável de control de k=0 até k=Nu

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
for j=1:N
    f_custo = f_custo + (Xr-X(:,j))'*Qx*(Xr-X(:,j)) ...
                 + dU(:,j)'*Qu*dU(:,j);
    Uant = U; % sinal de controle anterior
    U = U + dU(:,j);
    if(j==1)
        opti.subject_to(X(:,j) == f_modelo(Xk,U,Qk) + Xk -f_modelo(Xk1,Uant,Qk1) );
    elseif(j==2)
        opti.subject_to(X(:,j) == f_modelo(X(:,j-1),U,Qk) + X(:,j-1) - f_modelo(Xk,Uant,Qk) );
    else
        opti.subject_to(X(:,j) == f_modelo(X(:,j-1),U,Qk) + X(:,j-1) - f_modelo(X(:,j-2),Uant,Qk));
    end
    % Restrições de operação
    opti.subject_to(Xmin <= X(:,j) <= Xmax);
    opti.subject_to(Umin <= U <= Umax);
end

opti.minimize(f_custo);     % Declarando o objetivo: minimizar f_custo

%%% configurações do solver
solver = 'ipopt';   % Outro solver interessante: qrsqp

solver_opts = struct('print_level', 5, 'warm_start_init_point', 'yes');
casadi_opts = struct('print_time', false);
opti.solver(solver, casadi_opts, solver_opts);  % Construção do solver



%% simulação
% Laço de simulação em i
for k=2:nit
    %% simulação do modelo
    xs(k)= xs(k-1) + K1*qs(k-1) - K2*us(k-1)*sqrt(2*g*xs(k-1));   
    
    %%% Controlador NMPC com CasADi
    %% Soluciona o problema de otimização
    
    %%% definição dos parâmetros 
    opti.set_value(Xr,xr(k));
    opti.set_value(Xk1,xs(k-1));
    opti.set_value(Xk,xs(k));
    opti.set_value(Uk1,us(k-1));
    opti.set_value(Qk,qs0);
    opti.set_value(Qk1,qs0);
    
    % Inicialização X e dU com os valores atuais (warm start)
    opti.set_initial(dU,0);
    opti.set_initial(X,xs(k));

    sol = opti.solve();
    % Vetor com os valores do objetivo
    J(k) = sol.value(f_custo);
    % Vetor com os valores da entrada
    du(:,k) = sol.value(dU(:,1));
    us(:,k) = us(:,k-1) + du(:,k);
end

%% 
cores = gray(4);
cores = cores(1:end-1,:);

hf = figure
h=subplot(2,1,1)
plot(xs(1,2:nit)','LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(xr(:,2:nit)','-.','LineWidth',tamlinha,'Color',cores(2,:))
% ylim([0 2])
ylabel('Controladas (%)','FontSize', tamletra)
hl0 = legend('Saídas','Referência','Location','NorthEast')
% ylim([43 52])
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(us(1,2:nit)','LineWidth',tamlinha,'Color',cores(1,:))

ylabel('Manipuladas (%)','FontSize', tamletra)

grid on
set(h, 'FontSize', tamletra);
% ylim([28 35])


xlabel('Tempo (segundos)','FontSize', tamletra)
grid on
% ylim([-5 5])
set(h, 'FontSize', tamletra);

% h=subplot(3,1,3)
% plot(t,erro(1:nit-N2-1),'LineWidth',tamlinha,'Color',cores(1,:))
% title('Erro', 'FontSize', tamletra);
% xlabel('tempo', 'FontSize', tamletra);
% set(h, 'FontSize', tamletra);

hf.Position = tamfigura;

hl0.Position = [0.7310 0.6411 0.2054 0.1242]


hf.Position = tamfigura;
% print('exemplo_casadi_tanque_com_erro_pred','-depsc')
        