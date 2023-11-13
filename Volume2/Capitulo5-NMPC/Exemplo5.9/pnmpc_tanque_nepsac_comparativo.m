%%% Exemplo do tanque não linear comparando GPC linear, GPC com 
%%% resposta livre não linear e o NEPSAC

clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% parametros comuns aos controladores 
N1 = 1;
N2 = 30;
N = N2-N1+1;
Nu = 5;
Nq = 0;

lambda = 1; % ponderação do incremento de controle na função custo.
delta = 1; % ponderação do erro na função custo

Qe = delta*eye(N2-N1+1);
Qu = lambda*eye(Nu);

C = 1; % polinômio C do modelo da perturbação
D = [1 -1]; % polinômio D do modelo da perturbação

nc = size(C,2);
nd = size(D,2);

umax = 100; % valor máximo da abertura da válvula
umin = 0; % valor mínimo da abertura da válvula

%% parametros gerais da simulação
nit = 52; % numero de iterações
nin = 2;

%%% condições iniciais
h0 = 50;
u0 = 31.3209;
q0 = 20;

%%% parâmetros do modelo
global K1;
global K2;
global g;
K1 = 2;
K2 = 0.1;
g = 9.81;


%%% definição da referência
ref = h0*ones(1,nit);
ref(5:end) = 90;
ref(20:end) = 30;
ref(40:end) = 80;

%%% definição da perturbação
q = q0*ones(1,nit);
q(30:end) = 30;


%% Configuração do NEPSAC
T = tril(ones(Nu));
T1 = inv(T);

%% simulação do NEPSAC
saida= h0*ones(1,nit);
saidaModelo = h0*ones(1,nit);
entrada = u0*ones(1,nit);
du = zeros(1,nit);

eta = zeros(1,nit+N2);
e = zeros(1,nit+N2);



for k=nin:nit
    %%% simulacao do tanque
    saida(k) = modeloTanque(saida(k-1),entrada(k-1),q(k-1));
    
    tic
    %% controlador NEPSAC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(k) = modeloTanque(saida(k-1),entrada(k-1),q(k-1));
    
    %%% calculo dos etas futuros
    eta(k) = saida(k) - saidaModelo(k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:N2
        eta(k+j) = -D(2:end)*eta(k+j-1:-1:k+j-nd+1)' + C*e(k+j:-1:k+j-nc+1)';
    end
    
    
    ub = entrada(k-1)*ones(N2,1); % u base inicial
    dunorm = inf; % inicilizar norma de delta u como infinito
    
    rfut = ref(k)*ones(N2-N1+1,1);
    
    cont = 0;
    while (dunorm > 1e-3 && cont <= 20) % iterar até obter a resposta dentro do erro especificado
        %%% Obtenção da resposta base
        yb = zeros(N2,1); % resposta base

        yb(1) = modeloTanque(saida(k),ub(1),q(k-1)) + eta(k+1); % yb(k+1|k)
        for i=2:N2
            yb(i) = modeloTanque(yb(i-1),ub(i),q(k-1)) + eta(k+i);
        end


        %%% cálculo do controle ótimo
        udelta = ub;
        
        %%% cálculo de Ge
        Ge = zeros(N2,Nu);
        for j=1:Nu
            yot = zeros(N2,1);
            udelta = ub;
            if(j<Nu)
                udelta(j) = udelta(j)+1;
            else
                udelta(j:end) = udelta(j:end)+1;
            end
            
            yot(1) = modeloTanque(saida(k),udelta(1),q(k-1)) + eta(k+1);
            for i=2:N2
                yot(i) = modeloTanque(yot(i-1),udelta(i),q(k-1)) + eta(k+i);
            end
                        
            Ge(:,j) = yot - yb;
        end
        Ge = Ge(N1:N2,:);
        
        
        
        Hqp = 2*(Ge'*Qe*Ge+T1'*Qu*T1);
        fqp = 2*(Ge'*Qe*(yb-rfut)+T1'*Qu*T1*(-entrada(k-1)*ones(Nu,1)+ub(1:Nu)));
        
        X = quadprog(Hqp,fqp,[],[]);
        
        ub(1:Nu) = ub(1:Nu)+X;
        ub(Nu+1:end) = ub(Nu);
        
        dunorm = norm(X);
        cont = cont+1;
    end
    temp(k) = cont;
    
    entrada(k) = ub(1);
    x(k) = toc;

end

saidaNEPSAC = saida;
entradaNEPSAC = entrada;
duNEPSAC = du;
xNEPSAC = x;



%% Simulação do PNMPC
saida= h0*ones(1,nit);
saidaModelo = h0*ones(1,nit);
entrada = u0*ones(1,nit);
du = zeros(1,nit);

eta = zeros(1,nit+N2);
e = zeros(1,nit+N2);



for k=nin:nit
    %%% simulacao do tanque
    saida(k) = modeloTanque(saida(k-1),entrada(k-1),q(k-1));
    
    tic
    %% controlador PNMPC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(k) = modeloTanque(saida(k-1),entrada(k-1),q(k-1));
    
    %%% calculo dos etas futuros
    eta(k) = saida(k) - saidaModelo(k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:N2
        eta(k+j) = -D(2:end)*eta(k+j-1:-1:k+j-nd+1)' + C*e(k+j:-1:k+j-nc+1)';
    end
    
    %%% referências futuras
    rfut = ref(k)*ones(N2-N1+1,1);

    %%% calculo de G
    y0 = zeros(N2,1);
    y0(1) = modeloTanque(saida(k),entrada(k-1),q(k-1));
    for i=2:N2
        y0(i) = modeloTanque(y0(i-1),entrada(k-1),q(k-1));
    end

    %%% calculo de yi para obter G
    if(u0 == 0)
        epsilon = 0.001;
    else
        epsilon = u0/1000;
    end
    
    y = zeros(N2,Nu);
    for j=1:Nu
        u = entrada(k-1)*ones(N2,1);
        u(j:end) = u(j:end)+epsilon;

        y(1,j) = modeloTanque(saida(k),u(1),q(k-1));
        for i=2:N2
            y(i,j) = modeloTanque(y(i-1,j),u(i),q(k-1));
        end
    end


    G = (y-repmat(y0,1,Nu))./epsilon;
    
    %%% cálculo da resposta livre corrigida
    yc = zeros(N2,1);
    yc(1) = modeloTanque(saida(k),entrada(k-1),q(k-1)) + eta(k+1);
    for i=2:N2
        yc(i) = modeloTanque(yc(i-1),entrada(k-1),q(k-1)) + eta(k+i);
    end
    
    f = yc(N1:N2);
    
    Hqp = 2*(G'*Qe*G+Qu);
    fqp = -2*G'*Qe*(rfut-f);

    X = quadprog(Hqp,fqp,[],[]);
        
    du(k) = X(1);
    entrada(k) = entrada(k-1)+du(k); 
    x(k) = toc;
end

saidaPNMPC = saida;
entradaPNMPC = entrada;
duPNMPC = du;
xPNMPC = x;

%% Simulação do PNMPC
saida= h0*ones(1,nit);
saidaModelo = h0*ones(1,nit);
entrada = u0*ones(1,nit);
du = zeros(1,nit);

eta = zeros(1,nit+N2);
e = zeros(1,nit+N2);

No = 5; %% calcular a matriz G a cada 5 amostras.

for k=nin:nit
    %%% simulacao do tanque
    saida(k) = modeloTanque(saida(k-1),entrada(k-1),q(k-1));
    
    tic
    %% controlador PNMPC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(k) = modeloTanque(saida(k-1),entrada(k-1),q(k-1));
    
    %%% calculo dos etas futuros
    eta(k) = saida(k) - saidaModelo(k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:N2
        eta(k+j) = -D(2:end)*eta(k+j-1:-1:k+j-nd+1)' + C*e(k+j:-1:k+j-nc+1)';
    end
    
    %%% referências futuras
    rfut = ref(k)*ones(N2-N1+1,1);

    if(k == 2 || mod(k,No)==0)
        %%% calculo de G
        y0 = zeros(N2,1);
        y0(1) = modeloTanque(saida(k),entrada(k-1),q(k-1));
        for i=2:N2
            y0(i) = modeloTanque(y0(i-1),entrada(k-1),q(k-1));
        end

        %%% calculo de yi para obter G
        if(u0 == 0)
            epsilon = 0.001;
        else
            epsilon = u0/1000;
        end

        y = zeros(N2,Nu);
        for j=1:Nu
            u = entrada(k-1)*ones(N2,1);
            u(j:end) = u(j:end)+epsilon;

            y(1,j) = modeloTanque(saida(k),u(1),q(k-1));
            for i=2:N2
                y(i,j) = modeloTanque(y(i-1,j),u(i),q(k-1));
            end
        end


        G = (y-repmat(y0,1,Nu))./epsilon;
    end
    
    %%% cálculo da resposta livre corrigida
    yc = zeros(N2,1);
    yc(1) = modeloTanque(saida(k),entrada(k-1),q(k-1)) + eta(k+1);
    for i=2:N2
        yc(i) = modeloTanque(yc(i-1),entrada(k-1),q(k-1)) + eta(k+i);
    end
    
    f = yc(N1:N2);
    
    Hqp = 2*(G'*Qe*G+Qu);
    fqp = -2*G'*Qe*(rfut-f);

    X = quadprog(Hqp,fqp,[],[]);
        
    du(k) = X(1);
    entrada(k) = entrada(k-1)+du(k); 
    x(k) = toc;
end

saidaPNMPC2 = saida;
entradaPNMPC2 = entrada;
duPNMPC2 = du;
xPNMPC2 = x;

%%
%% Geração das figuras
cores = gray(5);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin);

hf = figure
h=subplot(2,1,1)
plot(t,saidaPNMPC(1,nin:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaNEPSAC(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaPNMPC2(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,ref(1,nin:nit),'-.','LineWidth',tamlinha,'Color',cores(4,:))
% ylim([0 100])
% xlim([0 nit-nin])
hl = legend('PNMPC','NEPSAC','PNMPC2','Referência','Location','SouthEast')
% hl.Position = [0.0756 0.5085 0.2375 0.2131];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on


h = subplot(2,1,2)
plot(t,entradaPNMPC(1,nin:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradaNEPSAC(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradaPNMPC2(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('Manipulada','FontSize', tamletra)
xlabel('Tempo (amostras)','FontSize', tamletra)
ylim([17 62])
% xlim([0 nit-nin])
grid on
set(h, 'FontSize', tamletra);
% ylim([-5 5])

% h=subplot(3,1,3)
% plot(t,erro(1:nit-N2-1),'LineWidth',tamlinha,'Color',cores(1,:))
% title('Erro', 'FontSize', tamletra);
% xlabel('tempo', 'FontSize', tamletra);
% set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7471 0.5088 0.2054 0.2371];
% print('pnmpc_tanque_nepsac_comparativo','-depsc')



%% função do modelo do sistema de tanque
function hk = modeloTanque(hk1,ak1,qk1)
    global K1;
    global K2;
    global g;

    %%% modelo discreto do tanque
    hk = hk1 + K1*ak1 -K2*qk1*sqrt(2*g*hk1);
    
    if(hk<0)
        hk = 0;
    end
end

