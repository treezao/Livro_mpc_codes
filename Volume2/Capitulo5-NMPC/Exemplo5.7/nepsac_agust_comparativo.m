%%% Exemplo do tanque não linear comparando GPC linear, GPC com 
%%% resposta livre não linear e o NEPSAC

clear all
close all
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

%% parametros comuns aos controladores 
N1 = 1;
N2 = 3;
N = N2-N1+1;
Nu = 2;
Nq = 0;

lambda = 10; % ponderação do incremento de controle na função custo.
delta = 1; % ponderação do erro na função custo

Qe = delta*eye(N2-N1+1);
Qu = lambda*eye(Nu);

C = 1; % polinômio C do modelo da perturbação
D = [1 -1]; % polinômio D do modelo da perturbação

nc = size(C,2);
nd = size(D,2);

%% parametros gerais da simulação
%%% condições iniciais
y0 = 0;
u0 = -0.25;
q0 = 0;

y0 = -0.67;
u0 = (y0-(4.*y0.*y0)./(1+y0.^2+y0.^2)-cos(0.5*(y0+y0)))./1.2;

dn = 0; % atraso nominal

%%% parâmetros do modelo

nin = 10;
nit = nin+200; % numero de iterações

%%% definição da referência
ref = y0*ones(1,nit);
ref(nin+20:end) = -2;
ref(nin+100:end) = 0.5;

%%% definição da perturbação
q = q0*ones(1,nit);
q(nin+50:end) = 1;
q(nin+150:end) = -1;


%% Configuração do GPC linear
%%% Modelo da planta linearizado no ponto de operação
a1 = (4.0*y0*(1+y0^2+y0^2)-4.0*y0*y0*2*y0)/(1+y0^2+y0^2)^2-0.5*sin(0.5*(y0+y0));
a2 = a1;
b1 = 1.2;

A = [1 -a1 -a2];
B = b1; 

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A


%%% Obtenção das matrizes do GPC linear
Btil = B; % incorporação do atraso no numerador B

G = zeros(N,N); % matriz dinâmica G
H = zeros(N,nb);

[E,F] = diofantina(conv(A,[1 -1]),N1,N2);

for i=N1:N2
    EjB = conv(E(i-N1+1,1:i),Btil);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

    
end
G = G(:,1:Nu);

Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;


%% simulação do GPC linear
saida= y0*ones(1,nit);
entrada = u0*ones(1,nit);
du = zeros(1,nit);

rfant = y0;


for k=nin:nit
    %%% simulacao do modelo
    saida(k) = modeloSim(saida(:,k-1:-1:k-2),entrada(:,k-1),q(:,k-1));
    
    %% controlador GPC
    % calculo das referencias
    
    R = ones(N,1)*ref(k);
    
    %%% cálculo da resposta livre;
    f = F*saida(1,k:-1:k-na)';
    
    if(~isempty(H))
        f = f+ H*du(1,k-1:-1:k-nb)'; % parcela dos incrementos de controle
    end
    
   
    %%% cálculo do incremento de controle ótimo
    duOti = quadprog(Hqp,fqp1*(R-f),[],[]);
    du(k) = duOti(1);
    
    %%% cálculo do sinal de controle ótimo
    entrada(k) = du(k)+entrada(k-1);
 

end

saidaGPC = saida;
entradaGPC = entrada;
duGPC = du;




%% simulação do GPC com resposta livre não linear
saida= y0*ones(1,nit);
entrada = u0*ones(1,nit);
du = zeros(1,nit);

saidaModelo = y0*ones(1,nit);

eta = zeros(1,nit+N2);
e = zeros(1,nit+N2);


for k=nin:nit
    %%% simulacao do CSTR
    saida(k) = modeloSim(saida(:,k-1:-1:k-2),entrada(:,k-1),q(:,k-1));

    %% controlador GPC com modelo não linear na resposta livre
    % calculo das referencias
    
    R = ones(N,1)*ref(k);

    
    %%% cálculo da resposta livre;
    saidaModelo(k) = modeloSim(saida(k-1:-1:k-2),entrada(k-1),q0);
    
    %%% calculo dos etas futuros
    eta(k) = saida(k) - saidaModelo(k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:N2
        eta(k+j) = -D(2:end)*eta(k+j-1:-1:k+j-nd+1)' + C*e(k+j:-1:k+j-nc+1)';
    end

    yb = zeros(N2,1); % resposta base
                        % yb=[y(k+1);y(k+2);...;y(k+N2)];
    ub = [entrada(k-dn:k-1), entrada(k-1)*ones(1,N2)]; % u base inicial

    yb(1) = modeloSim(saida(k:-1:k-1),ub(1),q0) + eta(k+1); % yb(k+1|k)
    yb(2) = modeloSim([yb(1),saida(k)],ub(2),q0) + eta(k+2); % yb(k+2|k)
    for i=3:N2
        yb(i) = modeloSim(yb(i-1:-1:i-2),ub(i),q0) + eta(k+i); %yb(k+i|k)
    end
    
    f = yb(N1:N2);
    
    %%% cálculo do incremento de controle ótimo
    duOti = quadprog(Hqp,fqp1*(R-f),[],[]);
    du(k) = duOti(1);
    
    %%% cálculo do sinal de controle ótimo
    entrada(k) = du(k)+entrada(k-1);
 

end
saidaGPCnlin = saida;
entradaGPCnlin = entrada;
duGPCnlin = du;



%% Configuração do NEPSAC
T = tril(ones(Nu));
T1 = inv(T);

%% simulação do NEPSAC
saida= y0*ones(1,nit);
entrada = u0*ones(1,nit);
du = zeros(1,nit);

saidaModelo = y0*ones(1,nit);

eta = zeros(1,nit+N2);
e = zeros(1,nit+N2);


for k=nin:nit
    %%% simulacao do CSTR
    saida(k) = modeloSim(saida(:,k-1:-1:k-2),entrada(:,k-1),q(:,k-1));
    
    %% controlador NEPSAC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(k) = modeloSim(saida(k-1:-1:k-2),entrada(k-1),q0);
    
    %%% calculo dos etas futuros
    eta(k) = saida(k) - saidaModelo(k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:N2
        eta(k+j) = -D(2:end)*eta(k+j-1:-1:k+j-nd+1)' + C*e(k+j:-1:k+j-nc+1)';
    end
    
    ub = [entrada(k-dn:k-1), entrada(k-1)*ones(1,N2)]; % u base inicial
                          %ub = [u(k-d) u(k-d+1) ... u(k) u(k+1)... u(k+N2-1)
    dunorm = inf; % inicilizar norma de delta u como infinito
   
    % calculo das referencias
    R = ones(N,1)*ref(k);

    
    cont = 0;
    while (dunorm > 1e-7 && cont <= 100) % iterar até obter a resposta dentro do erro especificado
        %%% Obtenção da resposta base
        yb = zeros(N2,1); % resposta base
                            % yb=[y(k+1);y(k+2);...;y(k+N2)];

        yb(1) = modeloSim(saida(k:-1:k-1),ub(1),q0) + eta(k+1); % yb(k+1|k)
        yb(2) = modeloSim([yb(1),saida(k)],ub(2),q0) + eta(k+2); % yb(k+2|k)
        for i=3:N2
            yb(i) = modeloSim(yb(i-1:-1:i-2),ub(i),q0) + eta(k+i); %yb(k+i|k)
        end

        
        %%% cálculo do controle ótimo

        %%% cálculo de Ge
        Ge = zeros(N2,Nu);
        u00 = entrada(k-1)/100; 
%         u00 = 1;
        for j=1:Nu
            yot = zeros(N2,1);
            udelta = ub;
            if(j<Nu)
                udelta(j+dn) = udelta(j+dn)+u00;
            else
                udelta(j+dn:end) = udelta(j+dn:end)+u00;
            end

            yot(1) = modeloSim(saida(k:-1:k-1),udelta(1),q0) + eta(k+1);
            yot(2) = modeloSim([yot(1),saida(k)],udelta(2),q0) + eta(k+2);
            for i=3:N2
                yot(i) = modeloSim(yot(i-1:-1:i-2),udelta(i),q0) + eta(k+i);
            end
                        
            Ge(:,j) = (yot - yb)/(u00);
        end
        Ge = Ge(N1:N2,:);
        
        
        
        Hqp = 2*(Ge'*Qe*Ge+T1'*Qu*T1);
        fqp = 2*(Ge'*Qe*(yb(N1:end)-R)+T1'*Qu*T1*(-entrada(k-1)*ones(Nu,1)+ub(1:Nu)'));
        
        X = quadprog(Hqp,fqp,[],[]);
        
        ub(dn+1:dn+Nu) = ub(dn+1:dn+Nu)+X';
        ub(dn+1+Nu:end) = ub(dn+Nu);
        
        dunorm = norm(X);
        cont = cont+1;
    end
%     temp(1,k) = cont;
%     temp(2,k) = dunorm;
    
    entrada(k) = ub(dn+1);
 

end



%%
%% Geração das figuras
cores = gray(5);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin);

hf = figure
h=subplot(2,1,1)
plot(t,saida(1,nin:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaGPC(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaGPCnlin(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,ref(1,nin:nit),'-.','LineWidth',tamlinha,'Color',cores(4,:))

hl = legend('NEPSAC','GPC linear','GPC não linear','Referência','Location','SouthEast')

ylabel('Saída','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on
xlim([0 t(end)])
ylim([-3.5 1.3])

h = subplot(2,1,2)
plot(t,entrada(1,nin:nit),'LineWidth',tamlinha,'Color',cores(3,:))
hold on
plot(t,entradaGPC(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(1,:))
plot(t,entradaGPCnlin(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Entrada','FontSize', tamletra)
xlabel('Tempo (amostras)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);
xlim([0 t(end)])



hf.Position = tamfigura;
hl.Position = [0.7417 0.5218 0.2536 0.2371];

% print('nepsac_agust_comp_completo','-depsc')


%% figura 1
hf = figure
h=subplot(2,1,1)
plot(t,saida(1,nin:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaGPC(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaGPCnlin(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,ref(1,nin:nit),'-.','LineWidth',tamlinha,'Color',cores(4,:))

hl = legend('NEPSAC','GPC linear','GPC não linear','Referência','Location','SouthEast')

ylabel('Saída','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on
xlim([0 99])
ylim([-2.5 -0.2])

h = subplot(2,1,2)
plot(t,entrada(1,nin:nit),'LineWidth',tamlinha,'Color',cores(3,:))
hold on
plot(t,entradaGPC(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(1,:))
plot(t,entradaGPCnlin(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Entrada','FontSize', tamletra)
xlabel('Tempo (amostras)','FontSize', tamletra)

grid on
set(h, 'FontSize', tamletra);
xlim([0 99])
ylim([-4.1 -1.9])



hf.Position = tamfigura;
hl.Position = [0.6774 0.7218 0.2536 0.2371];

% print('nepsac_agust_comp_parte1','-depsc')

%% figure 2
hf = figure
h=subplot(2,1,1)
plot(t,saida(1,nin:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaGPC(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaGPCnlin(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,ref(1,nin:nit),'-.','LineWidth',tamlinha,'Color',cores(4,:))

hl = legend('NEPSAC','GPC linear','GPC não linear','Referência','Location','SouthEast')

ylabel('Saída','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on
xlim([90 t(end)])
ylim([-3.5 1.3])

h = subplot(2,1,2)
plot(t,entrada(1,nin:nit),'LineWidth',tamlinha,'Color',cores(3,:))
hold on
plot(t,entradaGPC(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(1,:))
plot(t,entradaGPCnlin(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(2,:))

ylabel('Entrada','FontSize', tamletra)
xlabel('Tempo (amostras)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);
xlim([90 t(end)])
ylim([-4.1 1])


hf.Position = tamfigura;
hl.Position = [0.6881 0.5089 0.2536 0.2371];

% print('nepsac_agust_comp_parte2','-depsc')

%% Modelo do sistema
function [x] = modeloSim(y,u,q)
    yk1 = y(1); %y(k-1)
    yk2 = y(2); %y(k-2)
    
    uq = u+q;
    x = (4.0*yk1*yk2)/(1+yk1^2+yk2^2)+1.2*uq+cos(0.5*(yk1+yk2));
    
end
