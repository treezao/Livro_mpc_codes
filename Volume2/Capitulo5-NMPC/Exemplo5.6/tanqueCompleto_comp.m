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

lambda = 5; % ponderação do incremento de controle na função custo.
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
nin = 2;
nit = nin+60; % numero de iterações

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
ref(nin+5:end) = 90;
ref(nin+30:end) = 50;
ref(nin+45:end) = 80;

%%% definição da perturbação
q = q0*ones(1,nit);
q(15:end) = 30;



%% Configuração do GPC linear
%%% Modelo da planta linearizado no ponto de operação h0,u0,q0
A = [1 (-1+K2*q0*g/sqrt(2*g*h0))];
B = [K1]; 
Bq = [0 -K2*sqrt(2*g*h0)];

d=0;
dq = d;

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A
nbq = size(Bq,2)-1; % ordem do polinômio Bq


%%% Obtenção das matrizes do GPC linear
Btil = conv(B,[zeros(1,d) 1]); % incorporação do atraso no numerador B
Bqtil = conv(Bq,[zeros(1,dq) 1]); % incorporação do atraso no numerador Bq


G = zeros(N,N); % matriz dinâmica G
H = zeros(N,nb+d);
Gq = zeros(N,N-1); 
Hq = zeros(N,nbq+dq);

[E,F] = diofantina(conv(A,[1 -1]),N1,N2);

for i=N1:N2
    EjB = conv(E(i-N1+1,1:i),Btil);
    EjBq = conv(E(i-N1+1,1:i),Bq);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

    Gq(i-N1+1,1:i) = EjBq(i:-1:1);    
    Hq(i-N1+1,:) = EjBq(i+1:end);
    
end
G = G(:,1:Nu);
Gq = Gq(:,1:Nq);

Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;


%% simulação do GPC linear
saidaGPC= h0*ones(1,nit);
entradaGPC = u0*ones(1,nit);
duGPC = zeros(1,nit);


for k=2:nit
    %%% simulacao do tanque
    saidaGPC(k) = modeloTanque(saidaGPC(k-1),entradaGPC(k-1),q(k-1));
    
    %% controlador GPC
    R = ref(k)*ones(N2-N1+1,1);

    %%% cálculo da resposta livre;
    f = F*saidaGPC(1,k:-1:k-na)';
    
    if(~isempty(H))
        f = f+ H*duGPC(1,k-1:-1:k-nb-d)'; % parcela dos incrementos de controle
    end
    
    %%% cálculo do incremento de controle ótimo
    duOti = quadprog(Hqp,fqp1*(R-f),[],[]);
    duGPC(k) = duOti(1);
    
    %%% cálculo do sinal de controle ótimo
    entradaGPC(k) = duGPC(k)+entradaGPC(k-1);
 

end


%% simulação do GPC com resposta livre não linear
saidaGPCnlin= h0*ones(1,nit);
saidaModeloGPCnlin = h0*ones(1,nit);
entradaGPCnlin = u0*ones(1,nit);
duGPCnlin = zeros(1,nit);

etaGPCnlin = zeros(1,nit+N2);
eGPCnlin = zeros(1,nit+N2);


for k=2:nit
    %%% simulacao do tanque
    saidaGPCnlin(k) = modeloTanque(saidaGPCnlin(k-1),entradaGPCnlin(k-1),q(k-1));

    %% controlador GPC com modelo não linear na resposta livre
    R = ref(k)*ones(N2-N1+1,1);
    %%% cálculo da resposta livre;
    saidaModeloGPCnlin(k) = modeloTanque(saidaGPCnlin(k-1),entradaGPCnlin(k-1),q0);
    
    etaGPCnlin(k) = saidaGPCnlin(k) - saidaModeloGPCnlin(k); % cálculo do eta atual
    
    eGPCnlin(k) = -C(2:end)*eGPCnlin(k-1:-1:k-nc+1)'+ D*etaGPCnlin(k:-1:k-nd+1)';
    
    for j=1:N2
        etaGPCnlin(k+j) = -D(2:end)*etaGPCnlin(k+j-1:-1:k+j-nd+1)' + C*eGPCnlin(k+j:-1:k+j-nc+1)';
    end
    
    f = zeros(N2,1);
    f(1) = modeloTanque(saidaGPCnlin(k),entradaGPCnlin(k-1),q0) + etaGPCnlin(k);
    for i=2:N2
        f(i) = modeloTanque(f(i-1),entradaGPCnlin(k-1),q0) + etaGPCnlin(k+i);
    end
    
    f = f(N1:N2);
    
    %%% cálculo do incremento de controle ótimo
    duOti = quadprog(Hqp,fqp1*(R-f),[],[]);
    duGPCnlin(k) = duOti(1);
    
    %%% cálculo do sinal de controle ótimo
    entradaGPCnlin(k) = duGPCnlin(k)+entradaGPCnlin(k-1);
 

end



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



for k=2:nit
    %%% simulacao do tanque
    saida(k) = modeloTanque(saida(k-1),entrada(k-1),q(k-1));
    
    %% controlador NEPSAC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(k) = modeloTanque(saida(k-1),entrada(k-1),q0);
    
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

        yb(1) = modeloTanque(saida(k),ub(1),q0) + eta(k+1); % yb(k+1|k)
        for i=2:N2
            yb(i) = modeloTanque(yb(i-1),ub(i),q0) + eta(k+i);
        end


        %%% cálculo do controle ótimo
        udelta = ub;
        
        %%% cálculo de Ge
        Ge = zeros(N2,Nu);
        for j=1:Nu
            yot = zeros(N2,1);
            udelta = ub;
            epsilon = entrada(k-1)/100;
            if(j<Nu)
                udelta(j) = udelta(j)+epsilon;
            else
                udelta(j:end) = udelta(j:end)+epsilon;
            end
            
            yot(1) = modeloTanque(saida(k),udelta(1),q0) + eta(k+1);
            for i=2:N2
                yot(i) = modeloTanque(yot(i-1),udelta(i),q0) + eta(k+i);
            end
                        
            Ge(:,j) = (yot - yb)/epsilon;
        end
        Ge = Ge(N1:N2,:);
        
        
        
        Hqp = 2*(Ge'*Qe*Ge+T1'*Qu*T1);
        fqp1 = 2*(Ge'*Qe*(yb-rfut)+T1'*Qu*T1*(-entrada(k-1)*ones(Nu,1)+ub(1:Nu)));
        
        X = quadprog(Hqp,fqp1,[],[]);
        
        ub(1:Nu) = ub(1:Nu)+X;
        ub(Nu+1:end) = ub(Nu);
        
        dunorm = norm(X);
        cont = cont+1;
    end
    temp(k) = cont;
    
    entrada(k) = ub(1);
 

end

%%
%% Geração das figuras
cores = gray(5);
cores = cores(1:end-1,:);

t = (nin:nit)-nin;

hf = figure
h=subplot(2,1,1)
plot(t,saidaGPC(1,nin:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaGPCnlin(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saida(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,ref(1,nin:nit),'-.','LineWidth',tamlinha,'Color',cores(4,:))
ylim([45 100])
xlim([0 (nit-nin)])
hl = legend('GPC linear','GPC não linear','NEPSAC','Referência','Location','SouthEast')
ylabel('Controlada (%)','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t,entradaGPC(1,nin:nit),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradaGPCnlin(1,nin:nit),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entrada(1,nin:nit),':','LineWidth',tamlinha,'Color',cores(3,:))

ylabel('Manipulada (%)','FontSize', tamletra)
ylim([30 65])
grid on
set(h, 'FontSize', tamletra);
% ylim([-5 5])
xlim([0 (nit-nin)])
xlabel('Tempo (amostras)','FontSize', tamletra)

% h=subplot(3,1,3)
% plot(t,erro(1:nit-N2-1),'LineWidth',tamlinha,'Color',cores(1,:))
% title('Erro', 'FontSize', tamletra);
% xlabel('tempo', 'FontSize', tamletra);
% set(h, 'FontSize', tamletra);



hf.Position = tamfigura;
hl.Position = [0.7399 0.5057 0.2411 0.1806];
% print('tanqueCompleto_comp','-depsc')




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

