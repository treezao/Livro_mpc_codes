%%% Simulação comparativa do PNMPC, NEPSAC e xxx
clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')

import casadi.*
%% parametros 
%%% Modelo CSTR listado em DU et al 
nx = 2; %numero de estados
n = 1; % numero de saidas
m = 1; % numero de entradas
q = 1; % numero de perturbacoes

na = 1; % numero de amostras de y atrasadas usadas no calculo da saida
nb = 1; % numero de amostras de u atrasadas usadas no calculo da saida
nq = 1; % numero de amostras de q atrasadas usadas no calculo da saida

dnominal = 10;
dreal = 11;

umax = [0.5];
umin = [-0.5];

dumax = [.1];
dumin = -dumax;

y0 = 0;
u0 = -0.2500;
q0 = 0;

z = tf('z');
%% Sintonia PNMPC

N1 = 1;
N2 = 5; 
Ny = N2-N1+1;
    
Nu = [3]; % vetor contendo os horizontes de controle

delta = [1]; % vetor contendo as ponderacoes de saida
lambda = [5];

alpha = 0;

lf = 0
lf1 = 0.9
Fe = (1-lf)*z/(z-lf)

nafe = size(Fe.den{1},2)-1;

%% configuracaoes de simulacao
nin = 2+max(dnominal,dreal); %%% iteracao inicial
nit = nin+100; %%% iteracao final

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



%% inicializacao de algumas variaveis
saida = zeros(n,nit);
saida(1,1:nin) = y0;

entrada = zeros(m,nit);
entrada(1,:) = u0;

refs = zeros(n,nit+max(N2));
refs(1,:) = y0;
refs(1,nin+10:end) = 1;

perts = zeros(q,nit);
perts(1,:) = 0;
perts(1,nin+60:nit) = 0.1;


eta = zeros(n,nit+max(N2));
e = zeros(n,nit+max(N2));
saidaModelo = zeros(n,nit);


saidaPSF_nom = zeros(n,nit);
saidaPSF_nom(1,1:nin) = y0;

saidaPSF = zeros(n,nit);
saidaPSF(1,1:nin) = y0;

erroPSF = zeros(n,nit);
erroPSF_filtrado = zeros(n,nit);

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
brestLB0 = [];
brestUB0 = [];
Arest = [];
for i=1:m
    LBdu = [LBdu;ones(Nu(i),1)*dumin(i)];
    UBdu = [UBdu;ones(Nu(i),1)*dumax(i)];
    
    T = tril(ones(Nu(i)));
    Arest = blkdiag(Arest,T);
    
    brestLB0 = [brestLB0;ones(Nu(i),1)*umin(i)];
    brestUB0 = [brestUB0;ones(Nu(i),1)*umax(i)];

end

X = zeros(sum(Nu),1);

%% simulação 

for k=nin:nit
    %% modelo do CSTR
    saida(:,k) = modeloProcesso(saida(:,k-1:-1:k-2),entrada(:,k-dreal-1),perts(:,k-1));

    %% preditor de smith NL (NLFSP)
    saidaPSF_nom(:,k) = modeloNominal(saidaPSF_nom(:,k-1:-1:k-2),entrada(:,k-1),q0);
    erroPSF(:,k) = saida(:,k)-saidaPSF_nom(:,k-dnominal);
    
    for i=1:n
        erroPSF_filtrado(i,k) = -Fe.den{1}(2:end)*erroPSF_filtrado(i,k-1:-1:k-nafe)'+Fe.num{1}*erroPSF(i,k:-1:k-nafe)';
    end
    
    saidaPSF(:,k) = saidaPSF_nom(:,k) + erroPSF_filtrado(:,k);
    
    
    
    
    %% controlador PNMPC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(:,k) = modeloNominal(saidaPSF(:,k-1:-1:k-2),entrada(:,k-1),q0);
    %%% calculo dos etas futuros
    eta(:,k) = saidaPSF(:,k) - saidaModelo(:,k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:max(N2)
        eta(:,k+j) = -D(2:end)*eta(:,k+j-1:-1:k+j-nd+1)' + C*e(:,k+j:-1:k+j-nc+1)';
    end
    
    
    %%% cálculo da resposta livre corrigida
    yc = zeros(n,N2);
    yc(:,1) = modeloNominal(saidaPSF(:,k:-1:k-1),entrada(:,k-1),q0) + eta(:,k+1);
    yc(:,2) = modeloNominal([yc(:,1),saidaPSF(:,k)],entrada(:,k-1),q0) + eta(:,k+2);

    for i=3:N2
        yc(:,i) = modeloNominal(yc(:,i-1:-1:i-2),entrada(:,k-1),q0) + eta(:,k+i);
    end
    
    %%% passando para o formato vetorial correto
    f = zeros(sum(N2-N1)+1,1);
    for i=1:n
        f(sum(Ny(1:i-1))+1:sum(Ny(1:i)),1) = yc(i,N1(i):N2(i))';
    end

    
    %%% calculo de G
    % calculo de y0
    y00 = zeros(n,N2);
    y00(:,1) = modeloNominal(saidaPSF(:,k:-1:k-1),entrada(:,k-1),q0);
    y00(:,2) = modeloNominal([y00(:,1), saidaPSF(:,k)],entrada(:,k-1),q0);
    for i=3:N2
        y00(:,i) = modeloNominal(y00(:,i-1:-1:i-2),entrada(:,k-1),q0);
    end
    
    %%% passando para o formato vetorial correto
    y0v = zeros(sum(N2-N1)+1,1);
    for i=1:n
        y0v(sum(Ny(1:i-1))+1:sum(Ny(1:i)),1) = y00(i,N1(i):N2(i))';
    end

    %%% calculo de yi para obter G
    G = zeros(sum(Ny),sum(Nu));
    for i=1:m
        u00 = entrada(:,k-1);
        
        if(entrada(i,k-1) == 0)
            epsilon = 0.001;
        else
            epsilon = entrada(i,k-1)/1000;
        end

        for j=1:Nu
            yi = zeros(n,N2);

            ui = repmat(u00,[1 N2]);
            ui(i,j:end) = ui(i,j:end)+epsilon;

            yi(:,1) = modeloNominal(saidaPSF(:,k:-1:k-1),ui(:,1),q0);
            yi(:,2) = modeloNominal([yi(:,1),saidaPSF(:,k)],ui(:,2),q0);
            for jj=3:N2
                yi(:,jj) = modeloNominal(yi(:,jj-1:-1:jj-2),ui(:,jj),q0);
            end
            
            %%% passando para o formato vetorial correto
            yiv = zeros(sum(N2-N1)+1,1);
            for jj=1:n
                yiv(sum(Ny(1:jj-1))+1:sum(Ny(1:jj)),1) = yi(jj,N1(jj):N2(jj))';
            end
            
            G(:,sum(Nu(1:i-1))+j) = (yiv-y0v)/epsilon;
            
        end
    end

    %%% Referência futuras
    rf = alpha*rfant + (1-alpha)*refs(:,k);
    rfant = rf;
    rfut = [];
    for i=1:n
        rfut = [rfut;
                rf(i,1)*ones(N2(i)-N1(i)+1,1)];
    end

    Hoti = 2*(G'*Qe*G+Qu);
    foti = -2*G'*Qe*(rfut-f);
    
    temp1 = [];
    temp2 = [];
    for i=1:m
        temp1 = [temp1;-entrada(i,k-1)*ones(Nu(i),1)];
        temp2 = [temp2;-entrada(i,k-1)*ones(Nu(i),1)];
    end
    
    brestLB = brestLB0 + temp1;
    brestUB = brestUB0 + temp2;
    
    tic
    x_opt = QPsolver('h', Hoti, 'g', foti);
%     x_opt = QPsolver('h', Hoti, 'g', foti,'x0',X,'lbx', LBdu, 'ubx', UBdu,'a', Arest, 'lba',brestLB,'uba',brestUB);
    x2(k) = toc;
    
    X = full(x_opt.x);
%     X = quadprog(Hoti,foti,[Arest;-Arest],[brestUB;-brestLB],[],[],LBdu,UBdu);

    for i=1:m
        ind = sum(Nu(1:i-1));
        du(i,k) = X(ind+1);
    end
        
    entrada(:,k) = entrada(:,k-1)+du(:,k); 
    

end

saidaPNMPC_PSF = saida;
entradaPNMPC_PSF = entrada;
duPNMPC_PSF= du;

%% PSF com filtro diferente de 1
Fe = (1-lf1)*z/(z-lf1)

%%% inicializacao de algumas variaveis
saida = zeros(n,nit);
saida(1,1:nin) = y0;

entrada = zeros(m,nit);
entrada(1,:) = u0;

eta = zeros(n,nit+max(N2));
e = zeros(n,nit+max(N2));
saidaModelo = zeros(n,nit);

saidaPSF_nom = zeros(n,nit);
saidaPSF_nom(1,1:nin) = y0;

saidaPSF = zeros(n,nit);
saidaPSF(1,1:nin) = y0;

erroPSF = zeros(n,nit);
erroPSF_filtrado = zeros(n,nit);

du = zeros(m,nit);

rfant = saida(:,1);

%% simulação 

for k=nin:nit
    %% modelo do CSTR
    saida(:,k) = modeloProcesso(saida(:,k-1:-1:k-2),entrada(:,k-dreal-1),perts(:,k-1));

    %% preditor de smith NL (NLFSP)
    saidaPSF_nom(:,k) = modeloNominal(saidaPSF_nom(:,k-1:-1:k-2),entrada(:,k-1),q0);
    erroPSF(:,k) = saida(:,k)-saidaPSF_nom(:,k-dnominal);
    
    for i=1:n
        erroPSF_filtrado(i,k) = -Fe.den{1}(2:end)*erroPSF_filtrado(i,k-1:-1:k-nafe)'+Fe.num{1}*erroPSF(i,k:-1:k-nafe)';
    end
    
    saidaPSF(:,k) = saidaPSF_nom(:,k) + erroPSF_filtrado(:,k);
    
    
    
    
    %% controlador PNMPC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(:,k) = modeloNominal(saidaPSF(:,k-1:-1:k-2),entrada(:,k-1),q0);
    %%% calculo dos etas futuros
    eta(:,k) = saidaPSF(:,k) - saidaModelo(:,k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:max(N2)
        eta(:,k+j) = -D(2:end)*eta(:,k+j-1:-1:k+j-nd+1)' + C*e(:,k+j:-1:k+j-nc+1)';
    end
    
    
    %%% cálculo da resposta livre corrigida
    yc = zeros(n,N2);
    yc(:,1) = modeloNominal(saidaPSF(:,k:-1:k-1),entrada(:,k-1),q0) + eta(:,k+1);
    yc(:,2) = modeloNominal([yc(:,1),saidaPSF(:,k)],entrada(:,k-1),q0) + eta(:,k+2);

    for i=3:N2
        yc(:,i) = modeloNominal(yc(:,i-1:-1:i-2),entrada(:,k-1),q0) + eta(:,k+i);
    end
    
    %%% passando para o formato vetorial correto
    f = zeros(sum(N2-N1)+1,1);
    for i=1:n
        f(sum(Ny(1:i-1))+1:sum(Ny(1:i)),1) = yc(i,N1(i):N2(i))';
    end

    
    %%% calculo de G
    % calculo de y0
    y00 = zeros(n,N2);
    y00(:,1) = modeloNominal(saidaPSF(:,k:-1:k-1),entrada(:,k-1),q0);
    y00(:,2) = modeloNominal([y00(:,1), saidaPSF(:,k)],entrada(:,k-1),q0);
    for i=3:N2
        y00(:,i) = modeloNominal(y00(:,i-1:-1:i-2),entrada(:,k-1),q0);
    end
    
    %%% passando para o formato vetorial correto
    y0v = zeros(sum(N2-N1)+1,1);
    for i=1:n
        y0v(sum(Ny(1:i-1))+1:sum(Ny(1:i)),1) = y00(i,N1(i):N2(i))';
    end

    %%% calculo de yi para obter G
    G = zeros(sum(Ny),sum(Nu));
    for i=1:m
        u00 = entrada(:,k-1);
        
        if(entrada(i,k-1) == 0)
            epsilon = 0.001;
        else
            epsilon = entrada(i,k-1)/1000;
        end

        for j=1:Nu
            yi = zeros(n,N2);

            ui = repmat(u00,[1 N2]);
            ui(i,j:end) = ui(i,j:end)+epsilon;

            yi(:,1) = modeloNominal(saidaPSF(:,k:-1:k-1),ui(:,1),q0);
            yi(:,2) = modeloNominal([yi(:,1),saidaPSF(:,k)],ui(:,2),q0);
            for jj=3:N2
                yi(:,jj) = modeloNominal(yi(:,jj-1:-1:jj-2),ui(:,jj),q0);
            end
            
            %%% passando para o formato vetorial correto
            yiv = zeros(sum(N2-N1)+1,1);
            for jj=1:n
                yiv(sum(Ny(1:jj-1))+1:sum(Ny(1:jj)),1) = yi(jj,N1(jj):N2(jj))';
            end
            
            G(:,sum(Nu(1:i-1))+j) = (yiv-y0v)/epsilon;
            
        end
    end

    %%% Referência futuras
    rf = alpha*rfant + (1-alpha)*refs(:,k);
    rfant = rf;
    rfut = [];
    for i=1:n
        rfut = [rfut;
                rf(i,1)*ones(N2(i)-N1(i)+1,1)];
    end

    Hoti = 2*(G'*Qe*G+Qu);
    foti = -2*G'*Qe*(rfut-f);
    
    temp1 = [];
    temp2 = [];
    for i=1:m
        temp1 = [temp1;-entrada(i,k-1)*ones(Nu(i),1)];
        temp2 = [temp2;-entrada(i,k-1)*ones(Nu(i),1)];
    end
    
    brestLB = brestLB0 + temp1;
    brestUB = brestUB0 + temp2;
    
    tic
    x_opt = QPsolver('h', Hoti, 'g', foti);
%     x_opt = QPsolver('h', Hoti, 'g', foti,'x0',X,'lbx', LBdu, 'ubx', UBdu,'a', Arest, 'lba',brestLB,'uba',brestUB);
    x2(k) = toc;
    
    X = full(x_opt.x);
%     X = quadprog(Hoti,foti,[Arest;-Arest],[brestUB;-brestLB],[],[],LBdu,UBdu);

    for i=1:m
        ind = sum(Nu(1:i-1));
        du(i,k) = X(ind+1);
    end
        
    entrada(:,k) = entrada(:,k-1)+du(:,k); 
    

end

saidaPNMPC_PSF1 = saida;
entradaPNMPC_PSF1 = entrada;
duPNMPC_PSF1 = du;


%% PNMPC original
N1 = dnominal +N1;
N2 = dnominal + N2;

%%% inicializacao de algumas variaveis
saida = zeros(n,nit);
saida(1,1:nin) = y0;

entrada = zeros(m,nit);
entrada(1,:) = u0;

eta = zeros(n,nit+max(N2));
e = zeros(n,nit+max(N2));
saidaModelo = zeros(n,nit);

du = zeros(m,nit);

rfant = saida(:,1);

%% simulação 

for k=nin:nit
    %% modelo do CSTR
    saida(:,k) = modeloProcesso(saida(:,k-1:-1:k-2),entrada(:,k-dreal-1),perts(:,k-1));

    %% controlador PNMPC
    %%% cálculo da resposta do modelo nominal
    saidaModelo(:,k)= modeloNominal(saida(:,k-1:-1:k-2),entrada(:,k-dnominal-1),q0);
    %%% calculo dos etas futuros
    eta(:,k) = saida(:,k) - saidaModelo(:,k); % cálculo do eta atual
    
    e(k) = -C(2:end)*e(k-1:-1:k-nc+1)'+ D*eta(k:-1:k-nd+1)';
    
    for j=1:N2
        eta(:,k+j) = -D(2:end)*eta(:,k+j-1:-1:k+j-nd+1)' + C*e(:,k+j:-1:k+j-nc+1)';
    end
    
    
    %%% cálculo da resposta livre corrigida
    yc = zeros(n,N2);
    yc(:,1) = modeloNominal(saida(:,k:-1:k-1),entrada(:,k-dnominal),q0) + eta(:,k+1);
    yc(:,2) = modeloNominal([yc(:,1),saida(:,k)],entrada(:,k-dnominal+1),q0) + eta(:,k+2);

    for i=3:dnominal
        yc(:,i) = modeloNominal(yc(:,i-1:-1:i-2),entrada(:,k-dnominal-1+i),q0) + eta(:,k+i);
    end
    for i=dnominal+1:N2
        yc(:,i) = modeloNominal(yc(:,i-1:-1:i-2),entrada(:,k-1),q0) + eta(:,k+i);
    end
    
    %%% passando para o formato vetorial correto
    f = zeros(sum(N2-N1)+1,1);
    for i=1:n
        f(sum(Ny(1:i-1))+1:sum(Ny(1:i)),1) = yc(i,N1(i):N2(i))';
    end

    
    %%% calculo de G
    % calculo de y0
    y0 = zeros(n,N2);
    y0(:,1) = modeloNominal(saida(:,k:-1:k-1),entrada(:,k-dnominal),q0);
    y0(:,2) = modeloNominal([y0(:,1), saida(:,k)],entrada(:,k-dnominal+1),q0);
    for i=3:dnominal
        y0(:,i) = modeloNominal(y0(:,i-1:-1:i-2),entrada(:,k-dnominal-1+i),q0);
    end
    for i=dnominal+1:N2
        y0(:,i) = modeloNominal(y0(:,i-1:-1:i-2),entrada(:,k-1),q0);
    end
    
    %%% passando para o formato vetorial correto
    y0v = zeros(sum(N2-N1)+1,1);
    for i=1:n
        y0v(sum(Ny(1:i-1))+1:sum(Ny(1:i)),1) = y0(i,N1(i):N2(i))';
    end

    %%% calculo de yi para obter G
    G = zeros(sum(Ny),sum(Nu));
    for i=1:m
        u00 = entrada(:,k-1);
        
        if(entrada(i,k-1) == 0)
            epsilon = 0.001;
        else
            epsilon = entrada(i,k-1)/1000;
        end

        for j=1:Nu
            yi = zeros(n,N2);

            ui = repmat(u00,[1 N2]);
            ui(i,j:end) = ui(i,j:end)+epsilon;

            yi(:,1) = modeloNominal(saida(:,k:-1:k-1),entrada(:,k-dnominal),q0);
            yi(:,2) = modeloNominal([yi(:,1),saida(:,k)],entrada(:,k-dnominal+1),q0);
            for jj=3:dnominal
                yi(:,jj) = modeloNominal(yi(:,jj-1:-1:jj-2),entrada(:,k-dnominal-1+jj),q0);
            end
            for jj=dnominal+1:N2
                yi(:,jj) = modeloNominal(yi(:,jj-1:-1:jj-2),ui(:,jj-dnominal),q0);
            end
            
            %%% passando para o formato vetorial correto
            yiv = zeros(sum(N2-N1)+1,1);
            for jj=1:n
                yiv(sum(Ny(1:jj-1))+1:sum(Ny(1:jj)),1) = yi(jj,N1(jj):N2(jj))';
            end
            
            G(:,sum(Nu(1:i-1))+j) = (yiv-y0v)/epsilon;
            
        end
    end

    %%% Referência futuras
    rf = alpha*rfant + (1-alpha)*refs(:,k);
    rfant = rf;
    rfut = [];
    for i=1:n
        rfut = [rfut;
                rf(i,1)*ones(N2(i)-N1(i)+1,1)];
    end

    Hoti = 2*(G'*Qe*G+Qu);
    foti = -2*G'*Qe*(rfut-f);
    
    temp1 = [];
    temp2 = [];
    for i=1:m
        temp1 = [temp1;-entrada(i,k-1)*ones(Nu(i),1)];
        temp2 = [temp2;-entrada(i,k-1)*ones(Nu(i),1)];
    end
    
    brestLB = brestLB0 + temp1;
    brestUB = brestUB0 + temp2;
    
    tic
    x_opt = QPsolver('h', Hoti, 'g', foti);
%     x_opt = QPsolver('h', Hoti, 'g', foti,'x0',X,'lbx', LBdu, 'ubx', UBdu,'a', Arest, 'lba',brestLB,'uba',brestUB);
    x2(k) = toc;
    
    X = full(x_opt.x);
%     X = quadprog(Hoti,foti,[Arest;-Arest],[brestUB;-brestLB],[],[],LBdu,UBdu);

    for i=1:m
        ind = sum(Nu(1:i-1));
        du(i,k) = X(ind+1);
    end
        
    entrada(:,k) = entrada(:,k-1)+du(:,k); 
    

end

saidaPNMPC = saida;
entradaPNMPC = entrada;
duPNMPC = du;


%% Geração das figuras
cores = gray(5);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin);
ind = nin:nit;


hf= figure
hf.Position = tamfigura;
h=subplot(2,1,1)
plot(t,saidaPNMPC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidaPNMPC_PSF(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidaPNMPC_PSF1(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))



hl = legend('PNMPC','DTC-PNMPC F_e(z)=1','DTC-PNMPC F_e(z)\neq 1','Referência','Location','SouthEast')

% xlim([110 150])
ylim([-0.05 1.65])

ylabel('Controladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h=subplot(2,1,2)
plot(t,entradaPNMPC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradaPNMPC_PSF(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradaPNMPC_PSF1(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))


% xlim([110 150])
ylim([-0.3 0.4])


ylabel('Manipuladas','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

xlabel('Tempo (amostras)')


hl.Position = [0.6613 0.4992 0.3321 0.2758];

% print('dtc_nlin_pnmpc_plucenio_psf_final_erro','-depsc')




%%
%%% Modelo CSTR MIMO completo livro Bequette Process Dynamics
function [y] = modeloProcesso(yk,u,q)
    y = 2.5*yk(1)*yk(2)/(1+yk(1)^2+yk(2)^2)+1.2*(u(1)+q(1))+0.3*cos(0.5*(yk(1)+yk(2)));
end

%% Modelo CSTR nominal utilizado pelos controladores
function [y] = modeloNominal(yk,u,q)
    y = 2.5*yk(1)*yk(2)/(1+yk(1)^2+yk(2)^2)+1.2*(u(1)+q(1))+0.3*cos(0.5*(yk(1)+yk(2)));
end
