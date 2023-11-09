clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%

load('cstr_gdmc') %carrega os dados de simulação do GDMC1 e GDMC2, inclusive as sintonias


%% calculo do modelo linearizado
Ass = [0 0 0 ;
       0, -qi0/V-k0*exp(-ER/T0), -Ca0*k0*ER/(T0^2)*exp(-ER/T0);
       0,-dH/p/cp*k0*exp(-ER/T0), -qi0/V-(dH/p/cp)*Ca0*k0*ER/(T0^2)*exp(-ER/T0)];
Bss = [-1/Ac, 0 , 0;
       0, qi0/V, 0;
       0,0 -1/V];
Bqss = [1/Ac, 0;
        (Caf0-Ca0)/V, 0;
        (Ti0-T0)/V, qi0/V];

Css = eye(3);

Gss = ss(Ass,[Bss,Bqss],Css,[]);
Gs = tf(Gss);

Gssz = c2d(Gss,Ts,'zoh');
Gz0 = tf(Gssz);

Gz = Gz0(:,1:m);
Gzq = Gz0(:,m+1:end);

atrasoEntrada = totaldelay(Gz);


%% Controlador GPC - MIMO

%%% obtenção da representação MFD
[Amfd,Bmfd] = MFDredux(Gz.den,Gz.num);

%%%

Bqmfd = Bmfd(:,m+1:end);
Bmfd = Bmfd(:,1:m);

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

Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;

Kgpc = (G'*Qe*G+Qu)\G'*Qe;
Kgpc1 = [];
for i=1:m
    Kgpc1 = [Kgpc1; Kgpc(sum(Nu(1:i-1))+1,:)];
end



%% inicializacao dos estados e variaveis da simulação
saidas = repmat(x0,[1,nit]); % vetor da saída
entradas = repmat(u0,[1,nit]); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    saidas(:,k) = modeloCSTR(saidas(:,k-1),entradas(:,k-1),perts(:,k-1),Ts);
    
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
    %% Resolve o problema de otimização
%     du(:,k) = Kgpc1*(R-f);
    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasGPC = saidas;
entradasGPC = entradas;
duGPC = du;


%% simulação com o controlador SSMPC

A = Gssz.A;
B = Gssz.B(:,1:m);
C = Gssz.C;

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

%%% matrizes de ponderações
Qyi = diag(delta);

Qe = [];
for i=1:max(N2)
    Qe = blkdiag(Qe,Qyi);
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

indl = [];
for i=1:n
%     indl = [indl, (N1(i):N2(i))*n-1+(i-1)];
    indl = [indl, (N1(i)-1)*n+i:n:(N2(i)-1)*n+i];
end

indl = sort(indl);


F = F(indl,:);
Ii = Ii(indl,:);
G = G(indl,indc);

Qu = Qu(indc,indc);
Qe = Qe(indl,indl);

Kmpc = (G'*Qe*G+Qu)\(G'*Qe);
Kmpc1 = Kmpc(1:m,:);

Hqp = 2*(G'*Qe*G+Qu);
fqp1 = -2*G'*Qe;

%%% matrizes de restrições
I = eye(m);
T = [];
for i=1:max(Nu)
    T = [T;
          repmat(I,1,i),zeros(m,max(Nu)*m-i*m)];
end

T = T(indc,indc);

Rbar = [T;
        -T];


%% inicializacao dos estados e variaveis da simulação
estados = repmat(x0,[1,nit]);
saidas = repmat(x0,[1,nit]); % vetor da saída
entradas = repmat(u0,[1,nit]); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    estados(:,k) = modeloCSTR(estados(:,k-1),entradas(:,k-1),perts(:,k-1),Ts);
    saidas(:,k) = estados(:,k);
    
    %% -- Controlador GDMC 
    %%% referencias
    R = repmat(refs(:,k),max(N2),1);
    R = R(indl,:);
    
    %%% cálculo da resposta livre;
    f = Ii*saidas(:,k)...
        +F*[estados(:,k);
           estados(:,k-1)];

    %%% montagem das matrizes de restrição
    rbarmax = repmat([umax'-entradas(:,k-1)],max(Nu),1);
    rbarmin = repmat([-umin'+entradas(:,k-1)],max(Nu),1);
    rbar = [rbarmax(indc,:);
             rbarmin(indc,:)];
    
    fqp = fqp1*(R-f);
    %% Resolve o problema de otimização
    X = quadprog(Hqp,fqp,Rbar,rbar);
    du(:,k) = X(1:m,1);

    %     du(:,k) = Kmpc1*(R-f); % sem restrições

    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasSS = saidas;
entradasSS = entradas;
duSS = du;    



%% figuras
cores = gray(6);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin)*Ts;
ind = nin:nit;


hf = figure
h=subplot(3,2,1)
plot(t,saidasGDMC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGDMC2(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC(1,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,saidasSS(1,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(5,:))
grid on    
title('Controladas','FontSize', tamletra)
ylabel('h (m)','FontSize', tamletra)
xlim([0 600])
hl = legend('GDMC1','GDMC2','GPC','SSMPC','Referência','Location','SouthEast')
set(h, 'FontSize', tamletra);

h.YTickLabel = trocaponto(h.YTickLabel)


h=subplot(3,2,3)
plot(t,saidasGDMC(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGDMC2(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC(2,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,saidasSS(2,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(5,:))
ylabel('C_{a} (kmol/m^3)','FontSize', tamletra)
xlim([0 600])
set(h, 'FontSize', tamletra);
grid on
h.YTickLabel = trocaponto(h.YTickLabel)

h=subplot(3,2,5)
plot(t,saidasGDMC(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGDMC2(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC(3,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,saidasSS(3,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

plot(t,refs(3,ind),'-.','LineWidth',tamlinha,'Color',cores(5,:))
ylabel('T (K)','FontSize', tamletra)
xlabel('Tempo (segundos)','FontSize', tamletra)
xlim([0 600])
set(h, 'FontSize', tamletra);
grid on
h.YTickLabel = trocaponto(h.YTickLabel)


h=subplot(3,2,2)
plot(t,entradasGDMC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGDMC2(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC(1,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,entradasSS(1,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))


grid on
ylabel('q_o (m^3/s)','FontSize', tamletra)
xlim([0 600])
title('Manipuladas','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);

h.YTickLabel = trocaponto(h.YTickLabel)

h2=subplot(3,2,4)
plot(t,entradasGDMC(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGDMC2(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC(2,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,entradasSS(2,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

grid on
ylabel('C_{af} (kmol/m^3)','FontSize', tamletra)
xlim([0 600])
h2.YTick = [5 5.4];
set(h2,'YAxisLocation','right')
set(h2, 'FontSize', tamletra);

h2.YTickLabel = trocaponto(h2.YTickLabel)

h=subplot(3,2,6)
plot(t,entradasGDMC(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGDMC2(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC(3,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,entradasSS(3,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

grid on
ylabel('Q_h/(\rho c_p) (Km^3/s)','FontSize', tamletra)
xlim([0 600])
xlabel('Tempo (segundos)','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);

h.YTickLabel = trocaponto(h.YTickLabel)


% hf.Position = [680 558 560 420];

hf.Position = tamfigura3;
hl.Position = [0.3546 0.6373 0.2054 0.2310];
% hl2.Position = [0.8342 0.1739 0.1225 0.2952];

% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))

% print('cstr_ssmpc_comp_ref','-depsc')

%%
hf = figure
h=subplot(3,2,1)
plot(t,saidasGDMC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGDMC2(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,saidasSS(1,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(5,:))
grid on    
title('Controladas','FontSize', tamletra)
ylabel('h (m)','FontSize', tamletra)
xlim([600 900])
ylim([1.045 1.115])
set(h, 'FontSize', tamletra);

h.YTickLabel = trocaponto(h.YTickLabel)

h=subplot(3,2,3)
plot(t,saidasGDMC(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGDMC2(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,saidasSS(2,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(5,:))
ylabel('C_{a} (kmol/m^3)','FontSize', tamletra)
xlim([600 900])
ylim([0.98 1.07])
set(h, 'FontSize', tamletra);
grid on

h.YTickLabel = trocaponto(h.YTickLabel)

h=subplot(3,2,5)
plot(t,saidasGDMC(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGDMC2(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC(3,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,saidasSS(3,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

plot(t,refs(3,ind),'-.','LineWidth',tamlinha,'Color',cores(5,:))
ylabel('T (K)','FontSize', tamletra)
xlabel('Tempo (segundos)','FontSize', tamletra)
xlim([600 900])
ylim([385 405])
set(h, 'FontSize', tamletra);
grid on

h.YTickLabel = trocaponto(h.YTickLabel)

hl = legend('GDMC1','GDMC2','GPC','SSMPC','Referência','Location','SouthEast')


h=subplot(3,2,2)
plot(t,entradasGDMC(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGDMC2(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,entradasSS(1,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

h.YTickLabel = trocaponto(h.YTickLabel)

grid on
ylabel('q_o (m^3/s)','FontSize', tamletra)
xlim([600 900])
xlim([600 900])
ylim([0.00495 0.00535])
title('Manipuladas','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);


h2=subplot(3,2,4)
plot(t,entradasGDMC(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGDMC2(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,entradasSS(2,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

grid on
ylabel('C_{af} (kmol/m^3)','FontSize', tamletra)
xlim([600 900])
ylim([4.7 5.02])
h2.YTick = [5 5.4];
set(h2,'YAxisLocation','right')
set(h2, 'FontSize', tamletra);
h2.YTickLabel = trocaponto(h2.YTickLabel)

h=subplot(3,2,6)
plot(t,entradasGDMC(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGDMC2(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC(3,ind),':','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,entradasSS(3,ind),'-','LineWidth',tamlinha,'Color',cores(4,:))

grid on
ylabel('Q_h/(\rho c_p) (Km^3/s)','FontSize', tamletra)
xlim([600 900])
ylim([0.63 0.76])
xlabel('Tempo (segundos)','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);

h.YTickLabel = trocaponto(h.YTickLabel)


% hf.Position = [680 558 560 420];

hf.Position = tamfigura3;
hl.Position = [0.3565 0.6873 0.2054 0.2167];
% hl2.Position = [0.8342 0.1739 0.1225 0.2952];

% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))

% print('cstr_ssmpc_comp_pert','-depsc')




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
