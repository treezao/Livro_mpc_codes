clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%
cstrTito_integrador_gdmc_novo

saidas1 = saidas;
entradas1 = entradas;
du1 = du;

%% PARAMETROS DE AJUSTE DO CONTROLADOR

betaf = 0.4*[1 1 1]; %% polos dos filtros dos erros de predição por saída

%% definição do cenário de simulação

saidas = repmat(x0,[1,nit]); % vetor da saída
entradas = repmat(u0,[1,nit]); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle



%% Controlador GDMC - MIMO

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



%%% para as duas outras saídas basta remover o polo lento
pz = 0.8187;

for l=2:n
    F{l} = tf(0,1,Ts);
    
    for i=N1(l):N2(l)
        %%% monta sistema para obtenção dos parâmetros do filtro
        %%% primeira equação (z^i - F_i(z)) = 0 -> para z=pz
        %%% segunda equação F_i(1) = 0 -> força ganho unitário
        indf = i-N1(1)+1;
        Af = [1 1;
              pz^2, pz];
        bf = [(1-betaf(l))^2;
              pz^i*(pz-betaf(l))^2];
        X = Af\bf;
        F{l}(indf,1) = (X(1)*z^2+X(2)*z)/(z-betaf(l))^2;
        %%% teste da condição
%         zero(z^i - F{l}(indf,1))

%         F{l}(indf,1) = z^2/z^2;
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




% hf.Position = [680 558 560 420];

% hf.Position = tamfigura;
% hl.Position = [0.8056 0.5779 0.1811 0.3899];
% hl2.Position = [0.8342 0.1739 0.1225 0.2952];

% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))

% print('cstrTito_integrador_gdmc_novo_ref2','-depsc')

%%
hf = figure
h=subplot(3,2,1)
plot(t,saidas1(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidas(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
grid on    
title('Controladas','FontSize', tamletra)
ylabel('h (m)','FontSize', tamletra)
xlim([600 900])
ylim([1.045 1.115])
hl = legend('GDMC1','GDMC2','Referência','Location','NorthEast')
set(h, 'FontSize', tamletra);


h=subplot(3,2,3)
plot(t,saidas1(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
hold on
plot(t,saidas(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylabel('C_{a} (kmol/m^3)','FontSize', tamletra)
xlim([600 900])
ylim([0.98 1.07])
set(h, 'FontSize', tamletra);


h=subplot(3,2,5)
plot(t,saidas1(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
grid on
hold on
plot(t,saidas(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,refs(3,ind),'-.','LineWidth',tamlinha,'Color',cores(3,:))
ylabel('T (K)','FontSize', tamletra)
xlabel('Tempo (segundos)','FontSize', tamletra)
xlim([600 900])
ylim([385 405])
set(h, 'FontSize', tamletra);



h=subplot(3,2,2)
plot(t,entradas1(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradas(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))

grid on
ylabel('q_o (m^3/s)','FontSize', tamletra)
xlim([600 900])
ylim([0.00495 0.00535])
title('Manipuladas','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);


h2=subplot(3,2,4)
plot(t,entradas1(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradas(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))

grid on
ylabel('C_{af} (kmol/m^3)','FontSize', tamletra)
xlim([600 900])
ylim([4.7 5.02])
% h2.YTick = [5 5.4];
set(h2,'YAxisLocation','right')
set(h2, 'FontSize', tamletra);


h=subplot(3,2,6)
plot(t,entradas1(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradas(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))

grid on
ylabel('Q_h/(\rho c_p) (Km^3/s)','FontSize', tamletra)
xlim([600 900])
ylim([0.63 0.76])
xlabel('Tempo (segundos)','FontSize', tamletra)
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);




% hf.Position = [680 558 560 420];

% hf.Position = tamfigura;
hl.Position = [0.3547 0.7683 0.1768 0.1476];
% hl2.Position = [0.8342 0.1739 0.1225 0.2952];

% ax = get(gcf,'children');
% ind = find(isgraphics(ax,'Legend'));
% set(gcf,'children',ax([ind:end,1:ind-1]))

% print('cstrTito_integrador_gdmc_novo_pert2','-depsc')


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
