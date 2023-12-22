clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')

%%

Ts = 0.15;
z =tf('z',Ts);
s = tf('s');

%% caso 1
Ke = 1;
tau = 1.5;
G = Ke/(tau*s+1)
L1 = 10*Ts;

P = G*exp(-L1*s);

%%% ajuste PS - controlador primário
Czps = 3.70*(z-0.77)/(z-1);
Fzps = 0.23/(z-0.77);

Gnz = c2d(G,Ts,'zoh')
d = L1/Ts;
Lnz = z^-d;

Pr = Ke/(1*s+1);
Prz = c2d(Pr,Ts,'zoh');
dr = 12;
Prz = Prz*z^-dr;

dP = Prz/Gnz/Lnz-1;

%% ajuste PI
Cpi = 3.7*(z-0.77)/(z-1);

mf_pi = feedback(Cpi*Gnz,1);


%% ajuste GPC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

m = 1; % número de entradas da planta - manipuladas
n = 1; % número de saidas da planta - controladas

atrasoEntrada = L1/Ts;

Nu = 5; %horizontes de controle - dimensão 1 x m
N1 = 11; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = 25; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar
Nf = 40;% para gdmc somente;

delta = 1; %ponderação nos erros - dimensão 1 x n
lambda = 1; %ponderação nas ações de controle - dimensão 1 x m      

betaf = [0 0.3 0.6 0.8 0.9]; % para gdmc somente, polo do filtro

Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle


%% filtro equivalente PSF do GDMC
for kk = 1:size(betaf,2);
    modDegrauU = step(Gnz*Lnz,Ts:Ts:3*Nf*Ts);

    %%% Monta as matrizes offline que definem o problema DMC
    G = [];
    for i=1:Nu
        G(i:Ny,i) = modDegrauU(N1:N2-i+1);
    end

    Kgdmc = (G'*Qy*G+Qu)\G'*Qy; % utilizado apenas no caso irrestrito
    Kgdmc1 = Kgdmc(1,:);


    %%% cálculo dos filtros dos erros de predição (SISO) para o GDMC
    Fegdmc = tf(0,1,Ts);
    nf = 1;

    pz1 = pole(Gnz); % polo da FT da entrada

    for i=N1(1):N2(1)
        %%% monta sistema para obtenção dos parâmetros do filtro
        %%% primeira equação (z^i - F_i(z) = 0 -> para z=pz1
        %%% segunda equação (z^i - F_i(z) = 0 -> para z=pz2
        %%% terceira equação F_i(1) = 0 -> força ganho unitário
        indf = i-N1(1)+1;
        Af = [pz1 1;
              1 1];
        bf = [pz1^i*(pz1-betaf(kk));
              (1-betaf(kk))];
        X = Af\bf;
        Fegdmc(indf,1) = (X(1)*z+X(2))/(z-betaf(kk));
        %%% teste da condição
    %     pz1^i-(X(1)*pz1^2+X(2)*pz1)/(pz1-betaf)^2

        %%% armazena coeficientes gtil
        modDegrauUF{i} = filter(Fegdmc(indf,1).num{1},Fegdmc(indf,1).den{1},modDegrauU);

    end


    %%% calcula a matriz H para o cálculo da resposta livre no caso GDMC
    H1 = [];
    H2 = [];



    for i=N1(1):N2(1)
        H1 = [H1;modDegrauU(i+1:i+Nf)'];
        H2 = [H2;modDegrauUF{i}(1:Nf)'];

    end
    Hgdmc = H1-H2;

    zudmc = [];
    for k=1:size(Hgdmc,2)
        zudmc = [zudmc;z^(-k)];
    end


    %%% ceq normal do gdmc
    f0 = sum(Kgdmc1)
    
    Fke_num = 0;
    for i= 1:Ny
        Fke_num = Fke_num + Kgdmc1(i)*Fegdmc(i,1).num{1};
    end
    Fke = tf(Fke_num,Fegdmc(1).den{1},Ts)

    k1h = Kgdmc1*Hgdmc;

    k1hz = minreal(k1h*zudmc);
    Dc = (1-z^-1)*(1+k1hz);
    Cgdmc = minreal(Fke/Dc)

    mf_gdmc(kk) = feedback(Cgdmc*Gnz*Lnz,1);
    mf_gdmc2(kk) = 1/mf_gdmc(kk);

end

%% dados dos critérios de robustez
w = logspace(-1,log10(2*pi*1/2/Ts),200);
[mag,pha] = bode([mf_gdmc2,dP],w);

mag1 = reshape(mag(1,1,:),200,1);
mag2 = reshape(mag(1,2,:),200,1);
mag3 = reshape(mag(1,3,:),200,1);
mag4 = reshape(mag(1,4,:),200,1);
mag5 = reshape(mag(1,5,:),200,1);
mag6 = reshape(mag(1,6,:),200,1);

%%
cores = gray(7);
cores = cores(1:end-1,:);

hf = figure
h= subplot(1,1,1)
semilogx(w,20*log10(mag1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w,20*log10(mag2),'LineWidth',tamlinha,'Color',cores(2,:))
semilogx(w,20*log10(mag3),'LineWidth',tamlinha,'Color',cores(3,:))
semilogx(w,20*log10(mag4),'LineWidth',tamlinha,'Color',cores(4,:))
semilogx(w,20*log10(mag5),'LineWidth',tamlinha,'Color',cores(5,:))
semilogx(w,20*log10(mag6),':','LineWidth',tamlinha,'Color',cores(6,:))

semilogx(2*pi/2/Ts*[1 1],[-500 500],'LineWidth',tamlinha,'Color',cores(1,:))
grid on

hl = legend('z_f=0','z_f=0,3','z_f=0,6','z_f=0,8','z_f=0,9','\delta G(z)','Location','NorthWest')
ylabel('Magnitude (dB)','FontSize', tamletra)

% xlim([10^-2 3*10^1])
ylim([-30 35])

set(h, 'FontSize', tamletra);


xlabel('Frequência (rad/min)','FontSize', tamletra)


set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.1512 0.5446 0.1571 0.4468]; 
% print('trocador_PS_GPC_GDMC_bode2','-depsc')
