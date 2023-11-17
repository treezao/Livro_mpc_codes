clear all
close all
clc
addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%

% MIMO GPC - fracionador

m = 3; % número de entradas da planta - manipuladas
n = 3; % número de saidas da planta - controladas

s = tf('s');
% Definir os modelos que relacionam {saida x,entrada y};
G(1,1) = 4.05/(50*s+1)*exp(-27*s);
G(1,2) = 1.77/(60*s+1)*exp(-28*s);
G(1,3) = 5.88/(50*s+1)*exp(-27*s);
G(2,1) = 5.39/(50*s+1)*exp(-18*s);
G(2,2) = 5.72/(60*s+1)*exp(-14*s);
G(2,3) = 6.9/(40*s+1)*exp(-15*s);
G(3,1) = 4.38/(33*s+1)*exp(-20*s);
G(3,2) = 4.42/(44*s+1)*exp(-22*s);
G(3,3) = 7.2/(19*s+1);

% Definir os modelos que relacionam {saida x,perturbação z};

Ts = 1.0; %período de amostragem, OBS: manter atrasos múltiplos inteiros do Ts

%%% discretização do processo
Gz = c2d(G,Ts,'zoh');
z = tf('z',Ts);

atrasoEntrada = totaldelay(Gz);

dmin = min(atrasoEntrada')';
dminM = repmat(dmin,[1 m]);

atrasoGn = atrasoEntrada-dminM;

for i=1:n
    Lz(i,i) = 1/z^dmin(i);
end

%%% montagem de Gnz com os atrasos embutidos na FT
Gztemp = Gz;
Gztemp.inputdelay = 0;
Gztemp.outputdelay = 0;
Gztemp.iodelay = 0;
for i=1:n
    for j=1:m
        Gnz(i,j) = Gztemp(i,j)/z^atrasoGn(i,j);
    end
end

%%% definição do modelo com erro de modelagem
% Gr(1,1) = 4.05/(50*s+1)*exp(-27*s);
% Gr(1,2) = 1.77/(60*s+1)*exp(-28*s);
% Gr(1,3) = 5.88/(50*s+1)*exp(-27*s);
% Gr(2,1) = 5.39/(50*s+1)*exp(-18*s);
% Gr(2,2) = 5.72/(60*s+1)*exp(-14*s);
% Gr(2,3) = 6.9/(40*s+1)*exp(-15*s);
% Gr(3,1) = 4.38/(33*s+1)*exp(-20*s);
% Gr(3,2) = 4.42/(44*s+1)*exp(-22*s);
% Gr(3,3) = 7.2/(19*s+1);

Gr(1,1) = 4.05/(55*s+1)*exp(-27*s);
Gr(1,2) = 1.77/(54*s+1)*exp(-28*s);
Gr(1,3) = 5.88/(52*s+1)*exp(-27*s);
Gr(2,1) = 5.39/(50*s+1)*exp(-18*s);
Gr(2,2) = 5.72/(60*s+1)*exp(-14*s);
Gr(2,3) = 6.9/(40*s+1)*exp(-15*s);
Gr(3,1) = 4.38/(33*s+1)*exp(-20*s);
Gr(3,2) = 4.42/(44*s+1)*exp(-22*s);
Gr(3,3) = 7.2/(19*s+1);

Grz = c2d(Gr,Ts,'zoh');

atrasoEntradaR = totaldelay(Grz);

%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [5 5 5]; %horizontes de controle - dimensão 1 x m
N1 = [1 1 1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = [61 61 45]; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar

delta = [1 1 1]./Ny; %ponderação nos erros - dimensão 1 x n
lambda = [1 1 1]./Nu; %ponderação nas ações de controle - dimensão 1 x m       

Fe(1,1) = z/z;
Fe(2,2) = z/z;
Fe(3,3) = z/z;

nfe1 = size(Fe(1,1).den{1},2)-1;
nfe2 = size(Fe(2,2).den{1},2)-1;
nfe3 = size(Fe(3,3).den{1},2)-1;

usarRest = true; % flag para usar ou não restrição
umax = [0.8 0.8 0.8];
% umax = [10 10 10];
umin = -umax;

dumax = [0.1 0.1 0.1];
% dumax = [10 10 10];
dumin = -dumax;


Spsf = (eye(n)-Fe*Lz)*Gnz;

%% definição do cenário de simulação
nit =floor(450/Ts); %tempo de simulação
nin = 50; %número da iteração inicial
nit = nit+nin;


refs = zeros(n,nit);
perts = zeros(3,nit);

refs(1,nin+floor(10/Ts):end) = 1;
refs(2,nin+floor(100/Ts):end) = 0.5;
refs(3,nin+floor(170/Ts):end) = 0.2;

% perts(1,nin+floor(220/Ts):end) = 0.2;
% perts(2,nin+floor(350/Ts):end) = 0.5;
% perts(3,nin+floor(450/Ts):end) = -0.1;


%% Controlador GPC - MIMO

%%% obtenção da representação MFD

[Amfd,Bmfd] = MFDredux(Gz.den,Gz.num);


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
        Btil = conv(Bmfd{i,j}(2:end),[zeros(1,atrasoGn(i,j)) 1]); 
        
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
            H{i,j} = [H{i,j} zeros(Ny(i),nH(j)-nHi)];
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
Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

%%% matrizes do problema de otimização
Hoti = 2*(G'*Qy*G+Qu);
foti1 = -2*G'*Qy;

LB = [];
UB = [];
Aoti = [];
for i=1:m
    LB = [LB;repmat(dumin(i),Nu(i),1)];
    UB = [UB;repmat(dumax(i),Nu(i),1)];
    Aoti = blkdiag(Aoti,tril(ones(Nu(i))));
end
Aoti = [Aoti;-Aoti];


%% inicializacao dos estados e variaveis da simulação
% - obtencao dos numerados e denominadores
numR = cell(n,m);
denR = cell(n,m);

numS = cell(n,m);
denS = cell(n,m);

for i=1:n
	% obtencao das respostas ao degrau das entradas
	for j=1:m
		% obtencao dos numeradores e denominadores para o simulador
        [numR{i}{j}, denR{i}{j}] = tfdata(Grz(i,j),'v');
        [numS{i}{j}, denS{i}{j}] = tfdata(Spsf(i,j),'v');
    end
    
end


estadosEntr = cell(n,m);
estadosS_PSF = cell(n,m);
for i=1:n
    for j=1:m
		estadosEntr{i}{j} = zeros(1,nit);
        estadosS_PSF{i}{j} = zeros(1,nit);

    end
end

saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

saidasS_PSF = zeros(n,nit);
saidasFe = zeros(n,nit);
saidaPSF = zeros(n,nit);

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numR{i}{j},2);
            sizeD = size(denR{i}{j},2);

            str1 = numR{i}{j}*(entradas(j,k-atrasoEntradaR(i,j):-1:k+1-sizeN-atrasoEntradaR(i,j))'...
                               +perts(j,k-atrasoEntradaR(i,j):-1:k+1-sizeN-atrasoEntradaR(i,j))');
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end
        
        saidas(i,k) = aux ;%;+ 0.01*randn(1,1);
    end
    %% PSF
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numS{i}{j},2);
            sizeD = size(denS{i}{j},2);

            str1 = numS{i}{j}*(entradas(j,k:-1:k-sizeN+1)');
            str2 = -denS{i}{j}(1,2:end)*estadosS_PSF{i}{j}(k-1:-1:k-sizeD+1)';
            estadosS_PSF{i}{j}(k) = str1+str2;

            aux = aux + estadosS_PSF{i}{j}(k);
        end
        saidasS_PSF(i,k) = aux ;
        

        
        sizeN = size(Fe(i,i).num{1},2);
        
        saidasFe(i,k) = -Fe(i,i).den{1}(2:end)*saidasFe(i,k-1:-1:k-sizeN+1)' ...
                          +Fe(i,i).num{1}*saidas(i,k:-1:k-sizeN+1)';
                      
        saidaPSF(i,k) = saidasS_PSF(i,k) + saidasFe(i,k);
        
        
    end    
    
    
    %% -- Controlador GPC 
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),Ny(i),1)];
    end 
    
    %%% calculo da resposta livre
    du1 = [];
    for i=1:m
        du1 = [du1;du(i,k-1:-1:k-nH(i))'];
    end
    
    y1 = [];
    for i=1:n
        y1 = [y1;saidaPSF(i,k:-1:k-nF(i)+1)'];
    end
    
    f = F*y1 + H*du1;
    
    %%% montagem das matrizes de restrição
    boti = repelem(umax'-entradas(:,k-1),Nu');
    boti = [boti;
            repelem(-umin'+entradas(:,k-1),Nu')];
        
    
    foti = foti1*(R-f);
    if(usarRest)
        X = quadprog(Hoti,foti,Aoti,boti,[],[],LB,UB);
    else
        X = quadprog(Hoti,foti,[],[],[],[],[],[]);    
    end
    for i=1:m
        du(i,k) = X(sum(Nu(1:i-1))+1,1);
    end

    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasGPC1 = saidas;
entradasGPC1 = entradas;
duGPC1 = du;




%% GPC DTC com filtro 2
%% PARAMETROS DE AJUSTE DO CONTROLADOR
Fe(1,1) = z/z;
Fe(2,2) = z/z;
Fe(3,3) = z/z;

zf = 0.85;
zo1 = 0.9802;
A1 = [1 ,-1;
      zo1^2,-zo1]
b1 = [(1-zf)^2;
      zo1^dmin(1)*(zo1-zf)^2]
X1 = A1\b1;
Fe(1,1) = (X1(1)*z^2-X1(2)*z)/(z-zf)^2;

zo2 = 0.9802;
A2 = [1 ,-1;
      zo2^2,-zo2]
b2 = [(1-zf)^2;
      zo2^dmin(2)*(zo2-zf)^2]
X2 = A2\b2;
Fe(2,2) = (X2(1)*z^2-X2(2)*z)/(z-zf)^2;


nfe1 = size(Fe(1,1).den{1},2)-1;
nfe2 = size(Fe(2,2).den{1},2)-1;
nfe3 = size(Fe(3,3).den{1},2)-1;

Spsf = (eye(n)-Fe*Lz)*Gnz;

%% inicializacao dos estados e variaveis da simulação
% - obtencao dos numerados e denominadores
numS = cell(n,m);
denS = cell(n,m);

for i=1:n
	% obtencao das respostas ao degrau das entradas
	for j=1:m
		% obtencao dos numeradores e denominadores para o simulador
        [numS{i}{j}, denS{i}{j}] = tfdata(Spsf(i,j),'v');
    end
    
end


estadosEntr = cell(n,m);
estadosS_PSF = cell(n,m);
for i=1:n
    for j=1:m
		estadosEntr{i}{j} = zeros(1,nit);
        estadosS_PSF{i}{j} = zeros(1,nit);

    end
end

saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

saidasS_PSF = zeros(n,nit);
saidasFe = zeros(n,nit);
saidaPSF = zeros(n,nit);

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numR{i}{j},2);
            sizeD = size(denR{i}{j},2);

            str1 = numR{i}{j}*(entradas(j,k-atrasoEntradaR(i,j):-1:k+1-sizeN-atrasoEntradaR(i,j))'...
                               +perts(j,k-atrasoEntradaR(i,j):-1:k+1-sizeN-atrasoEntradaR(i,j))');
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end
        
        saidas(i,k) = aux ;%;+ 0.01*randn(1,1);
    end
    %% PSF
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numS{i}{j},2);
            sizeD = size(denS{i}{j},2);

            str1 = numS{i}{j}*(entradas(j,k:-1:k-sizeN+1)');
            str2 = -denS{i}{j}(1,2:end)*estadosS_PSF{i}{j}(k-1:-1:k-sizeD+1)';
            estadosS_PSF{i}{j}(k) = str1+str2;

            aux = aux + estadosS_PSF{i}{j}(k);
        end
        saidasS_PSF(i,k) = aux ;
        

        
        sizeN = size(Fe(i,i).num{1},2);
        
        saidasFe(i,k) = -Fe(i,i).den{1}(2:end)*saidasFe(i,k-1:-1:k-sizeN+1)' ...
                          +Fe(i,i).num{1}*saidas(i,k:-1:k-sizeN+1)';
                      
        saidaPSF(i,k) = saidasS_PSF(i,k) + saidasFe(i,k);
        
        
    end    
    
    
    %% -- Controlador GPC 
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),Ny(i),1)];
    end 
    
    %%% calculo da resposta livre
    du1 = [];
    for i=1:m
        du1 = [du1;du(i,k-1:-1:k-nH(i))'];
    end
    
    y1 = [];
    for i=1:n
        y1 = [y1;saidaPSF(i,k:-1:k-nF(i)+1)'];
    end
    
    f = F*y1 + H*du1;
    
    %%% montagem das matrizes de restrição
    boti = repelem(umax'-entradas(:,k-1),Nu');
    boti = [boti;
            repelem(-umin'+entradas(:,k-1),Nu')];
        
    
    foti = foti1*(R-f);
    if(usarRest)
        X = quadprog(Hoti,foti,Aoti,boti,[],[],LB,UB);
    else
        X = quadprog(Hoti,foti,[],[],[],[],[],[]);    
    end   
    for i=1:m
        du(i,k) = X(sum(Nu(1:i-1))+1,1);
    end

    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasGPC2 = saidas;
entradasGPC2 = entradas;
duGPC2 = du;



%% GPC Normal
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

N1 = dmin'+N1; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = dmin'+N2; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar

%% Controlador GPC - MIMO

%%% obtenção da representação MFD
[Amfd,Bmfd] = MFDredux(Gz.den,Gz.num);

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
            H{i,j} = [H{i,j} zeros(Ny(i),nH(j)-nHi)];
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
Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle

%%% matrizes do problema de otimização
Hoti = 2*(G'*Qy*G+Qu);
foti1 = -2*G'*Qy;

LB = [];
UB = [];
Aoti = [];
for i=1:m
    LB = [LB;repmat(dumin(i),Nu(i),1)];
    UB = [UB;repmat(dumax(i),Nu(i),1)];
    Aoti = blkdiag(Aoti,tril(ones(Nu(i))));
end
Aoti = [Aoti;-Aoti];


%% inicializacao dos estados e variaveis da simulação
estadosEntr = cell(n,m);
for i=1:n
    for j=1:m
		estadosEntr{i}{j} = zeros(1,nit);
    end
end


saidas = zeros(n,nit); % vetor da saída
entradas = zeros(m,nit); % vetor do sinal de controle
du = zeros(m,nit); % vetor dos incrementos de controle

for k = nin:nit
    %% -- simulador, não alterar
    %%% Simulador MIMO que utiliza estados
    for i=1:n
        aux = 0;
        for j=1:m
            sizeN = size(numR{i}{j},2);
            sizeD = size(denR{i}{j},2);

            str1 = numR{i}{j}*(entradas(j,k-atrasoEntradaR(i,j):-1:k+1-sizeN-atrasoEntradaR(i,j))'...
                               +perts(j,k-atrasoEntradaR(i,j):-1:k+1-sizeN-atrasoEntradaR(i,j))');
            str2 = -denR{i}{j}(1,2:end)*estadosEntr{i}{j}(k-1:-1:k-sizeD+1)';
            estadosEntr{i}{j}(k) = str1+str2;

            aux = aux + estadosEntr{i}{j}(k);
        end
        
        saidas(i,k) = aux ;%;+ 0.01*randn(1,1);
    end
    %% -- Controlador GPC 
    %%% referencias
    R = [];
    for i=1:n
        R = [R;
             repelem(refs(i,k),Ny(i),1)];
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
    boti = repelem(umax'-entradas(:,k-1),Nu');
    boti = [boti;
            repelem(-umin'+entradas(:,k-1),Nu')];
        
    
    foti = foti1*(R-f);
    if(usarRest)
        X = quadprog(Hoti,foti,Aoti,boti,[],[],LB,UB);
    else
        X = quadprog(Hoti,foti,[],[],[],[],[],[]);    
    end   
    for i=1:m
        du(i,k) = X(sum(Nu(1:i-1))+1,1);
    end

    entradas(:,k) = entradas(:,k-1)+du(:,k);

end

saidasGPC0 = saidas;
entradasGPC0 = entradas;
duGPC0 = du;



%% Gera gráficos
cores = gray(5);
cores = cores(1:end-1,:);

t = ((nin:nit)-nin)*Ts;
ind = nin:nit;


hf= figure
h=subplot(3,2,1)
plot(t,saidasGPC0(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGPC1(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC2(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))
grid on
xlim([0 450])
ylim([-0.1 1.5])
ylabel('y_1')
set(h, 'FontSize', tamletra);
title('Controladas','FontSize', tamletra)



h=subplot(3,2,3)
plot(t,saidasGPC0(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGPC1(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC2(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))
grid on
xlim([0 450])
ylim([-0.1 0.55])
ylabel('y_2')
set(h, 'FontSize', tamletra);



h=subplot(3,2,5)
plot(t,saidasGPC0(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGPC1(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC2(3,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(3,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))
grid on
xlim([0 450])
ylim([-0.1 0.4])
ylabel('y_3')
set(h, 'FontSize', tamletra);
xlabel('Tempo (min)')

hl = legend('GPC','DTC-GPC1','DTC-GPC2','Referências','Location','SouthEast')



h=subplot(3,2,2)
plot(t,entradasGPC0(1,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGPC1(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC2(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

grid on
xlim([0 450])
% ylim([-0.1 0.9])
ylabel('u_1')
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);
title('Manipuladas','FontSize', tamletra)



h=subplot(3,2,4)
plot(t,entradasGPC0(2,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGPC1(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC2(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

grid on
xlim([0 450])
% ylim([-0.75 0.1])
ylabel('u_2')
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);



h=subplot(3,2,6)
plot(t,entradasGPC0(3,ind),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,entradasGPC1(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,entradasGPC2(3,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

grid on
% xlim([0 220])
% ylim([-0.25 0.25])
ylabel('u_3')
set(h,'YAxisLocation','right')
set(h, 'FontSize', tamletra);
xlabel('Tempo (min)')




% hf.Position = tamfigura;
hl.Position = [0.3332 0.2913 0.2179 0.1750];


% if(usarRest)
%     print('fracionador_gpc_dtc_final_erro_rest','-depsc')
% else
%     print('fracionador_gpc_dtc_final_erro','-depsc')
%     
% end

%%
hf= figure
h=subplot(3,1,1)
plot(t(1:200),saidasGPC0(1,ind(1:200)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGPC1(1,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC2(1,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(1,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))
grid on
xlim([0 450])
ylim([-0.1 1.5])
ylabel('y_1')
set(h, 'FontSize', tamletra);
title('Controladas','FontSize', tamletra)



h=subplot(3,1,2)
plot(t(1:200),saidasGPC0(2,ind(1:200)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGPC1(2,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC2(2,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(2,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))
grid on
xlim([0 450])
ylim([-0.1 0.6])
ylabel('y_2')
set(h, 'FontSize', tamletra);



h=subplot(3,1,3)
plot(t(1:200),saidasGPC0(3,ind(1:200)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,saidasGPC1(3,ind),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,saidasGPC2(3,ind),':','LineWidth',tamlinha,'Color',cores(3,:))

plot(t,refs(3,ind),'-.','LineWidth',tamlinha,'Color',cores(4,:))
grid on
xlim([0 450])
ylim([-0.1 0.4])
ylabel('y_3')
set(h, 'FontSize', tamletra);
xlabel('Tempo (min)')

hl = legend('GPC','DTC-GPC1','DTC-GPC2','Referências','Location','SouthEast')




% hf.Position = tamfigura;
hl.Position = [0.7636 0.3389 0.2179 0.1750];


% if(usarRest)
%     print('fracionador_gpc_dtc_final_erro2_rest','-depsc')
% else
%     print('fracionador_gpc_dtc_final_erro2','-depsc')
%     
% end