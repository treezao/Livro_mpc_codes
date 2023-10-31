%%% cálculo dos filtros dos erros de predição (SISO) para o GDMC

F = tf(0,1,Ts);
nf = 2; %ordem do filtro

pz = pole(Gz); % obtem o polo indesejado (instável) de malha aberta

for i=N1(1):N2(1)
    %%% monta sistema para obtenção dos parâmetros do filtro
    %%% primeira equação (z^i - F_i(z) = 0 -> para z=pz
    %%% segunda equação F_i(1) = 0 -> força ganho unitário
    indf = i-N1(1)+1;
    Af = [pz^2 pz;
          1 1];
    bf = [pz^i*(pz-betaf)^2;
          (1-betaf)^2];
    X = Af\bf;
    F(indf,1) = (X(1)*z^2+X(2)*z)/(z-betaf)^2;
    %%% teste da condição
%     pz^i-(X(1)*pz+X(2))/(pz-betaf)^2

    %%% armazena coeficientes gtil
    modDegrauUF{i} = filter(F(indf,1).num{1},F(indf,1).den{1},Gcoef);

end


%%% calcula a matriz H para o cálculo da resposta livre no caso GDMC
H1 = [];
H2 = [];

for i=N1(1):N2(1)
    H1 = [H1;Gcoef(i+1:i+Nf)'];
    H2 = [H2;modDegrauUF{i}(1:Nf)'];
    
end
H = H1-H2


%% inicialização vetores
entradas = Cfbar*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = Cobar*ones(nit,1); % vetor com as saídas do sistema

erro = zeros(nit,1); % vetor de erros

yfilt = Cobar*ones(nit,N); % saídas filtradas

%% simulação com referencia futura
for i = nin:nit
    %% simulação do modelo processo
    saidas(i) = simModelo(saidas(i-1),entradas(i-1)+perts(i-1),Ts);
    
    erro(i) = refs(i)-saidas(i);
    %% Controlador
    %%% referencias
    R = refs(i)*ones(N,1);
    
    %%% calculo da resposta livre
    for k=1:N(1)
        yfilt(i,k) = -F(k,1).den{1}(2:end)*yfilt(i-1:-1:i-nf,k) + F(k,1).num{1}*saidas(i:-1:i-nf);
    end
    
    f = H*du(i-1:-1:i-Nf) + yfilt(i,:)';
    
    %%% calculo do incremento de controle ótimo    
    % com restrições
    fqp = fqp1*(R-f);
    
    rbar = [repelem(umax-entradas(i-1),Nu)';
             repelem(entradas(i-1)-umin,Nu)';
             ymax-f;
             -ymin+f];
    [X,FVAL,EXITFLAG] = quadprog(Hqp,fqp,Rbar,rbar,[],[],LB,UB);
    
    %%% caso dê infactível, remover a restrição na saída
    if(EXITFLAG==-2)
        rbar = [repelem(umax-entradas(i-1),Nu)';
                 repelem(entradas(i-1)-umin,Nu)'];
        X = quadprog(Hqp,fqp,Rbar(1:2*Nu,:),rbar,[],[],LB,UB);
    end

    du(i) = X(1);
    
    
    % sem restrições
%     du(i) = Kdmc1*(R-f);
    
    %%% aplicar no sistema
    entradas(i) = entradas(i-1)+du(i);    
    
end

