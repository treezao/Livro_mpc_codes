function u = gdmc_simulink(ent)
% r -> referência
% yfilt -> saídas filtradas
% uant -> sinal de controle passado

global Kdmc1 Nf H N
persistent dup uant

if(isempty(dup))
    dup = zeros(Nf,1);
end

if(isempty(uant))
    uant = 0;
end

r = ent(1);
yfilt = ent(2:end);


%%% resposta livre
f = H*dup + yfilt;
    
%% Resolve o problema de otimização
du = Kdmc1*(r*ones(N,1) - f);
u = uant + du;

%% atualização de variáveis
uant = u;
dup = [du; dup(1:end-1)];