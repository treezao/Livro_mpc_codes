clc
clear all
close all
%% DEFINIÇÃO DOS MODELOS
N_entradas = 1; % número de entradas da planta - manipuladas
N_saidas = 1; %número de saidas da planta - controladas
N_perturbacoes =0; %número de perturbações mensuráveis
usar_feed_forward = 0;
flag_restricao_du = 1;
flag_restricao_u = 0;
flag_restricao_y = 1;
INTERVALHO_RT_FALHA = 3; %Intervalo de falta de tempo para a oimização 
s = tf('s');
% Definir os modelos que relacionam {saida x,entrada y};
Ts = 0.1;
sys = 4/(s*(s+1)^2);
sysd{1,1} = c2d(sys,Ts,'zoh');
T = [1 -1];
A{1,1} = sysd{1,1}.den{1};
Atio{1,1} = conv(A{1,1},T)
B{1,1} = sysd{1,1}.num{1}(2:end);
B_p{1,1} = [0.5];
dimB = size(B{1,1},2);
dimB_p = size(B_p{1,1},2);
dimA = size(A{1,1},2);

%% PARAMETROS DE AJUSTE DO CONTROLADOR

Nu = [25]; %horizontes de controle - dimensão 1 x N_entradas
N1 = [1]; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x N_saidas
N2 = [90]; % fim dos horizontes de predição - dimensão 1 x N_saidas
N = N2-N1+1; %horizonte de predição - não editar
solver = 'gpad'; %escolha do solver

for i = 1:N_saidas
    for j = 1:N_entradas
        g{i,j} = step(sysd{i,j},N2(i));%Ts:Ts:N(i)*Ts);
    end
end
           
pondY = [1]; %ponderação nos erros - dimensão 1 x N_saidas
pondU = [1]; %ponderação nas ações de controle - dimensão 1 x N_entradas       

Qy = diag(repelem(pondY,1,N)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(pondU,1,Nu)); %montagem da matriz de ponderação da ação de controle

%Matrizes de restrições 

%% PARAMETROS DA SIMULAÇÃO
Tini = dimA+1+5;
Tsim = 60; %tempo de simulação
degraus_ref{1} = [0 1]; %Amplitudes dos degraus no ensaio para cada saída
degraus_pert{1} = [0 0.05];
degraus_t{1} = [2 Tsim-2];
degraus_pert_t{1} = [20 Tsim-20];

Tsim = Tsim + Tini;
for i = 1:N_saidas
    ref{i} = repelem(degraus_ref{i},1,[Tini N(i)+N1(i)]+degraus_t{i});%vetor de referências - [valores dos degraus] - [repetição do respectivo valor]  
end 
for i = 1:N_perturbacoes
    pert{i} = repelem(degraus_pert{i},1,[Tini N(i)+N1(i)]+degraus_pert_t{i});%vetor de perturbações - [valores dos degraus] - [repetição do respectivo valor]
end

 %% RODA SIMULAÇÃO
 %inicializa variaveis
max_iter = 0;
for i = 1:4
    for j = 1:Tini-1
        y{i}(j) = 0;
        du{i}(j) = 0;
        u{i}(j) = 0;
    end
end
fma{1} = zeros(N2(1)+dimA,1);
ysim{1} = zeros(N2(1),1);


%% Monta as matrizes offline que definem o problema MIMO DMC
for i=1:Nu(1)
    Gm{1,1}(i:N(1),i) = g{1,1}(N1(1):N2(1)-i+1);
end
G = cell2mat(Gm);
H = 2*(G'*Qy*G+Qu);
invH = inv(H);

R = [];
r = [];
%restrição no incremeto de controle
if flag_restricao_du
    du_lim = [10.05 -10.05]; %limites por ação de controle ex:2 [sup_1 sup_2 inf_1 inf_2]
    R = [eye(sum(Nu));-eye(sum(Nu))];
    R_inc = [eye(sum(Nu));-eye(sum(Nu))];
    r_inc = [repelem(du_lim(1:N_entradas),1,Nu)';repelem(-1*ones(1,N_entradas).*du_lim(N_entradas+1:end),1,Nu)'];
end

if flag_restricao_u
    R_con = [tril(ones(sum(Nu)));-tril(ones(sum(Nu)))];
    R = [R; R_con];
    u_lim = [0.15 0];
end

if flag_restricao_y
%    R_y = [G;-G];
%    R = [R; R_y];
    R = [R; G];
    y_lim = [1 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:4
    controle_n = 0 
    for k=Tini:Tsim
     %%%%%%%%%roda a planta
        ysim{p} = 0;
        ysim{p} = -fliplr(Atio{1,1}(2:dimA+1))*y{p}(k-dimA:k-1)';
        ysim{p} = ysim{p} + fliplr(B{1,1}(1:dimB))*du{p}(k-dimB:k-1)';
        if N_perturbacoes
            ysim{p} = ysim{p} + -fliplr(B_p{1,1}(1:dimB_p))*(pert{1}(k-dimB_p:k-1)-pert{1}(k-dimB_p-1:k-2))';
        end
        y{p}(k) = ysim{p}(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Calcula as matrizes online do GPC
        f_aux = [];
        
        du_aux{1} = zeros(dimB+N2(1),1);
        du_aux{1}(1:dimB-1)= du{p}(k-dimB+1:k-1)';
        
        fma{1}(1:dimA) = y{p}(k-dimA+1:k)';
        for n = 1:N2
            fma{1}(dimA+n) = -fliplr(Atio{1,1}(2:dimA+1))*fma{1}(n:dimA+n-1);
            fma{1}(dimA+n) = fma{1}(dimA+n) + fliplr(B{1,1}(1:dimB))*du_aux{1}(n:dimB+n-1);
        end
        f_aux = [f_aux; fma{1}(dimA+N1(1):N2(1)+dimA)];
        f = f_aux;
    
        %refm(1:N(1)) = ref{1}(k+N1(1):k+N1(1)+N(1)-1);
        refm(1:N(1)) = ones(1,N(1))*ref{1}(k);
        b = (2*(f-refm')'*Qy*G)';

        r = [];
        if flag_restricao_du
            r = r_inc;
        end
        if flag_restricao_u
            u_lim_p = u_lim(1)*ones(sum(Nu));
            u_lim_n = u_lim(2)*ones(sum(Nu));
            r_p=[];
            r_m= [];
            for i = 1:N_entradas
                r_p = [r_p; repelem(u_lim_p(i)-u{i}(k-1),1,Nu(i))'];
                r_m = [r_m; repelem(-u_lim_n(N_entradas+i)+u{i}(k-1),1,Nu(i))'];
            end
            r_con = [r_p; r_m];
            r = [r; r_con]
        end
        
        if flag_restricao_y
            %y_lim_p = y_lim(1)*ones(sum(N),1);
            %y_lim_n = y_lim(2)*ones(sum(N),1);

 %           r_p=[];
 %           r_m= [];
 %           for i = 1:N_entradas
 %               r_p = [r_p; repelem(y_lim_p(i)-f,1,Nu(i))'];
 %               r_m = [r_m; repelem(-y_lim_n(N_entradas+i)+f,1,Nu(i))'];
  %          end
            %r_y = [y_lim_p-f; -y_lim_n+f];
            r_y_p = [refm(1:N(1))'-f];
            r_y_n = [refm(1:N(1))'-f];
%            if ref{1}(k) == y{1}(k)
 %               r_y = [];
 %           elseif ref{1}(k) > y{1}(k)
                r_y = refm(1:N(1))'-f;
 %           else
  %              r_y = -f+refm(1:N(1))';
 %           end

            r = [r; r_y];
        end

        %% Atualiza as variáveis de controle
        controle_n = controle_n+1;
        switch p
            case 1 %caso sem perda
                if ref{1}(k) == y{1}(k)
                    delta_u = quadprog(H,b,[],[]);
                else
                    delta_u = quadprog(H,b,R,r);
                end
                du{p}(k) = delta_u(1);
            case 2 % caso repete anterior
                if controle_n == INTERVALHO_RT_FALHA 
                    du{p}(k) = 0;
                    controle_n = 0;
                else
                    [delta_u,iter] = gpad_f(invH,b,R,r,1000)
                    du{p}(k) = delta_u(1);
                end

            case 3 %caso subótimo
                if controle_n == INTERVALHO_RT_FALHA 
                    [delta_u,iter] = gpad_f(invH,b,R,r,10);
                    controle_n = 0;
                else
                    [delta_u,iter] = gpad_f(invH,b,R,r,1000);
                end
                du{p}(k) = delta_u(1);
                
            case 4 %caso saturação
                if controle_n == INTERVALHO_RT_FALHA 
                    %[delta_u,iter] = gpad_f(invH,b,R_inc,1000*ones(2*Nu(1),1),1000);
                    % if delta_u(1) > du_lim(1)
                    %     delta_u(1) = du_lim(1);
                    % end
                    % if delta_u(1) < du_lim(2)
                    %     delta_u(1) = du_lim(2);
                    % end
                    delta_u = quadprog(H,b,[],[]);
                     controle_n = 0;
                else
                    [delta_u,iter] = gpad_f(invH,b,R,r,1000);
                end
                du{p}(k) = delta_u(1);
        end
        u{p}(k) = u{p}(k-1) + du{p}(k);
        if flag_restricao_u
            if (u{p}(k) >u_lim(1)) && (p==4)  
                u{p}(k) = u_lim(1);
                du{p}(k) = u{p}(k)- u{p}(k-1); 
            end
        end
    end
end
    
%%PLOTA OS GRÁFICOS
% parametros

t = Ts*(Tini:Tsim);
tamletra=11; %tamanho da letra
tamlinha = 1.5; % tamanho da linha
tamfigura = [200 200 560 310];
cores = gray(6);
cores = cores(1:end-1,:);

%tamnum=20;
fh = figure();
subplot(2,1,1);

%plot(t,ref{1}(Tini:Tsim),':',t,y{1}(Tini:Tsim),'-',t,y{2}(Tini:Tsim),'--',t,y{3}(Tini:Tsim),'-.',t,y{4}(Tini:Tsim),':*','LineWidth',2,'color','k')
plot(t,ref{1}(Tini:Tsim),'-','LineWidth',tamlinha,'Color',cores(5,:))
hold on
plot(t,y{1}(Tini:Tsim),'-','LineWidth',tamlinha,'Color',cores(1,:))
plot(t,y{2}(Tini:Tsim),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,y{3}(Tini:Tsim),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,y{4}(Tini:Tsim),':','LineWidth',tamlinha,'Color',cores(4,:))
legend('Referência','GPC','GPC 1','GPC 2','GPC 3','Location','SouthEast')
%axis([Tini*Ts (Tsim-Tini)*Ts -0.3 1.2]);
 %   ylabel('Saídas e Referências do Sistemas', 'FontSize', tamletra);
xlabel('Tempo (s)', 'FontSize', tamletra);
ylabel('Controladas', 'FontSize', tamletra);
set(gca,'FontSize',tamletra);
grid on
%cd joao
%print -depsc results_con_ex_49_1.eps

%figure
subplot(2,1,2);
%plot(t,u{1}(Tini:Tsim),'-',t,u{2}(Tini:Tsim),'--',t,u{3}(Tini:Tsim),'-.',t,u{4}(Tini:Tsim),':*','LineWidth',2,'color','k')
plot(t,u{1}(Tini:Tsim),'-','LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,u{2}(Tini:Tsim),'--','LineWidth',tamlinha,'Color',cores(2,:))
plot(t,u{3}(Tini:Tsim),'-.','LineWidth',tamlinha,'Color',cores(3,:))
plot(t,u{4}(Tini:Tsim),':','LineWidth',tamlinha,'Color',cores(4,:))
%title(['y1' ' e Referência do Sistema'], 'FontSize', tamletra)
%legend('U_1','U_2','U_3','U_4','Location','NorthEast')
%axis([Tini*Ts (Tsim-Tini)*Ts 0 0.3]);
 %   ylabel('Saídas e Referências do Sistemas', 'FontSize', tamletra);
xlabel('Tempo (s)', 'FontSize', tamletra);
ylabel('Manipuladas', 'FontSize', tamletra);
set(gca,'FontSize',tamletra);
grid on
%fh.WindowState = 'maximized';

fh.Position = tamfigura;

print -depsc tempo_insuficiente2.eps