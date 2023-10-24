clc, 
clear all, 
close all
run('../../../Bibliotecas/parametrosFiguras.m')


%% Plot da função custo e seus valores para diferentes vetores solução

nit = 200; % número de pontos
t = linspace(0,2,nit); % vetor de pontos

[u1,u2] = meshgrid(t,t); % combinação de pontos 

n = size(u1,1);
m = size(u1,2);

%%% definição da função custo e cálculo de seu valor para diferentes pontos
x0 = 1.5;
y0 = 0.5;
xy0 = [x0;y0];

xa = 1;
yb = .25;

teta = deg2rad(30);
R = [cos(teta), -sin(teta);
     sin(teta),cos(teta)];

H = [ 1/xa^2 0; 
     0 1/yb^2];
H = inv(R)'*H*inv(R);

b = (-2*H'*xy0);
c = xy0'*xy0+3;

J = zeros(n,m);
imin = 1;
jmin = 1;
Jmin = Inf;
for i=1:n
    for j=1:m
        v = [u1(i,j);u2(i,j)];
        J(i,j) = v'*H*v+b'*v+c;
        
        if(J(i,j)<Jmin)
            Jmin = J(i,j);
            imin = i;
            jmin = j;
        end
        
    end
end

%% Soluções

%%% solução ótima sem restrições (J1)
xo = inv(2*H)*-b

%%% solução saturada manualmente (J3)
xsat = xo;
if(xsat(1)<0)
    xsat(1)=0;
elseif(xsat(1)>1)
    xsat(1)=1;
end
if(xsat(2)<0)
    xsat(2)=0;
elseif(xsat(2)>1)
    xsat(2)=1;
end

xsat
Xsatval = xsat'*H*xsat+b'*xsat+c

%%% solução ótima com restrições (J2)
[X,Jval]=quadprog(2*H,b,[],[],[],[],[0;0],[1;1]);
X
Xval = Jval+c

%%
figure
hs1 = subplot(1,2,1)
deci=8;
hs = mesh(u1(1:deci:end,1:deci:end),u2(1:deci:end,1:deci:end),J(1:deci:end,1:deci:end))
colormap('gray')
zlim([0,5])

xlabel('\Delta u(k)')
hx1 = ylabel('\Delta u(k+1)')
hx1.Position = [-0.3084 1.0763 -0.0149];
zlabel('J')

hold on
hold on
plot3([1 1],[0 1],[0 0],'--','LineWidth',2,'Color',[0.50,0.54,0.53])
plot3([0 1],[1 1],[0 0],'--','LineWidth',2,'Color',[0.50,0.54,0.53])

hold on
plot3(xo(1),xo(2),0,'rx','LineWidth',5)
plot3(xsat(1),xsat(2),0,'rx','LineWidth',5)
plot3(X(1),X(2),0,'rx','LineWidth',5)

t10 = text(xo(1),xo(2),0,'J_1','FontSize',15)
t20 = text(X(1),X(2),0,'J_2','FontSize',15)
t30 = text(xsat(1),xsat(2),0,'J_3','FontSize',15)

t10.Position =  [1.5910 0.4883 0.0746];
t20.Position =  [1.1253 0.2001 -0.0187];
t30.Position =  [0.9871 0.5548 0.3357];

hs1.FontSize = 12;

h = gca
h.CameraPosition = [-3.0098 -14.7312 17.5939];

hs2 = subplot(1,2,2)

linhas = [0, ...
           logspace(log10(min(min(J))),log10(Xval),3), ...
           mean([Xval,Xsatval])*0.95,...
          logspace(log10(Xsatval),log10(max(max(J))),10)];
contour(J,sort(linhas),'LineWidth',2)
colormap('gray');




h = gca;
lt = round(linspace(0,2,5),1);
lt2 = linspace(0,nit,5);
for i=1:5
    h.XTick(i) = lt2(i); 
    h.XTickLabel{i} = num2str(lt(i),'%0.1f');
end

lt = round(linspace(0,2,5),1);
lt2 = linspace(0,nit,5);

h.YTick = [];
h.YTickLabel = {};
for i=1:5
    h.YTick(i) = lt2(i); 
    h.YTickLabel{i} = num2str(lt(i),'%0.1f');
end



fact = 1.5;

temp = h.Colormap;
for i=1:length(h.Colormap)
    temp(i) = temp(i)*fact;
    if(temp(i)>1)
        temp(i)=1;
    end
    
end
    
hold on
plot([1 1]*nit/2,[0 1]*nit/2,'--','LineWidth',2,'Color',[0.50,0.54,0.53])
plot([0 1]*nit/2,[1 1]*nit/2,'--','LineWidth',2,'Color',[0.50,0.54,0.53])


hold on
plot(xo(1)/2*nit,xo(2)/2*nit,'rx','LineWidth',5)
plot(xsat(1)/2*nit,xsat(2)/2*nit,'rx','LineWidth',5)
plot(X(1)/2*nit,X(2)/2*nit,'rx','LineWidth',5)

t1 = text(xo(1)/2*nit,xo(2)/2*nit,'J_1','FontSize',15)
t2 = text(X(1)/2*nit,X(2)/2*nit,'J_2','FontSize',15)
t3 = text(xsat(1)/2*nit,xsat(2)/2*nit,'J_3','FontSize',15)

t1.Position =  [148.7850 57.5840 0];
t2.Position =  [75.6384 26.9409 0];
t3.Position =  [74.6472 56.8713 0];

xlabel('\Delta u(k)')
ylabel('\Delta u(k+1)')
title('Projeção de J')

hs2.FontSize = 12;

ha = gca;
for i=1:size(ha.YTickLabel,1)
    ha.YTickLabel{i} = strrep(ha.YTickLabel{i},'.',',');
end

for i=1:size(ha.XTickLabel,1)
    ha.XTickLabel{i} = strrep(ha.XTickLabel{i},'.',',');
end

hx1 = annotation('textarrow', [0.7504 0.6361], [0.5538 0.8376])
hx2 = annotation('textbox', [0.6500 0.7490 0.2029 0.0695],'String','Aumento de J','FitBoxToText', 'on','Fontsize',tamletra,'LineStyle','none')

