close all
clear all

run('../../../Bibliotecas/parametrosFiguras.m')

%%

Ts = 1; % em minutos
G = tf(0.4,[1 -0.9],Ts);
z = tf('z',Ts); % variável z da Transformada-Z

Nss = 66; % horizonte de modelo

g = step(G,Ts:Ts:Nss*Ts); % obtendo os coeficientes da resposta ao degrau unitário
h = g - [0;g(1:end-1)] ; % obtendo os coeficientes da resposta ao impulso unitário


Gfir = tf(0,1,Ts); % montagem do modelo FIR
for i=1:Nss
    Gfir = Gfir + h(i)*z^(-i);
end


%% diagrama de bode dos modelos

w = logspace(-3, log10(2*pi/Ts/2),200); % faixa de frequencia

[mG,pG] = bode(G,w); %magnitude e fase de Gu
[mGfir,pGfir] = bode(Gfir,w); % magnitude e fase de Gfir


cores = gray(3);
cores = cores(1:end-1,:);

hf = figure
h = subplot(2,1,1)
semilogx(w,20*log10(mG(:)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w,20*log10(mGfir(:)),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx([2*pi/Ts/2, 2*pi/Ts/2],[-100 100],'k') 

ylim([-15 15])
ylabel('Magnitude (dB)','FontSize', tamletra)
hl = legend('G_u(z)','G_{FIR}(z)','Location','NorthEast','FontSize', tamletra)
grid on

h2 = subplot(2,1,2)
semilogx(w,(pG(:)),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w,(pGfir(:)),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx([2*pi/Ts/2, 2*pi/Ts/2],[-200 200],'k') 

ylim([-180 0])
ylabel('Fase (graus)','FontSize', tamletra)
xlabel('Frequência (rad/min)','FontSize', tamletra)
grid on

hf.Position = tamfigura;
% hf.Position = [0.1527 0.1769 0.2054 0.1806];

set(h, 'FontSize', tamletra);
set(h2, 'FontSize', tamletra);