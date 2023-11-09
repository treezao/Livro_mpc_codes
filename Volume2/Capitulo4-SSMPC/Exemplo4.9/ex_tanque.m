clear all, 
close all, 
clc


%% modelo do sistema
%%% parâmetros do modelo do tanque
a1 = 0.1;
a2 = 0.4;
a3 = 0.5;
a4 = a3;

%%% matrizes do modelo
A = [(1-a2), 0;
     a2, 1]
B = diag([a1,a3])
Bq = [0;-a4]
C = eye(2)
Cq = zeros(2,1);


%% obtenção das matrizes do controlador SSMPC para a sintonia 1

G = [C*B zeros(2,2);
     C*(eye(2)+A)*B, C*B;
     C*(eye(2)+A+A^2)*B,C*B]
 
 
F = [C*(A), -C*(A);
     C*(A+A^2), -C*(A+A^2);
     C*(A+A^2+A^3), -C*(A+A^2+A^3)]
 


Qei = diag([0.8,1]);
Qui = diag([1,0.6]);

Qe = blkdiag(Qei,Qei,Qei);
Qu = blkdiag(Qui,Qui);

K = (G'*Qe*G+Qu)\G'*Qe

Gq = [Cq zeros(2,1) zeros(2,1);
      C*Bq+Cq, Cq, zeros(2,1);
      C*(eye(2)+A)*Bq+Cq,C*Bq+Bq,Cq]


Fq = [C*Bq+Cq,-C*Bq;
      C*(eye(2)+A)*Bq+Cq,-C*(eye(2)+A)*Bq;
      C*(eye(2)+A+A^2)*Bq+Cq,-C*(eye(2)+A+A^2)*Bq]

 
  
%%  obtenção do controlador equivalente e da malha fechada
%%% controlador equivalente
K1 = K(1:2,:)
KF = K(1:2,:)*F
KF0 = KF(:,1:2)
KF1 = KF(:,3:4)



KFq = K1*Fq
KFq0 = KFq(:,1)
KFq1 = KFq(:,2)

Kr = K1*repmat(eye(2),[3,1])

Ky = Kr*C

Ac = [eye(2) -KF1,-KFq1;
      zeros(2,5);
      zeros(1,5)]
Bc = [Kr, -KF0-Ky,-KFq0;
      zeros(2,2),eye(2,2),zeros(2,1);
      zeros(1,2+2),eye(1)]
Cc = [eye(2),-KF1,-KFq1]
Dc = [Kr,-KF0-Ky,-KFq0]


%%% malha fechada

Amf = [A-B*(KF0+Kr*C),B*Cc;
       [-KF0-Kr*C;eye(2);zeros(1,2)], Ac]
Bmf = [B*Kr,(Bq-B*KFq0);
       Kr,-KFq0;
       zeros(2,3);
       zeros(1,2),eye(1)]
Cmf = [C, zeros(2,5);
       -KF0-Kr,Cc;
       -C,zeros(2,5)]
Dmf = [zeros(2,2) Cq;
       Kr, -KFq0;
       eye(2),-Cq]
   
sysmf1 = ss(Amf,Bmf,Cmf,Dmf,1) % sistema completo em espaço de estados
