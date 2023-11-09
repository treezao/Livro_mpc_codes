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


%% obtenção das matrizes do controlador SSMPC para a sintonia 2
 
G2 = [C*B zeros(2,1);
     C*(eye(2)+A)*B, C*B(:,1);
     C(2,:)*(eye(2)+A+A^2)*B,C(2,:)*B(:,1)]
 


 
F2 = [C*(A), -C*(A);
     C*(A+A^2), -C*(A+A^2);
     C(2,:)*(A+A^2+A^3), -C(2,:)*(A+A^2+A^3)]
  

Qe2 = blkdiag(Qei,Qei,1);
Qu2 = blkdiag(Qui,1);


K2 = (G2'*Qe2*G2+Qu2)\G2'*Qe2
