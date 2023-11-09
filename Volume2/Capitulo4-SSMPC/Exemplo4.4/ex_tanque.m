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
 
Gq = [Cq zeros(2,1) zeros(2,1);
      C*Bq+Cq, Cq, zeros(2,1);
      C*(eye(2)+A)*Bq+Cq,C*Bq+Bq,Cq]


Fq = [C*Bq+Cq,-C*Bq;
      C*(eye(2)+A)*Bq+Cq,-C*(eye(2)+A)*Bq;
      C*(eye(2)+A+A^2)*Bq+Cq,-C*(eye(2)+A+A^2)*Bq]

