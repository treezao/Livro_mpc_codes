%% simulação das equações não lineares do sistema CSTR
function [out] = simModelo(Co0,Cf,Ts)

global Qf Vc k1 k2 


dydt = @(t,y) (Qf/Vc)*(Cf-y) - k1*y/(k2*y+1)^2;


[tout,yout] = ode45(dydt,[0 Ts],Co0);

out = yout(end);
end


