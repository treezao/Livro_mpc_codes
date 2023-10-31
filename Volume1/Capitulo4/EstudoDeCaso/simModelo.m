function out = simModelo(ae,as,hant)
%%% parametros do modelo
a1 = 2;
a2 = 0.2;
A = 10;
g = 9.8;

%%% saturação do sinal de controle
if(ae>1)
    ae=1;
elseif(ae<0)
    ae = 0;
end


out = hant +(a1/A)*ae - (a2/A)*as*sqrt(2*g*hant);

%%% saturação do nível do tanque
if(out>5)
    out = 5;
elseif(out<0)
    out = 0;
end

end


