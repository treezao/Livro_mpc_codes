%% função de saturação
%%% satura os elementos de e entre os valores umin e umax
function [u] = sat(e,umin,umax)

u = e;

for i=1:size(e,1)
    for j=1:size(e,2)
        if(u(i,j)>umax)
            u(i,j)=umax;
        elseif(u(i,j)<umin)
            u(i,j)=umin;
        end
    end
end
