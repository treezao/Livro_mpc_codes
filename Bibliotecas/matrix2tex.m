function txt = matrix2tex(G,varargin)
    
    if(nargin>1)
        N = varargin{1};
    else
        N = 4;
    end

    txt = [newline];
    nl = size(G,1);
    nc = size(G,2);
    
    
    for i=1:nl
        for j=1:nc-1
            txt = [txt,num2str(G(i,j),N),'  &  '];
        end
        
        if(i<nl)
            txt = [txt,num2str(G(i,end),N), '\\',newline];
        else
            txt = [txt,num2str(G(i,end),N), newline];
        end
    end
    
    txt = strrep(txt,'.',',');
end