function out = trocaponto(h)

for i=1:size(h,1)
    out{i} = strrep(h{i},'.',',');
end


end
