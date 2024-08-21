function [hang1,hang2] = motaishu(g,mshu)
    hang1=1;
    hang2=0;
    for shu=1:g
        if shu==1
            hang2=hang2+mshu{shu};
        else
            hang1=hang1+mshu{shu-1};     
            hang2=hang2+mshu{shu};
        end
    end
end


