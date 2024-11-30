function res = res_deal(Ibin, fgmask)

    res = zeros(size(Ibin));
    L = bwlabel(Ibin);
    stats = regionprops(Ibin); 
    Ar = cat(1, stats.Area);  
    ind = find(Ar ==max(Ar));
    res(find(L==ind)) = 255;
    res = uint8(res);
   
    for i=1:length(Ar)
        mid_res = fgmask .* res;
        if ~any(mid_res(:))
        res = zeros(size(Ibin));
        sorted_Ar = sort(Ar);
        ind = find(Ar == sorted_Ar(end-1));
        res(find(L==ind)) = 255;
        else
            break;
        end
    end
end