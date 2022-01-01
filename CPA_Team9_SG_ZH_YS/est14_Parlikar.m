function [uncalib_co_values, full_taus] = est14_Parlikar(MAP,Period,Pdias,PP)
    

    alpha = 2;
    Period = Period/(60*125); 
%     del_Pi = diff(PP);
%     del_P = [PP(1); del_Pi];
    del_P = PP;
    beta = zeros(size(MAP,1),1);
    window = 13;
    full_taus = zeros(size(MAP,1),1);
    den_y = []; num_y = [];

    for i = 1: size(MAP,1)
        if i < ((window + 1)/2) || i >= size(MAP,1) - ((window + 1)/2)
            continue
        else

            low = (i - ((window-1)/2));

            u = (i+((window-1)/2));

            x = MAP(low:u);
            den_y = Period(low:u);
            num_y = alpha*(MAP(low:u) - Pdias(low:u)) - del_P(low:u);
            y = num_y./den_y;

            full_taus(i) = 1/(sum(x.*y)./sum(x.^2));
        end
    end


    nearestfun = @(taus) interp1(find(taus),taus(taus~=0),(1:length(taus))','nearest','extrap');
    full_taus = 0.5*(nearestfun(full_taus) + flip(nearestfun(flip(full_taus))));

    uncalib_co_values = ((del_P./Period) + (MAP./full_taus));