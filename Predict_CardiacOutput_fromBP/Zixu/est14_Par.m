function x=est14_Par(MAP,Period,Pdias,PP,F1)
%alpha = 2;
%MAP = F1(:,6);
%Period = F1(:,7);
alpha = 2;
%Pdias = F1(:,4);
del_Pi = diff(PP);
del_P = [del_Pi; del_Pi(end)];

% num_tau = MAP.*Tn;
% den_tau = alpha*(MAP - DAP) - del_P;
%% 

window = 21;
%taus = zeros(size(r_feat,1),1);%initialize
% num_y = zeros(size(r_feat,1),1);
% den_y = zeros(size(r_feat,1),1);
c_end=size(F1,1)-(window-1)/2;
c_1=1+(window-1)/2;
no_window=c_end-c_1+1;
m=size(F1,1);

%the windows except for the last window
    for i=1:no_window-1
    c(i)=i+(window-1)/2;
    for j=1:window
        for a=0:window-1
        x(j)=MAP(c(i)-(window-1)/2+a);
        den_y(j) = Period(c(i)-(window-1)/2+a);
        num_y(j) = alpha*(MAP(c(i)-(window-1)/2+a) - Pdias(c(i)-(window-1)/2)+a) - del_P(c(i)-(window-1)/2+a);
        y(j) = num_y(j)/den_y(j);
        end
        tau(j)=sum(x.^(2))./sum(x.*y);
    end
    taus(i)=tau(1);
    end
    
    %the last window
    for i=(length(taus)+1):m
        for j=1:window
        for a=0:window-1
        x(j)=MAP(c_end-(window-1)/2+a);
        den_y(j) = Period(c_end-(window-1)/2+a);
        num_y(j) = alpha*(MAP(c_end-(window-1)/2+a) - Pdias(c_end-(window-1)/2)+a) - del_P(c_end-(window-1)/2+a);
        y(j) = num_y(j)/den_y(j);
        end
        tau(j)=sum(x.^(2))./sum(x.*y);
    end
    taus(i)=tau(1);
    end
    %taus_memory=tall(taus);
    
for i=1:m
uncalib_co_values (i)= ((del_P(i)/Period(i)) + (MAP(i)/taus(i)));
TPR(i)=taus(i)./4.6023;
end

x=zeros(m,2);
x(:,1)=uncalib_co_values;
x(:,2)=TPR;