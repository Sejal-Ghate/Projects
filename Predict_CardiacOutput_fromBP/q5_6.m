clc;
clear all;
close all;

% run this command only the first time
%addpath(genpath('PhysioToolkitCardiacOutput_MatlabCode'));
cd('C:\Users\sejal\Desktop\ICM_Physiome\Project\CPA_Team9_SG_ZH_YS');
%% Question 5:
% Enter the subsaject name here: input should be string

%%%%%%%%%% In question 1: subject number 20  %%%%%%%%%%%%%
subj_name = "s00020";

% Need ABP file of subject
path_abp = dir((fullfile(subj_name,'*_ABP.txt')));
ABP = table2array(readtable(((fullfile(subj_name,path_abp.name)))));

% Need n file of subject
path_n = dir((fullfile(subj_name,'*n.txt')));
n_data = readtable((fullfile(subj_name,path_n.name)));



hours = 12;
end_time = hours*60*60; %the time of first 12 hours (in seconds)
end_index = find(ABP(:,1)==end_time); %the index where the time (the first column) is the end of first 12 hours



abp_hrs = ABP(1:end_index,:);%ABP from first 12 hours
n_data_hrs_idx = find(n_data.ElapsedTime == end_time); %find the index when 12 hour ends
n_data_hrs = n_data(1:n_data_hrs_idx,:); % subjdata in first 12 hours


onset_times = wabp(abp_hrs(:,2));
r_feat = abpfeature(abp_hrs(:,2)',onset_times);


%%

MAP = r_feat(:,6);
Tn = r_feat(:,7);
alpha = 2;
DAP = r_feat(:,4);
del_Pi = diff(r_feat(:,5));
del_P = [del_Pi; del_Pi(end)];

window = 21;
taus = zeros(size(r_feat,1),1);

y = zeros(window,1);
x = zeros(window,1);

for i = 1: size(r_feat,1)
    if i < ((window - 1)/2) || i >= size(r_feat,1) - ((window - 1)/2)
        continue
    else
        for j = 1 : window   
            den_y(j) = Tn(i+j-1);
            num_y(j) = alpha*(MAP(i+j-1) - DAP(i+j-1)) - del_P(i+j-1);
            y(j) = num_y(j)/den_y(j);
            x(j) = MAP(i+j-1);
        end
        
    end
    b = sum(x.*y)/sum(x.^2);
    taus(i) = 1/b;
    
end

%% 

nearestfun = @(taus) interp1(find(taus),taus(taus~=0),(1:length(taus))','nearest','extrap');
full_taus = 0.5*(nearestfun(taus) + flip(nearestfun(flip(taus))));

uncalib_co_values = ((del_P./Tn) + (MAP./full_taus));

%calib_co = Cn .* (uncalib_co_values);

%% Calibration

cotd_idx_hrs = find(n_data_hrs.CO ~= 0); %find index of cotd that are not equal to zero in first 12 hours
cotd_hrs = n_data_hrs.CO(cotd_idx_hrs); %values of non-zero cotd in first 12 hours
cotd_hrs_time = n_data_hrs.ElapsedTime(cotd_idx_hrs); %time when cotd was measured in first 12 hours
cotd_12_time_hr = (cotd_hrs_time)./3600; %time when cotd measured in first 12 hours are converted to hours

% use cotd_hrs_time to check time indexes for non zero tco values


















