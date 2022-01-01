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
beat_q = jSQI(r_feat,onset_times,abp_hrs);

%%

MAP = r_feat(:,6);
Tn = r_feat(:,7)/(60*125);
alpha = 2;
DAP = r_feat(:,4);
del_Pi = diff(r_feat(:,5));
del_P = [r_feat(1,5); del_Pi];
beta = zeros(size(r_feat,1),1);
window = 21;
taus = zeros(size(r_feat,1),1);
den_y = []; num_y = [];



for i = 1: size(r_feat,1)
    if i < ((window + 1)/2) || i >= size(r_feat,1) - ((window + 1)/2)
        continue
    else
        
        l = (i - ((window-1)/2));
        %disp(l);
        u = (i+((window-1)/2));
        %disp(u);
        x = MAP(l:u);
        den_y = Tn(l:u);
        num_y = alpha*(MAP(l:u) - DAP(l:u)) - del_P(l:u);
        y = num_y./den_y;

        taus(i) = 1/(sum(x.*y)./sum(x.^2));
    end
    
end

nearestfun = @(taus) interp1(find(taus),taus(taus~=0),(1:length(taus))','nearest','extrap');
full_taus = 0.5*(nearestfun(taus) + flip(nearestfun(flip(taus))));

betas = 1./full_taus;
uncalib_co_values = ((del_P./Tn) + (MAP./full_taus));
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


[co, to, told, fea] = estimateCO_v3(onset_times,r_feat,beat_q,5,0);

to_s = to*60;

index_time = find(to_s == cotd_hrs_time); %the index of the first co_uncalibrated from co when the first tco occurs in n_data
est_co_uncalibrated = co(index_time); %the value of uncalibrated co when the firsted cotd is measured

k = tco1/est_co_uncalibrated; %get calibration factor k=first cotd measured/the corresponding uncalibrated estimated co at the same time

est_co_calibrated = k*co; %input k to the L-function to obtain calibrated co
