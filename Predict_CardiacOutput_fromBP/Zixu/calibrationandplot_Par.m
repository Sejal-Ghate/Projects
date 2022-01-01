%% calibration of estimated_co based on Par alg

clc;
clear all;
close all;

% run this command only the first time
%addpath(genpath('PhysioToolkitCardiacOutput_MatlabCode'));
addpath(genpath('PhysioToolkitCardiacOutput_MatlabCode'));
addpath 'C:\Users\lenovo\Downloads';
load('taus.mat');
subj_name = "s00020";

%% ABP file and n txt file in first 12 hours
% Need ABP file of subject
    path_abp = dir((fullfile(subj_name,'*_ABP.txt')));
    ABP = table2array(readtable(((fullfile(subj_name,path_abp.name)))));

    % Need n file of subject
    path_n = dir((fullfile(subj_name,'*n.txt')));    
    n_data = readtable((fullfile(subj_name,path_n.name)));
    
    
    end_time = 12*60*60; %the time of first 12 hours (in seconds)
    end_index = find(ABP(:,1)==end_time); %the index where the time (the first column) is the end of first 12 hours
    abp_hrs = ABP(1:end_index,:);%ABP from first 12 hours
    %% Applying estimateCO_v3 on data

    % to on
    t_on = wabp(abp_hrs(:,2)); %the onset times of first 12 hours of ABPs

    % getting features
    features = abpfeature(abp_hrs(:,2),t_on);%get features of ABP_12hours

    % getting beats
    beats = jSQI(features,t_on,abp_hrs(:,2));%beats quality

    %using the fifth estimator
    est=14; 

    %set filter order as zero
    filt_order = 0;

    %get uncalibrated estimated co
    [co, to, told, fea,TPR] = estimateCO_v3(t_on,features,beats,est,filt_order);
    %% Calibration 
    
    %Values of N_data that are non Zero and ABP within 12 hours:
    to_s = round(to*60);
    n_data_hrs_idx = find(n_data.ElapsedTime == end_time); %find the index when 12 hour ends
    n_data_hrs = n_data(1:n_data_hrs_idx,:); % subjdata in first 12 hours

    cotd_idx_hrs = find(n_data_hrs.CO ~= 0); %find index of cotd that are not equal to zero in first 12 hours
    cotd_hrs = n_data_hrs.CO(cotd_idx_hrs); %values of non-zero cotd in first 12 hours
    cotd_hrs_time = n_data_hrs.ElapsedTime(cotd_idx_hrs); %time when cotd was measured in first 12 hours
    cotd_12_time_hr = (cotd_hrs_time)./3600; %time when cotd measured in first 12 hours are converted to hours
    
    % non zero tco values table

    non0_idx = find(n_data_hrs.CO ~= 0 );%find the index where tco is not equal to zero
    non0_values = n_data_hrs(non0_idx,:);%subjdata.n when tco were measured
    
    %for stem plot-pp, map and hr suggests the value where tco were
    %measured at the same time
    non0_pp = non0_values.ABPSys - non0_values.ABPDias;% at the above calculated time array->pulse pressure =ABPS-ABPD
    non0_map = non0_values.ABPMean;% at the same time as tco mearements, MAP=ABPmean
    non0_hr = non0_values.HR;%at the same time as tco mearements,HR directly shown in n data
    
    
    %find the index of time in to_s series corresponding to cotd and the
    %calibratoin factor in each co-cotd pair
    L=length(cotd_hrs_time);%how mant times cotd was measured in first 12 hours
    for i=1:L
        index_time(i)=find(to_s==cotd_hrs_time(i),1);
        calibration_factor(i)=cotd_hrs(i)/co(index_time(i));
        i=i+1;
    end
    
    calibration_factor_mean=mean(calibration_factor);
    MAP_mean=mean(non0_map);
    
    for i=1:L
        beta1_num=(non0_map(i)-MAP_mean)*(calibration_factor(i)-calibration_factor_mean);
        beta1_de=(non0_map(i)-MAP_mean)^(2);
    end
    beta1=sum(beta1_num)/(beta1_de);
    betazero=calibration_factor_mean-beta1*MAP_mean;
    final_calibration_factor=betazero+beta1*MAP_mean;
    
    est_co_calibrated=final_calibration_factor.*co;
    
    
    %%  plot Cotd in first 12 hours with a function of time in hour
    figure();
    %tiledlayout(4,1);
    %ax1=nexttile;
    
    subplot(4,1,1);
    % continuous plots
    plot((to_s)./3600,est_co_calibrated(1:length(to_s),1)); %plot calibrated co with a function of time in first hours
    hold on;
    % stem plot
    stem(cotd_12_time_hr,cotd_hrs,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','#D95319','LineWidth',1.5); %stem plot of cotd
    
    title('Estimated CO versus thermodilution CO measurements');
    ylabel('CO');
    hold off;
    
    %% plot of pulse pressure

    %ax2=nexttile;
    subplot(4,1,2);

    %continuous plot-data from feature matrix estimated by "abpfeature"
    plot((to_s)./3600,fea(:,5));%in the output of "abpfeature", the fifth column refers to pulse pressure
    hold on;
    %stem plot
    stem(cotd_12_time_hr,non0_pp,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','#D95319','LineWidth',1.5);

    title('Estimated pulse pressure vs PP measurements');
    ylabel('PP');
    hold off;

    %% plot of MAP

    %ax3=nexttile;
    subplot(4,1,3)

    %continuous plot
    plot((to_s)./3600,fea(:,6));%%in the output of "abpfeature", the 6th column refers to MAP
    hold on;
    %stem plot
    stem(cotd_12_time_hr,non0_map,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','#D95319','LineWidth',1.5);

    title('Estimated mean arterial pressure vs. MAP measurements');
    ylabel('MAP');
    hold off;

    %% plot of HR

    %ax4=nexttile;
    subplot(4,1,4);

    %continuous plot
    plot((to_s)./3600,60*125./fea(:,7));%xcalculating HR based on beat period
    hold on;
    %stem plot
    stem(cotd_12_time_hr,non0_hr,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','#D95319','LineWidth',1.5);    

    title('Estimated heart rate vs HR measurements');
    xlabel('time[hours]');
    ylabel('HR');
    hold off;
    
    %% extra plot of TPR
    figure();
    plot((to_s)./3600,TPR(1:length(to_s),1));
    title('TPR');
    
   

