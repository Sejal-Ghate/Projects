function [Cn, final_calibration_factor] = parlikar(hours, subj_name)

    % Need ABP file of subject
    path_abp = dir((fullfile(subj_name,'*_ABP.txt')));
    ABP = table2array(readtable(((fullfile(subj_name,path_abp.name)))));

    % Need n file of subject
    path_n = dir((fullfile(subj_name,'*n.txt')));
    n_data = readtable((fullfile(subj_name,path_n.name)));
    
    %hours = 12;
    end_time = hours*60*60; %the time of first 12 hours (in seconds)
    end_index = find(ABP(:,1)==end_time); %the index where the time (the first column) is the end of first 12 hours

    abp_hrs = ABP(1:end_index,:);%ABP from first 12 hours
    n_data_hrs_idx = find(n_data.ElapsedTime == end_time); %find the index when 12 hour ends
    n_data_hrs = n_data(1:n_data_hrs_idx,:); % subjdata in first 12 hours

    onset_times = wabp(abp_hrs(:,2));
    r_feat = abpfeature(abp_hrs(:,2)',onset_times);
    beat_q = jSQI(r_feat,onset_times,abp_hrs);

    MAP = r_feat(:,6);
    Period = r_feat(:,7);
    Pdias = r_feat(:,4);
    PP = r_feat(:,5);
    
    [co_uncal, to_par, told_par, fea_par, taus_tpr] = estimateCO_v3_Parlikar(onset_times,r_feat,beat_q,14,15);
    
    %% Values of N_data that are non Zero and ABP within 12 hours:

    n_data_hrs_idx = find(n_data.ElapsedTime == end_time); %find the index when 12 hour ends
    n_data_hrs = n_data(1:n_data_hrs_idx,:); % subjdata in first 12 hours

    cotd_idx_hrs = find(n_data_hrs.CO ~= 0); %find index of cotd that are not equal to zero in first 12 hours
    cotd_hrs = n_data_hrs.CO(cotd_idx_hrs); %values of non-zero cotd in first 12 hours
    cotd_hrs_time = n_data_hrs.ElapsedTime(cotd_idx_hrs); %time when cotd was measured in first 12 hours
    cotd_12_time_hr = (cotd_hrs_time)./3600; %time when cotd measured in first 12 hours are converted to hours
   
    %% find first non 0 value of TCO and time corresponding to first TCO

    tco1 = cotd_hrs(1);%the first non-zero cotd value (used for calibration)
    tco1_time = n_data.ElapsedTime(n_data.CO == tco1);% the time when the first measurement of cotd was conducted
%% Calibration of estimated co values from estimate function

    to_s = round(to_par*60); %convert time to seconds: to is orignally in minutes

    index_time = find(to_s == tco1_time); %the index of the first co_uncalibrated from co when the first tco occurs in n_data
    est_co_uncalibrated = co_uncal(index_time); %the value of uncalibrated co when the firsted cotd is measured
    %% Calculation of normal calibration

   Cn = tco1/est_co_uncalibrated; %get calibration factor k=first cotd measured/the corresponding uncalibrated estimated co at the same time
   est_co_calibrated = Cn*co_uncal; %input k to the L-function to obtain calibrated co
   
   
    %% non zero tco values table

    non0_idx = find(n_data_hrs.CO ~= 0 );%find the index where tco is not equal to zero
    non0_values = n_data_hrs(non0_idx,:);%subjdata.n when tco were measured
    
    %for stem plot-pp, map and hr suggests the value where tco were
    %measured at the same time
    non0_pp = non0_values.ABPSys - non0_values.ABPDias;% at the above calculated time array->pulse pressure =ABPS-ABPD
    non0_map = non0_values.ABPMean;% at the same time as tco mearements, MAP=ABPmean
    non0_hr = non0_values.HR;%at the same time as tco mearements,HR directly shown in n data
    
    
    %% find the index of time in to_s series corresponding to cotd and the
    %calibratoin factor in each co-cotd pair
    L=length(cotd_hrs_time);%how mant times cotd was measured in first 12 hours
    for i=1:L
        index_time(i)=find(to_s==cotd_hrs_time(i),1);
        calibration_factor(i)=cotd_hrs(i)/co_uncal(index_time(i));
        map_fea(i,1) = non0_map(i)/fea_par(index_time(i),6);
    end
    
    calibration_factor_mean=mean(calibration_factor);
%     MAP_mean=mean(non0_map);
    MAP_mean = mean(map_fea);
    for i=1:L
%         beta1_num=(non0_map(i)-MAP_mean)*(calibration_factor(i)-calibration_factor_mean);
%         beta1_de=(non0_map(i)-MAP_mean)^(2);
        beta1_num=(map_fea(i)-MAP_mean)*(calibration_factor(i)-calibration_factor_mean);
        beta1_de=(map_fea(i)-MAP_mean).^(2);

    end
    

    
   %% Calculation of robust calibration
   
    beta1=sum(beta1_num)./(beta1_de);
    betazero=calibration_factor_mean-beta1.*MAP_mean;
    final_calibration_factor=betazero+beta1.*MAP_mean;
    
    est_co_gamma=final_calibration_factor.*co_uncal;
    
   %% Calculate TPR
   
   %map values
   map_tpr = fea_par(:,6);
   
   % pp values
   PP_tpr = fea_par(:,5);
   PP_diff = diff(PP_tpr);
   PP_tpr = [fea_par(1);PP_diff];
   
   % period values
   Tn_tpr = fea_par(:,7);
   
   % tpr formula
   tpr_val = map_tpr./(est_co_gamma - (final_calibration_factor.*(PP_tpr./Tn_tpr)));
   
   %graph
   
   tpr_val = tpr_val.*(60/1000);
   figure(1);
   title(sprintf("Subject # %s",subj_name));
   xlabel("Time in hrs");
   ylabel("TPR(mmhg/(mL/s))");
   plot(to_par./60,tpr_val);
   
    %%  plot Cotd in first 12 hours with a function of time in hour for C2 and gamma Calibration
      
    figure(2);
    tiledlayout(2,1);
    ax1=nexttile;
    
    % continuous plots
    plot(ax1,(to_s)./3600,est_co_calibrated); %plot calibrated co with a function of time in first hours
    hold on;
    % stem plot
    stem(ax1,cotd_12_time_hr,cotd_hrs,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','#D95319','LineWidth',1.5); %stem plot of cotd
    
    title(ax1,'Estimated CO versus thermodilution CO measurements - C2 Calibration');
    ylabel(ax1,'CO');
    hold off;
    
    ax2=nexttile;
    
    % continuous plots
    plot(ax2,(to_s)./3600,est_co_gamma); %plot calibrated co with a function of time in first hours
    hold on;
    % stem plot
    stem(ax2,cotd_12_time_hr,cotd_hrs,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','#D95319','LineWidth',1.5); %stem plot of cotd
    
    title(ax2,'Estimated CO versus thermodilution CO measurements - Gamma calibration');
    ylabel(ax2,'CO');
    hold off;   
    
end