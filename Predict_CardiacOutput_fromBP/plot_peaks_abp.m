%%%%%%%%%%%%%%%%% ICM Question 1 and 2 %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Sejal Ghate, Zixu Han, Yongzhi Sun %%%%%%%%%%%%%%%%

function [] = plot_peaks_abp(hours, subj_name)

    % Need folder of Physionet Toolkit
    

    % Need ABP file of subject
    path_abp = dir((fullfile(subj_name,'*_ABP.txt')));
    ABP = table2array(readtable(((fullfile(subj_name,path_abp.name)))));

    % Need n file of subject
    path_n = dir((fullfile(subj_name,'*n.txt')));    
    n_data = readtable((fullfile(subj_name,path_n.name)));

    % 
    time_s = hours*10*10;
    start_dur = find(ABP(:,1)>=time_s);
    t1 = start_dur(1);
    new_abp = ABP(t1:t1+1990,2);
    figure();
    plot(new_abp);
    title(sprintf("First 20 peaks of ABP wave starting at the %d th hour", hours));
    xlabel("Time (s)");
    ylabel("Amplitude");

    %% Obtain onset time (x hours)

    r = wabp(new_abp);
    r_feat = abpfeature(new_abp',r);

    figure();
    plot(1:1:1991,new_abp,'LineWidth',1.5);
    hold on;
    sz = 100;

    for i = 1:length(r)
        if i < length(r)
            scatter(r(i,1),new_abp(r(i)),sz,'*','k','LineWidth',2);
            scatter(r_feat(i,9),new_abp(r_feat(i,9)),sz,'x','r','LineWidth',2);
            scatter(r_feat(i,11),new_abp(r_feat(i,11)),sz,'o','b','LineWidth',2);
        end
    end

    xlabel('Time Index','FontSize',12,'FontWeight','bold');
    ylabel('ABP (mmHg)','FontSize',12,'FontWeight','bold');
    title(sprintf("Visualization of ABP at %d hours and analyzed features",hours),'FontSize',14,'FontWeight','bold');
    legend('ABP Waveform','Onset point','End of systole (0.3 beatperiod ^{0.5} method)','End of systole (lowest nonnegative slope method)');

end