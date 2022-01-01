%import static training data (700 patients) - this is the basis of the simple model
load('static_data_training.mat');
%The header variable contains the meaning of each column of static_train
%generate simple glm
%define Y = observations which should be loaded from clinical table
Y = static_train(:,2);

%define X = covariate matrix by taking features from table. 
%This currently only uses Gender as a covariate. % uses all?
X = static_train(:,3:7);

%compute glm
[B,dev,stats] = glmfit(X,Y,'binomial');
%construct phat from parameters and X 
Phat = 1./(1+exp(-[ones(size(X,1),1) X]*B)); %equivalent way to compute Phat
%Phat is the estimated probability of sepsis occurence for patients


%plot phat versus patient along with its confidence bounds (1.96*stats.se)
% % Phat_UB = ???
% % Phat_LB = ???

% plot first 30 patients prediction, uncertainty and labels.
figure(1)
plot(Phat(1:30))
hold on
% plot(Phat_LB,'b-')
% hold on
% plot(Phat_UB,'b-')
% hold on
plot(Y(1:30),'r*')
title('Models for Each Patient')

%test performance of models
% [threshold] = test_performance(Phat, Y);
% l=Y;
% [X,Y,T,AUC] = perfcurve(l,Phat, 1);
% figure;
% plot(X,Y)
