% frequency-wise envelope correlation changes throughout recovery
% Step 1: Load in Calibration data for all patients; first 5 and last 5 runs
% Step 2: Compute power envelopes for individual frequencies ~1-25 Hz
% Step 3: Calculate envelope correlations between ipsi and contralateral
% electrodes
% Step 4: Plot average correlation at each frequency early and late for
% each patient (10 plots) and average all patients together at 2 timepoints
% (1 plot to summarize)

clear
close all
%% load data
%preprocessed calibration data (intro rest period)
load('D:\Dinal\OrganizedCalibrationData\FrontalCal.mat');
load('D:\Dinal\OrganizedCalibrationData\MotorCal.mat');
load('D:\Dinal\OrganizedCalibrationData\TemporalCal.mat');

%% convolve wavelets with signal, calculate correlations
subjects = [1,3,7,9,11,12,13,17,20,22];
feature_channel = [6,4,6,6,6,4,4,4,6,4];
progressbar('Subjects','Runs')

s_i = 1;
early = zeros(10,30,5);
late = zeros(10,30,5);

% populate early matrix and late matrix
% i.e. early = patient x frequency x run

%find appropriate runs
for cur_subj = subjects
    trim_i = 1;
   for run = 1:size(FullRawCalMat{cur_subj},3)
       if (sum(squeeze(FullRawCalMat{cur_subj}(:,1,run)))==0) || (sum(squeeze(FullRawCalMat{cur_subj}(:,2,run)))==0)
           fprintf('Run %d in Subject %d Removed - all zeros in Motor electrode\n',run,cur_subj)
       elseif(sum(squeeze(FrontalRawCalMat{cur_subj}(:,1,run)))==0) || (sum(squeeze(FrontalRawCalMat{cur_subj}(:,2,run)))==0)
           fprintf('Run %d in Subject %d Removed - all zeros in Frontal electrode\n',run,cur_subj)
       elseif (sum(squeeze(TemporalRawCalMat{cur_subj}(:,1,run)))==0) || (sum(squeeze(TemporalRawCalMat{cur_subj}(:,2,run)))==0)
           fprintf('Run %d in Subject %d Removed - all zeros in Temporal electrode\n',run,cur_subj)
       else
           trimCal{cur_subj}(:,:,trim_i) = FullRawCalMat{cur_subj}(:,:,run);
           trimFrontal{cur_subj}(:,:,trim_i) = FrontalRawCalMat{cur_subj}(:,:,run);
           trimTemporal{cur_subj}(:,:,trim_i) = TemporalRawCalMat{cur_subj}(:,:,run);
           trim_i = trim_i+1;
       end
   end
end

%Following commented out blocks to be used if bad runs are not the same
%across all electrodes

% for cur_subj = subjects
%     trim_i = 1;
%    for run = 1:size(FullRawCalMat{cur_subj},3)
%        if (sum(squeeze(FullRawCalMat{cur_subj}(:,1,run)))==0) || (sum(squeeze(FullRawCalMat{cur_subj}(:,2,run)))==0)
%            fprintf('Run %d in Subject %d Removed - all zeros\n',run,cur_subj)
%        else
%            trimTemporal{cur_subj}(:,:,trim_i) = TemporalRawCalMat{cur_subj}(:,:,run);
%            trim_i = trim_i+1;
%        end
%    end
% end

% for cur_subj = subjects
%     trim_i = 1;
%    for run = 1:size(FrontalRawCalMat{cur_subj},3)
%        if (sum(squeeze(FrontalRawCalMat{cur_subj}(:,1,run)))==0) || (sum(squeeze(FrontalRawCalMat{cur_subj}(:,2,run)))==0)
%            fprintf('Run %d in Subject %d Removed - all zeros\n',run,cur_subj)
%        else
%            trimFrontal{cur_subj}(:,:,trim_i) = FrontalRawCalMat{cur_subj}(:,:,run);
%            trim_i = trim_i+1;
%        end
%    end
% end

%  band-limited effects

% beta = allFreqs(:,13:30);
% avgBetaCorr = mean(beta,2);
% ARATchange = [13 1 8 6 2 6 1 7 14 6];
% ARAT_initial = [16.5 6 4 8 32 14 5 5 29.5 10.5]';
% ARAT_change = [12.5 3 8 8 2 7 1 7 13.5 5.5]';
paperARATChange = [12.5 0.5 8 6 2 6 1 7 13.5 5.5]';

% [h p] = corr(avgBetaCorr, ARAT_change);

MDiffs = CorrChange(trimCal);
FDiffs = CorrChange(trimFrontal);
TDiffs = CorrChange(trimTemporal);

Mdelta = mean(MDiffs(:,1:3),2);
Fdelta = mean(FDiffs(:,1:3),2);
Tdelta = mean(TDiffs(:,1:3),2);

Mtheta = mean(MDiffs(:,4:7),2);
Ftheta = mean(FDiffs(:,4:7),2);
Ttheta = mean(TDiffs(:,4:7),2);

Malpha = mean(MDiffs(:,8:12),2);
Falpha = mean(FDiffs(:,8:12),2);
Talpha = mean(TDiffs(:,8:12),2);

Mbeta = mean(MDiffs(:,13:30),2);
Fbeta = mean(FDiffs(:,13:30),2);
Tbeta = mean(TDiffs(:,13:30),2);

figure

subplot(1,3,1)
plot(Malpha,paperARATChange,'ro')
title('Motor - Alpha')
xlabel('Change in Alpha Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Malpha,paperARATChange);
sprintf('Motor Electrode Alpha Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

subplot(1,3,2)
plot(Falpha,paperARATChange,'bo')
title('Frontal - Alpha')
xlabel('Change in Alpha Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Falpha,paperARATChange);
sprintf('Frontal Electrode Alpha Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

subplot(1,3,3)
plot(Talpha,paperARATChange,'go')
title('Temporal - Alpha')
xlabel('Change in Alpha Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Talpha,paperARATChange);
sprintf('Temporal Electrode Alpha Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

figure

subplot(1,3,1)
plot(Mbeta,paperARATChange,'ro')
title('Motor - Beta')
xlabel('Change in Beta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Mbeta,paperARATChange);
sprintf('Motor Electrode Beta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

subplot(1,3,2)
plot(Fbeta,paperARATChange,'bo')
title('Frontal - Beta')
xlabel('Change in Beta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Fbeta,paperARATChange);
sprintf('Frontal Electrode Beta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

subplot(1,3,3)
plot(Tbeta,paperARATChange,'go')
title('Temporal - Beta')
xlabel('Change in Beta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Tbeta,paperARATChange);
sprintf('Temporal Electrode Beta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

figure

subplot(1,3,1)
plot(Mdelta,paperARATChange,'ro')
title('Motor - Delta')
xlabel('Change in Delta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Mdelta,paperARATChange);
sprintf('Motor Electrode Delta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

subplot(1,3,2)
plot(Fdelta,paperARATChange,'bo')
title('Frontal - Delta')
xlabel('Change in Delta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Fdelta,paperARATChange);
sprintf('Frontal Electrode Delta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

subplot(1,3,3)
plot(Tdelta,paperARATChange,'go')
title('Temporal - Delta')
xlabel('Change in Delta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Tdelta,paperARATChange);
sprintf('Temporal Electrode Delta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

figure

subplot(1,3,1)
plot(Mtheta,paperARATChange,'ro')
title('Motor - Theta')
xlabel('Change in Theta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Mtheta,paperARATChange);
sprintf('Motor Electrode Theta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

subplot(1,3,2)
plot(Ftheta,paperARATChange,'bo')
title('Frontal - Theta')
xlabel('Change in Theta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Ftheta,paperARATChange);
sprintf('Frontal Electrode Theta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

subplot(1,3,3)
plot(Ttheta,paperARATChange,'go')
title('Temporal - Theta')
xlabel('Change in Theta Connectivity')
ylabel('Change in ARAT')
[r p] = corr(Ttheta,paperARATChange);
sprintf('Temporal Electrode Theta Connectivity Change Correlation with ARAT Change: \nr = %5.4f, p = %5.4f',r,p)

% function [ span ] = fwhm2span( cf,fwhm )
% %FWHM2SPAN Summary of this function goes here
% %   Detailed explanation goes here
% 
% span= (cf./(2*pi))./(fwhm./ (2*sqrt(2*log(2))));
% 
% 
% end

function corrDiffs = CorrChange(signalMat)

subjects = [1,3,7,9,11,12,13,17,20,22];
feature_channel = [6,4,6,6,6,4,4,4,6,4];
progressbar('Subjects','Runs')

s_i = 1;
early = zeros(10,30,5);
late = zeros(10,30,5);

% Gabor wavelets
fwhm = 0.8;
cf = 1:30; %center frequencies
span = fwhm2span(cf,fwhm);

for cur_subj = subjects
    runNum = size(signalMat{cur_subj},3);
    runsToCount = horzcat(1:5,(runNum-4):runNum);
    run_i = 1;
    for run = runsToCount
        cur_sig = squeeze(signalMat{cur_subj}(:,:,run));
        out = gabor_response_span(cur_sig,cf,span, 256);
        outEnv = squeeze((abs(out(:,:,:))).^2);
        
        if run_i<=5
            early(cur_subj,:,run_i) = diag(corr(squeeze(outEnv(:,:,1)'),squeeze(outEnv(:,:,2)')));
        else
            late(cur_subj,:,run_i-5) = diag(corr(squeeze(outEnv(:,:,1)'),squeeze(outEnv(:,:,2)')));
        end
        
        progressbar([],run_i/10)
        run_i = run_i + 1;
    end
    progressbar(s_i/10)
    s_i=s_i+1;
end

early_z = 0.5*(log(1+early) - log(1-early));
late_z = 0.5*(log(1+late) - log(1-late));
%% Plot for each patient (and average and total)

% figure
q=1;

for cur_subj = subjects
%     subplot(2,5,q)
    
    eAvg = mean(squeeze(early_z(cur_subj,:,:)),2)';
    lAvg = mean(squeeze(late_z(cur_subj,:,:)),2)';
    
    eStdev = std(squeeze(early_z(cur_subj,:,:)),0,2)';
    lStdev = std(squeeze(late_z(cur_subj,:,:)),0,2)';
    
%     boundedline(cf,eAvg,eStdev,'-bx','alpha');
%     hold on
%     boundedline(cf,lAvg,lStdev,'-rx','alpha');
%     hold off
    
%     alpha(cur_subj,:) = lAvg(8:12) - eAvg(8:12);
    allFreqs(cur_subj,:) = lAvg - eAvg;
    
    q=q+1;
end

allFreqs = allFreqs(subjects,:);
corrDiffs = allFreqs;
end