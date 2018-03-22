clear all
close all
%% load data
%preprocessed calibration data (intro rest period)
load('D:\Dinal\OrganizedCalibrationData\FrontalCal.mat');
load('D:\Dinal\OrganizedCalibrationData\MotorCal.mat');
load('D:\Dinal\OrganizedCalibrationData\TemporalCal.mat');

%% convolve wavelets with signal, calculate correlations
subjects = [1,3,7,9,11,12,13,17,20,22];
feature_channel = [6,4,6,6,6,4,4,4,6,4]; %represents the channel recording from ipsilateral ortex (1=4 & 2=6)
ipsichan = [2,1,2,2,2,1,1,1,2,1];
contrachan = [1,2,1,1,1,2,2,2,1,2];

s_i = 1;
early = zeros(10,30,5);
late = zeros(10,30,5);

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

CorrChange(trimCal);
% FDiffs = CorrChange(trimFrontal);
% TDiffs = CorrChange(trimTemporal);

% figure 
% 
% subplot (1,1,1)
% plot(Malpha,r)
% title('cross correlational analysis alpha in motor cortex')
% xlabel('Time(s)')
% % ylabel('Ampltidue')
% 
% subplot (1,1,2)
% plot(Falpha,r)
% title('cross correlational analysis alpha in frontal cortex')
% xlabel('Time(s)')
% % ylabel('Ampltidue')
% 
% subplot (1,1,3)
% plot(Talpha,r)
% title('cross correlational analysis alpha in temporal cortex')
% xlabel('Time(s)')
% ylabel('Ampltidue')

function [] = CorrChange(signalMat)

subjects = [1,3,7,9,11,12,13,17,20,22];
feature_channel = [6,4,6,6,6,4,4,4,6,4];
ipsichan = [2,1,2,2,2,1,1,1,2,1];
contrachan = [1,2,1,1,1,2,2,2,1,2];

s_i = 1;
% early = zeros(10,30,5);
% late = zeros(10,30,5);

% Gabor wavelets
fwhm = 4;
cf = 5:8; %change center frequencies for each band
span = fwhm2span(cf,fwhm);

%loop through valid runs

progressbar('Subjects','Runs')

% lag_final = nan(numel(subjects),10);

for cur_subj = 1:numel(subjects)
    disp(num2str(cur_subj));
    ptID = subjects(cur_subj);
    runNum = size(signalMat{ptID},3);
    runsToCount = horzcat(1:5,(runNum-4):runNum);
    runCountExceptions
    run_i = 1;
    
    runID = 0;
    for run = runsToCount
        runID = runID+1;
        cur_sig = squeeze(signalMat{ptID}(:,:,run));  %define current signal
        out = gpu_gabor_response_span(cur_sig,cf,span,256); %gabor convolution
        outEnv = squeeze((abs(out(:,:,:)))); %amplitude envelope calculation. for power, this would be squared

        ipsienv = outEnv(1,:,ipsichan(cur_subj)); %amplitude envelope for ipsi channel
        contraenv = outEnv(1,:,contrachan(cur_subj)); %amplitude envelope for contra channel

        lag = nan(size(outEnv,1),numel(ipsienv)+numel(contraenv)-1);

        for freq = 1:size(outEnv,1)
            ipsienv = outEnv(freq,:,ipsichan(cur_subj)); %amplitude envelope for ipsi channel
            contraenv = outEnv(freq,:,contrachan(cur_subj)); %amplitude envelope for contra channel
            [r, lagTemp] = xcorr(ipsienv,contraenv); %calculate the cross correlation
            lag(freq,:) = r; %
        end
      
        temp = nanmean(lag,1);
        [~,I] = max(abs(temp),[],2);
        
        figure
        
        plot(lagTemp,temp)
        title(sprintf('Cross correlation for Subject %i and Run #%i for Alpha band in Motor Cortex', cur_subj, run_i))
        xlabel('')
        ylabel('Lag indices')
        legend(sprintf('Maximum is at lag %d', lagTemp(I)))
             
  
%         if run_i<=5
%             [r, lag] = xcorr(ipsienv,contraenv);
%             early(cur_subj,run_i) = max(abs(lag));
% %             early(cur_subj,:,run_i) = diag(corr(squeeze(outEnv(:,:,1)'),squeeze(outEnv(:,:,2)')));
%         else 
%             late(cur_subj,:,run_i-5) = max(abs(lag));
% %             late(cur_subj,:,run_i-5) = diag(corr(squeeze(outEnv(:,:,1)'),squeeze(outEnv(:,:,2)')));
%         end
%         
        progressbar([],run_i/10)
        run_i = run_i + 1;

        
    end
    progressbar(s_i/10)
    s_i=s_i+1;
end

end