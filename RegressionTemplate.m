clear
close all
%% load data
load('D:\Dinal\OrganizedCalibrationData\MotorCal.mat');

%% create wavelets

%example: alpha band
cf = 10;
fwhm = 4;
span = fwhm2span(cf,fwhm);

%% convolve alpha wavelets with signal, calculate power

%initialize some values to use later
subjects = [1,3,7,9,11,12,13,17,20,22];
feature_channel = [6,4,6,6,6,4,4,4,6,4];
feature_frequency = [15,16,19,11,11,9,15,11,11,17];
ARATchange = [12.5 0.5 8 6 2 6 1 7 13.5 5.5]';
contra_channel = [1,2,1,1,1,2,2,2,1,2];

% standard script to remove bad runs
for cur_subj = subjects
    trim_i = 1;
   for run = 1:size(FullRawCalMat{cur_subj},3)
       if (sum(squeeze(FullRawCalMat{cur_subj}(:,1,run)))==0) || (sum(squeeze(FullRawCalMat{cur_subj}(:,2,run)))==0)
           fprintf('Run %d in Subject %d Removed - all zeros in Motor electrode\n',run,cur_subj)
       else
           trimCal{cur_subj}(:,:,trim_i) = FullRawCalMat{cur_subj}(:,:,run);
           trim_i = trim_i+1;
       end
   end
end

% power calculated here for left motor electrode; replace with whatever value of interest
for cur_subj = subjects
    for run = 1:size(trimCal{cur_subj},3)
        cur_sig = squeeze(trimCal{cur_subj}(:,:,run));
        out = gabor_response_span(cur_sig,cf,span, 256);
        outEnv = squeeze((abs(out(:,:,:))).^2);
        freqPower{cur_subj}(run,:) = median(outEnv(:,1),1);
    end
end

% at this point, you should have 1 value per subject, per run, and per
% channel (depending on the metric)

%% Calculate and Plot Individual Level Effects
figure
q = 1;
for cur_subj = subjects
    subplot(2,5,q) %assuming 10 subjects
    y = freqPower{cur_subj}'; %y is value of interest
    y = y(isnan(y)==0);
    x = (1:length(y)); %x is number of runs
    plot(x,y,'o') %plot value at each run for that patient
    lsline %line of best fit (least squares error)
    %define linear regression model
    mdl = fitlm(x,y); % requires statistics and machine learning toolbox, can alternatively use polyfit(x,y,1)
    Rsq(cur_subj) = mdl.Rsquared; %r squared error
    CoefStruct = table2struct(mdl.Coefficients,'ToScalar',true); %pull slope value from table
    RegressionSlope(cur_subj) = CoefStruct.Estimate(2); %2nd coefficient is slope
%     temp=polyfit(x,y,1);
%     b1(cur_subj) = temp(1);
    %plotting labels
    temptitletxt = sprintf('Patient %d, ARATchange=%d',cur_subj,ARATchange(q));
    xlabel('Number of Runs'); ylabel('Value of Interest');
    title(temptitletxt);
    q=q+1;
end

%% Calculate and Plot Group Level Effects
figure
%compare regression slope to motor recovery
RegressionSlope = RegressionSlope(subjects);
plot(RegressionSlope,ARATchange,'bo')
[r p] = corr(RegressionSlope',ARATchange,'type','Spearman')
