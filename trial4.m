clear 
close all
%% load data
data_dir = '\\bciserver2\Kenny\data\Stroke-motor\2000Samples_14kForCalibration\grouped';
D = dir(data_dir);

% files = dir('*/mat');

matfiles = dir('*.mat') ;
N = length(matfiles) ;
filenames = cell(N,1) ;

for i = 1:N
   filenames{i} = matfiles(i).name ;
   data = load(matfiles(i).name) ;
 
end


%% Data Block Format
% 1st = pt
% 2nd = Time
% 3rd = TotalTime
% 4th = lobe (frontal, motor, temporal)
% 5th = stage (calibration, task, etc)
% 6th = Ipsi (0, 1)
% 7th = Run
% 8th - MaxRun

%% Separating data based on lobe

lobeList = {'Frontal','Temporal','Motor'};
stageList = {'Calibration','Task','Rest'};
ipsiList = {'0','1'};
subjects = [1 3 7 9 11 12 13 17 22];

cellArray = cell(numel(subjects),numel(lobeList),numel(stageList));

for ptInd = 1:numel(subjects)
    
    ptID = subjects(ptInd);
    progressbar('Subjects')
    s_i = 1;
    
    for lobeInd = 1:numel(lobeList)
        % find which file is relevant
        for i = 1:numel(data)
            if strcmp(data(i).Lobe,LobeList{lobeInd})
                load (data(i));
            end
        end
        
        % load the relevant file (output = data)
        for stageInd = 1:numel(stageList)
            % get data for this stage
            stageData = [];
            for i = 1:numel(data)
                if strcmp(data(i).Stage,stageList{stageInd})
                    stageData = [stageData data(i)];
                end
            end
            
            % filter by ipsi
            ipsiData = [];
            for i = 1:numel(stageData)
                if stageData(i).Ipsi == 1
                    ipsiData = [ipsiData stageData(i)];
                end
            end
            
            % filter by contra
            
            contraData = [];
            for i = 1:numel(stageData)
                if stageData(i).Ipsi == 0
                    contraData = [contraData stageData(i)];
                end
            end
            
            % run lag analysis for each run
            lagArray = nan(10,1);
            maxRunNum = stageData.MaxRun;
            runsToCount = horzcat(1:5,(maxRunNum-4):maxRunNum);
            
            % Gabor wavelets
            fwhm = 2; %alpha = 4, beta = 8, delta = 2
            cf = 2; %change center frequencies for each band. alpha = 10, beta = 21, delta = 2
            span = fwhm2span(cf,fwhm);
                      
            for runInd = runsToCount
                % get the actual run number (output = run)
                                
                % find the ipsi that is that run
                ipsiRun = ipsiData.runNum;
                % find the contra that is that run
                contraRun = contraData.runNum;
                % lag analysis on ipsi and contra
                if (ipsiData.Pt == contraData.Pt)
                    
                    if (ipsiRun(runInd) == contraRun(runInd))
                    
                        cur_sig = squeeze(signalMat{ptID}); %define current signal
                        ipsiEnv = squeeze((gpu_gabor_response_span(cur_sig,cf,span,256))^.2);
                        contraEnv = squeeze((gpu_gabor_response_span(cur_sig,cf,span,256)).^2);
                    
                        [r, lagTemp] = xcorr(ipsiEnv,contraEnv);
                   
                        % save to an array
                    
                        [~,I] = max(abs(r));
                        lagArray(runInd) = lagTemp(I) ;
                    end
                end
            end
            
            % save the lag values to the big cell array
            cellArray{ptInd,lobeInd,stageInd} = lagArray;
            
            progressbar(s_i/10)
            s_i=s_i+1;
        end
    end
end

%% Histogram plot of number of patients with specific lag times

figure
histogram(lag_store(:,1:5),'BinWidth',5)
% title('Beta Band in Motor Cortex')
xlabel('lag with bin width of 5')
ylabel('No. of runs')
hold on

histogram(lag_store(:,6:10),'BinWidth',5)
legend('First 5 runs','Last 5 runs')
hold off

[H,P,CI,STATS] = ttest(lag_store(1:50),lag_store(51:100));

%% Plot change in ARAT against change in lag

figure 
x = lag_store(:,numel(subjects));
y = paperARATChange;
% plot(lag_store(:,numel(subjects)),paperARATChange,'bo');
%lsline

temp = polyfit(x,y,1);
f = polyval(temp,x);
plot(x,y,'o',x,f,'-') 
m = temp(1);
c = temp(2);

title('Alpha Band in Frontal Cortex')
xlabel('Median change in lag')
ylabel('Change in ARAT')
legend('data','linear fit') 

%% Calculating Lag Differences
% 
% function store = LagDiffs(signalMat)
% 
% subjects = [1,3,7,9,11,12,13,17,20,22];
% feature_channel = [6,4,6,6,6,4,4,4,6,4];
% ipsichan = [2,1,2,2,2,1,1,1,2,1];
% contrachan = [1,2,1,1,1,2,2,2,1,2];
% lag_store = [];
% lobe = {'Frontal','Temporal','Motor'};
% stage = {'Calibration','Rest','Task'};
% 
% %arrays storing stage specific information
% 
% dataArray = [];
% % CalArray = [];
% % RestArray = [];
% % TaskArray = [];
% 
% % 1 = Calibration
% % 2 = Rest
% % 3 = Task
% 
% s_i = 1;
% 
% % Gabor wavelets
% fwhm = 2; %alpha = 4, beta = 8, delta = 2
% cf = 2; %change center frequencies for each band. alpha = 10, beta = 21, delta = 2
% span = fwhm2span(cf,fwhm);
% 
% progressbar('Subjects','Runs')
% 
% for i = 1:numel(signalMat)  
%     
%             if contains (signalMat(block).Stage,stage{1})
%             dataArray(1,:) = [dataArray signalMat(block)];
%             end
%             
%             if contains (signalMat(block).Stage,stage{2})
%             dataArray(2,:) = [dataArray signalMat(block)];
%             end
%             
%             if contains (signalMat(block).Stage,stage{3})
%             dataArray(3,:) = [dataArray signalMat(block)];
%             end                          
% end
% 
% for i = 1:size(dataArray,1) 
%     %ptID = data(block).Pt;
%     
% for cur_subj = 1:numel(subjects)   
%     ptID = subjects(cur_subj);
%     
%     runNum = data(block).Run;
%     maxRun = data(block).MaxRun;
%     
%     runsToCount = horzcat(1:5,(maxRun-4):maxRun);
%     runCountExceptions
%     run_i = 1;
%     
%     runID = 0;
%     for run = runsToCount
%         runID = runID+1;
%         cur_sig = squeeze(signalMat{ptID}(:,:,run)); %define current signal
%         out = gpu_gabor_response_span(cur_sig,cf,span,256); %gabor convolution
%         outEnv = squeeze((out(:,:,:))); %amplitude envelope calculation. for power, this would be squared
% 
%         ipsienv = outEnv(:,ipsichan(cur_subj)); %amplitude envelope for ipsi channel
%         contraenv = outEnv(:,contrachan(cur_subj)); %amplitude envelope for contra channel
% 
% 
%         ipsienv = outEnv(:,ipsichan(cur_subj)); %amplitude envelope for ipsi channel
%         contraenv = outEnv(:,contrachan(cur_subj)); %amplitude envelope for contra channel
%         [r, lagTemp] = xcorr(ipsienv,contraenv); %calculate the cross correlation
%         
%         [~,I] = max(abs(r));
%         lagDiff = lagTemp(I);
%            
% %         figure
% %         
% %         plot(lagTemp,r)
% %         title(sprintf('Cross correlation for Subject %i and Run #%i for Beta band', cur_subj, run_i));
% %         xlabel('Time')
% %         ylabel('Cross correlation')
% %         ytickformat('%.1f')
% %         legend(sprintf('Maximum lag is at %0.2f', lagDiff))
%         
% %         lag_final(cur_subj,runID,:) = nanmean(lag,1); %average the lags
% %         for one subject and one run
% 
%         lag_store(cur_subj, run_i) = lagDiff;
%          
%         progressbar([],run_i/10)
%         run_i = run_i + 1;
% 
%     end
%     lag_store(cur_subj,11) = median(lag_store(s_i,:));
%     progressbar(s_i/10)
%     s_i=s_i+1;
% end
% store = lag_store;
% end
% end