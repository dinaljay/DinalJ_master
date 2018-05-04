close all

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
dataDir = '\\bciserver2\Kenny\data\Stroke-motor\2000Samples_14kForCalibration\grouped';
saveDir = '\\bciserver2\Dinal\Lag Array Data';
saveFile = 'lagArrayAmpEnv.mat';
subjects = [1 3 7 9 11 12 13 17 22];

cellArray = cell(numel(subjects),numel(lobeList),numel(stageList),10);


cf = [2, 10, 21];
fwhm = [2, 4, 8];
meta.Order = {'Patient','Lobe','Stage','Run'};


for ptInd = 1:numel(subjects)
    
    ptID = subjects(ptInd);
%     progressbar('Subjects')
%     progressbar('Lobe Data Upload')
    disp(num2str(ptID));
    for lobeInd = 1:numel(lobeList)
        lobe = lobeList{lobeInd};
        disp(['  Lobe ' lobe]);
        lobeData = [];
        
        % find which file is relevant
        
        D = dir(dataDir); D(1:2) = [];
        for i = 1:numel(D)
            underInd = strfind(D(i).name,'_');
            filePt = D(i).name(underInd+1:end-4);
            filePt = str2double(filePt);
            if filePt == ptID && contains(D(i).name,lobe)
                load(fullfile(dataDir,D(i).name)); % output = data
                lobeData = [lobeData data];
            end
        end

        % load the relevant file (output = data)
        for stageInd = 1:numel(stageList)
            disp(['    Stage # ' num2str(stageInd)]);
            stageT = tic;
            stageData = [];
            
            % get data for this stage
            
            for i = 1:numel(lobeData)
                if strcmp(lobeData(i).Stage,stageList{stageInd})
                   stageData = [stageData lobeData(i)];
                end
            end

            
            % filter by ipsi
            ipsiData = [];
            contraData = [];
            
             for i = 1:numel(stageData)
                if stageData(i).Ipsi
                   ipsiData = [ipsiData stageData(i)];
                else
                   contraData = [contraData stageData(i)];
                end
             end

            % run lag analysis for each run
            maxRunNum = ipsiData(1).MaxRun;
            runsToCount = horzcat(1:5,(maxRunNum-4):maxRunNum);
            
            if ptID == 1
                runsToCount(1:5) = [1,3,4,5,6];
            elseif ptID == 9
                runsToCount(1:5) = [1,2,3,5,6];
            elseif ptID == 11
                runsToCount(1:5) = [1,3,4,6,7];
            elseif ptID == 12
                runsToCount(6:10) = [58,59,60,61,63];
            elseif ptID == 17
                runsToCount(6:10) = [69,70,71,72,74];
            end
             
            for runInd = 1:10
                % get the actual run number (output = run)
                run = runsToCount(runInd);
                
                % find the ipsi that is that run
                ipsiRun = [];
                for i = 1:numel(ipsiData)
                    if ipsiData(i).Run == run
                        ipsiRun = [ipsiRun ipsiData(i)];
                    end
                end
                
                % find the contra that is that run
                contraRun = [];
                for i = 1:numel(contraData)
                    if contraData(i).Run == run
                        contraRun = [contraRun contraData(i)];
                    end
                end
                
                lagArray = nan(numel(ipsiRun),numel(cf));
                for block = 1:numel(ipsiRun)
                    % find the contra block that corresponds to the ipsi
                    % block
                    ipsiRunBlock = ipsiRun(block);
                    startTime = ipsiRun(block).Time(1);
                    for contraBlock = 1:numel(contraRun)
                        if startTime == contraRun(contraBlock).Time(1)
                            contraRunBlock = contraRun(contraBlock);
                        end
                    end
                    
                    % lag analysis on ipsi and contra
                    for band = 1:numel(cf)
                        f = FrequencyDataObtainer();
                        f.FrequencyCenter = cf(band);
                        f.FrequencyCenterRange = fwhm(band);
                        ipsiBandData = f.obtain(ipsiRunBlock);
                        contraBandData = f.obtain(contraRunBlock);
                        ipsiBandData = ipsiBandData.raw();
                        contraBandData = contraBandData.raw();
                        
                        ipsiEnv = (abs(ipsiBandData).^2);        %Contralateral amplitude envelope calculation
                        contraEnv = (abs(contraBandData).^2);    %Ipsilateral amplitude envelope calculation
                        
                        [r, lagTemp] = xcorr(ipsiEnv,contraEnv); %cross correlational analysis of envelopes
                        [~,I] = max(abs(r));
                        lagArray(block,band) = lagTemp(I)/ipsiRunBlock.samplingRate();
                    end
                end
                cellArray{ptInd,lobeInd,stageInd,runInd} = lagArray;

            end
            
            stageT = toc(stageT);
            disp(['    Took ' num2str(stageT) ' seconds.']);
        end
    end
end

% save(fullfile(saveDir,saveFile),'cellArray','cf','fwhm','meta');

%% Create mean and median arrays

for ptInd = 1:numel(subjects)
    for lobeInd = 1:numel(lobeList)
        for stageInd = 1:numel(stageList)
            for runInd = 1:10
                                
                meanDeltaArray{ptInd, lobeInd, stageInd, runInd} = mean(cellArray{ptInd,lobeInd,stageInd,runInd}(:,1));
                meanAlphaArray{ptInd, lobeInd, stageInd, runInd} = mean(cellArray{ptInd,lobeInd,stageInd,runInd}(:,2));
                meanBetaArray{ptInd, lobeInd, stageInd, runInd} = mean(cellArray{ptInd,lobeInd,stageInd,runInd}(:,3));
                
                medianDeltaArray{ptInd, lobeInd, stageInd, runInd} = median(cellArray{ptInd,lobeInd,stageInd,runInd}(:,1));
                medianAlphaArray{ptInd, lobeInd, stageInd, runInd} = median(cellArray{ptInd,lobeInd,stageInd,runInd}(:,2));
                medianBetaArray{ptInd, lobeInd, stageInd, runInd} = median(cellArray{ptInd,lobeInd,stageInd,runInd}(:,3));
                
            end
        end
        
    end
end

% Histogram plot of number of patients with specific lag times
pVal = nan(numel(lobeList),numel(stageList));
for lobeInd = 1:numel(lobeList)
    for stageInd = 1:numel(stageList)
        
        figure
        val1 = [];
        for pt = 1:9
            for run = 1:5
                val1 = [val1 meanDeltaArray{pt,lobeInd,stageInd,run}];
            end
        end
        histogram(val1,'BinWidth',0.05)
        title('Beta Band in Motor Cortex')
        title(sprintf('%s Cortex in %s Stage using mean lag time across trials in delta band', lobeList{lobeInd}, stageList{stageInd}))
        xlabel('lag with bin width of 0.001 seconds')
        ylabel('No. of runs')
        
        hold on
        
        val2 = [];
        for pt = 1:9
            for run = 6:10
                val2 = [val2 meanDeltaArray{pt,lobeInd,stageInd,run}];
            end
        end
        
        text(0.1,0.8,['p = ' num2str(p)],'Units','normalized');
        histogram(val2,'BinWidth',0.05)
        legend('First 5 runs','Last 5 runs')
        hold off
        
        [h,p] = ttest2(val1,val2);
        %         p = ranksum(val1,val2);
        pVal(lobeInd,stageInd) = p;
    end
    
end

pVal = nan(numel(lobeList),numel(stageList));
for lobeInd = 1:numel(lobeList)
    for stageInd = 1:numel(stageList)
        
        figure
        val1 = [];
        for pt = 1:9
            for run = 1:5
                val1 = [val1 meanAlphaArray{pt,lobeInd,stageInd,run}];
            end
        end
        histogram(val1,'BinWidth',0.5)
        title('Beta Band in Motor Cortex')
        title(sprintf('%s Cortex in %s Stage using mean lag time across trials in delta band', lobeList{lobeInd}, stageList{stageInd}))
        xlabel('lag with bin width of 0.001 seconds')
        ylabel('No. of runs')
        
        hold on
        
        val2 = [];
        for pt = 1:9
            for run = 6:10
                val2 = [val2 meanAlphaArray{pt,lobeInd,stageInd,run}];
            end
        end
        
        text(0.1,0.8,['p = ' num2str(p)],'Units','normalized');
        histogram(val2,'BinWidth',0.5)
        legend('First 5 runs','Last 5 runs')
        hold off
        
        [h,p] = ttest2(val1,val2);
        %         p = ranksum(val1,val2);
        pVal(lobeInd,stageInd) = p;
        
    end
    
end

pVal = nan(numel(lobeList),numel(stageList));
for lobeInd = 1:numel(lobeList)
    for stageInd = 1:numel(stageList)
        
        figure
        val1 = [];
        for pt = 1:9
            for run = 1:5
                val1 = [val1 meanBetaArray{pt,lobeInd,stageInd,run}];
            end
        end
        histogram(val1,'BinWidth',0.01)
        title('Beta Band in Motor Cortex')
        title(sprintf('%s Cortex in %s Stage using mean lag time across trials in delta band', lobeList{lobeInd}, stageList{stageInd}))
        xlabel('lag with bin width of 0.01 seconds')
        ylabel('No. of runs')
        
        hold on
        
        val2 = [];
        for pt = 1:9
            for run = 6:10
                val2 = [val2 meanBetaArray{pt,lobeInd,stageInd,run}];
            end
        end
        
        [h,p] = ttest2(val1,val2);
        %         p = ranksum(val1,val2);
        pVal(lobeInd,stageInd) = p;
        %
        text(0.1,0.8,['p = ' num2str(p)],'Units','normalized');
        histogram(val2,'BinWidth',0.01)
        legend('First 5 runs','Last 5 runs')
        hold off
        
    end
    
end
