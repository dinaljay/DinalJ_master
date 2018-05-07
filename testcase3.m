close all

% cellArray{ptInd,lobeInd,stageInd,runInd}
lobeList = {'Frontal','Temporal','Motor'};
stageList = {'Calibration','Task','Rest'};
freqList = {'Delta','Alpha','Beta'};
ipsiList = {'0','1'};  
paperARATChange = [12.5 0.5 8 6 2 6 1 7 5.5];
subjects = [1 3 7 9 11 12 13 17 22];

for ptInd = 1:numel(subjects)
    for lobeInd = 1:numel(lobeList)
        for stageInd = 1:numel(stageList)
                       
            meanDeltaArray{ptInd, lobeInd, stageInd, runInd} = mean(cellArray{ptInd,lobeInd,stageInd,runInd}(:,1));
            meanAlphaArray{ptInd, lobeInd, stageInd, runInd} = mean(cellArray{ptInd,lobeInd,stageInd,runInd}(:,2));
            meanBetaArray{ptInd, lobeInd, stageInd, runInd} = mean(cellArray{ptInd,lobeInd,stageInd,runInd}(:,3));
            
            medianDeltaArray{ptInd, lobeInd, stageInd, runInd} = median(cellArray{ptInd,lobeInd,stageInd,runInd}(:,1));
            medianAlphaArray{ptInd, lobeInd, stageInd, runInd} = median(cellArray{ptInd,lobeInd,stageInd,runInd}(:,2));
            medianBetaArray{ptInd, lobeInd, stageInd, runInd} = median(cellArray{ptInd,lobeInd,stageInd,runInd}(:,3));
            
            
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
        title(sprintf('%s Cortex in %s Stage using mean lag time across trials in Delta band', lobeList{lobeInd}, stageList{stageInd}))
        xlabel('lag with bin width of 0.05 seconds')
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
        title(sprintf('%s Cortex in %s Stage using mean lag time across trials in Alpha band', lobeList{lobeInd}, stageList{stageInd}))
        xlabel('lag with bin width of 0.5 seconds')
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
        title(sprintf('%s Cortex in %s Stage using mean lag time across trials in Beta band', lobeList{lobeInd}, stageList{stageInd}))
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
