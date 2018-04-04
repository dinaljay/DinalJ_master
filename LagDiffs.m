function store = LagDiffs(signalMat)

subjects = [1,3,7,9,11,12,13,17,20,22];
feature_channel = [6,4,6,6,6,4,4,4,6,4];
ipsichan = [2,1,2,2,2,1,1,1,2,1];
contrachan = [1,2,1,1,1,2,2,2,1,2];
lag_store = zeros(10,11);

s_i = 1;

% Gabor wavelets
fwhm = 2; %alpha = 4, beta = 8, delta = 2
cf = 2; %change center frequencies for each band. alpha = 10, beta = 21, delta = 2
span = fwhm2span(cf,fwhm);

progressbar('Subjects','Runs')

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
        cur_sig = squeeze(signalMat{ptID}(:,:,run)); %define current signal
        out = gpu_gabor_response_span(cur_sig,cf,span,256); %gabor convolution
        outEnv = squeeze((out(:,:,:))); %amplitude envelope calculation. for power, this would be squared

        ipsienv = outEnv(:,ipsichan(cur_subj)); %amplitude envelope for ipsi channel
        contraenv = outEnv(:,contrachan(cur_subj)); %amplitude envelope for contra channel

%         lag = nan(size(outEnv,1),numel(ipsienv)+numel(contraenv)-1);
        
  %      for freq = 1:size(outEnv,1)
            ipsienv = outEnv(:,ipsichan(cur_subj)); %amplitude envelope for ipsi channel
            contraenv = outEnv(:,contrachan(cur_subj)); %amplitude envelope for contra channel
            [r, lagTemp] = xcorr(ipsienv,contraenv); %calculate the cross correlation
        %   lag(freq,:) = r; 
    %    end
        
      %  temp = nanmean(lag,1);
        
        [~,I] = max(abs(r));
        lagDiff = lagTemp(I);
           
%         figure
%         
%         plot(lagTemp,r)
%         title(sprintf('Cross correlation for Subject %i and Run #%i for Beta band', cur_subj, run_i));
%         xlabel('Time')
%         ylabel('Cross correlation')
%         ytickformat('%.1f')
%         legend(sprintf('Maximum lag is at %0.2f', lagDiff))
        
%         lag_final(cur_subj,runID,:) = nanmean(lag,1); %average the lags
%         for one subject and one run

        lag_store(cur_subj, run_i) = lagDiff;
         
        progressbar([],run_i/10)
        run_i = run_i + 1;

    end
    lag_store(cur_subj,11) = median(lag_store(s_i,:));
    progressbar(s_i/10)
    s_i=s_i+1;
end
store = lag_store;
end