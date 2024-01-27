clear all

%% LFP ALIGNED TO DELTA PEAKS HEATMAP

%% Set data paths
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';

% Initialize result arrays
aCOHs = [];
aSPs = [];
aCOHs2 = [];
aPHs = [];
aNT = [];
aAs = [];
aAs2 = [];
Atrace = [];
mm = 0;

for stim_freq = 1:2
    if stim_freq == 1
        % 10Hz
        allSES = {'\608448\01062020\608448', '\608451\01062020\608451', '\608452\01062020\608452', '\608450\03092020\608450', ...
            '\602101\03092020\602101', '\611311\03092020\611311', '\608448\01172020\608448', '\608451\07082020\608451', ...
            '\608452\07082020\608452', '\611111\07082020\611111', '\602101\07082020\602101', '\612535\07082020\612535', ...
            '\615883\02052021\10Hz\615883', '\615883\03122021\10Hz\615883'};
    else
        % 140Hz
        allSES = {'\608448\01102020\608448', '\608451\01102020\608451', '\608452\01102020\608452', '\602101\03112020\602101', ...
            '\611311\03112020\611311', '\611111\03112020\611111', '\608450\01172020\608450', '\608451\07012020\608451', ...
            '\608452\07012020\608452', '\611111\07012020\611111', '\612535\07022020\612535', '\602101\07022020\602101', ...
            '\615883\02052021\145Hz\615883', '\615883\03122021\145Hz\615883'};
    end

    if stim_freq == 1
        ses_sel = 2:14;
    else
        ses_sel = 1:11;
        ses_sel = [ses_sel, 13, 14];
    end

    for ses = ses_sel
        ses
        if stim_freq == 1
            cd([datapath allSES{ses} '_10Hz_AudioVisual\motion_corrected']);
        else
            cd([datapath allSES{ses} '_145Hz_AudioVisual\motion_corrected']);
        end

        % Load motion and trace data
        motion = h5read('processed_motion.h5', '/raw_speed_trace');
        tracesF = h5read('processed_trace.h5', '/trace');
        traces = h5read('processed_trace.h5', '/onset_binary_trace');
        load('LFP_ts.mat')
        stim_onsets = LFP_data.Stim_onset;
        stim_offsets = LFP_data.Stim_offset;

        % Preprocess motion data
        v = motion;
        idx = find(~isnan(v));
        v = v(~isnan(v));
        traces = traces(idx, :);
        tracesF = tracesF(idx, :);


        Sampling_freq=20;
        % Mean across all neuron
        mean_traces=zscore(nanmean(traces,2));  %% Here is use median instead of mean as a better population average
        % High-speed periods
        moving_period=h5read('processed_motion.h5','/moving_period');
        moving_period=moving_period(idx);  % moving period in 0 and 1

        % Define constants and parameters
        Sampling_freq = 20;
        time_vect1 = 0:1/Sampling_freq:2099.9;
        signal1 = v;
        time_vect2 = 0:1/1000:2099.9;
        stim_vec = zeros(1, size(tracesF, 1));

        % Create stimulus vectors
        for i = 1:length(stim_onsets)
            timsel = (stim_onsets(i) - 0:stim_onsets(i) + 1200 - 1) - idx(1);
            stim_vec(timsel) = 1;
        end

        % Load and preprocess LFP data
        stim_vecLFP = zeros(1, size(LFP_data.LFP, 1));
        for i = 1:length(LFP_data.LFP_stim_onset)
            timsel = (LFP_data.LFP_stim_onset(i) - 0:LFP_data.LFP_stim_offset(i));
            stim_vecLFP(timsel) = 1;
        end
        LFP = LFP_data.LFP;
        vv = LFP(1:50:end);
        vv = vv(idx);

        % Align LFP and motion
        start_frame = LFP_data.Start_Imaging;
        delay_frames = start_frame + idx(1);
        delay_frame_LFP = ceil(delay_frames * 1000 / 20);
        shifted_stim_onsets = (stim_onsets - idx(1)) / 20 * 1000;
        shifted_stim_offsets = (stim_offsets - idx(1)) / 20 * 1000;
        aligned_LFP = LFP(delay_frame_LFP:end);
        stim_vecLFP = stim_vecLFP(delay_frame_LFP:end);

        % Interpolate motion data
        v_Intp = interp1(time_vect1, v(1:length(time_vect1)), time_vect2);
        moving_period = interp1(time_vect1, moving_period(1:length(time_vect1)), time_vect2);
        aligned_LFP = aligned_LFP(1:length(v_Intp));
        stim_vec = interp1(time_vect1, stim_vec(1:length(time_vect1)), time_vect2);
        v_Intp = v_Intp(1:5:end);
        moving_period = moving_period(1:5:end);
        aligned_LFP = aligned_LFP(1:5:end);
        stim_vec = stim_vec(1:5:end);
        stim_vecLFP = stim_vecLFP(1:5:end);

        % Perform further processing
        if 1
            % Set sampling frequency and filter parameters
            FS = 200;
            Fn = FS / 2;
            FB = [7 8];
            [B, A] = butter(2, [min(FB) / Fn max(FB) / Fn]);

            % Apply bandpass filter to aligned_LFP
            thetaS = filtfilt(B, A, aligned_LFP);

            % Z-score normalize aligned_LFP
            aligned_LFP = zscore(aligned_LFP);

            % Calculate deviations
            devianz = abs(aligned_LFP) > 5;

            % Initialize lfp_all structure
            lfp_all = struct();
            lfp_all.trial{1}(1, :) = aligned_LFP;
            lfp_all.trial{1}(2, :) = v_Intp;
            lfp_all.time{1} = (1:size(lfp_all.trial{1}, 2)) / FS;
            lfp_all.label = {'LFP', 'motion'};

            % Filter and process motion data (deltP1)
             deltP1=lfp_all.trial{1}(2,:);
            FB = [3 4];
            [B, A] = butter(2, [min(FB) / Fn max(FB) / Fn]);

            % Calculate delta phase (deltPh) and delta power (deltP)
            deltPh = angle(hilbert(filtfilt(B, A, deltP1)));
            deltP = fastsmooth(abs(hilbert(filtfilt(B, A, deltP1))), 100, 1, 1);

            % Set threshold for delta power
            thresD = prctile(deltP(moving_period == 1), 0);

            % Detect delta triggers (dtrigs)
            dtrigs = zeros(1, length(deltPh));
            timer1 = 0;

            for xi = 1:length(deltPh)
                timer1 = timer1 + 1;

                if deltPh(xi) > 0 && deltPh(xi) < 0.2 && timer1 > 40
                    dtrigs(xi) = 1;
                    timer1 = 0;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Configure parameters for frequency analysis using wavelet transform
            cfg = struct();
            cfg.method = 'wavelet';       % Method for frequency analysis
            cfg.output = 'fourier';       % Output type
            cfg.taper = 'hanning';        % Tapering function
            cfg.keeptapers = 'yes';       % Keep tapers
            cfg.keeptrials = 'yes';       % Keep trials
            cfg.trials = 'all';           % Use all trials
            cfg.tapsmofrq = 6;            % Taper smoothing frequency
            cfg.channel = 'all';          % All channels
            cfg.foi = 8:2:100;            % Frequencies of interest
            cfg.t_ftimwin = ones(1, length(cfg.foi)) * 0.15;  % Time window for each frequency
            cfg.toi = lfp_all.time{1};     % Time of interest
            cfg.width = 6;                % Width parameter

            % Perform frequency analysis using wavelet transform
            freq2 = ft_freqanalysis(cfg, lfp_all);

            % Extract the Fourier coefficients (complex values) for each channel and frequency
            cwt_out = freq2.fourierspctrm;

            % Define the channel of interest (cha)
            cha = 1;
            wavA = abs(squeeze(cwt_out(1, cha, :, :)));      % Amplitude
            % Thresholding 
              
            wavA(:, ~(moving_period   )) = NaN;  %
           
            % Define window parameters
            lwin = 100;
            rwin = 100;

            % Calculate theta-triggered waveM
            [data] = theta_trig_waveM(wavA, dtrigs, lwin, rwin);

            % Calculate theta-triggered traceCFC
            [dataTr] = theta_trig_traceCFC(deltP1, dtrigs, lwin, rwin);


            mm = mm + 1;

            % Calculate A and C matrices for baseline and stimulation
            A = squeeze(nanmean(data, 3)); % Baseline
            C1 = bsxfun(@minus, A, nanmean(A, 2));
            C2 = bsxfun(@plus, A, nanmean(A, 2));
            aMAP(:, :,  mm) = smooth2a(C1 ./ C2, 1, 1);

          

            % Store additional results
            Atrace = [Atrace, nanmean(dataTr, 2)];



        end
    end
end

% Save results
%savepath = 'C:\Users\hanlab\Dropbox (BOSTON UNIVERSITY)\Paper - Sudi Striatum Audio Visual Stim\Eric_figures\plots';


savepath='C:\Users\sudiksha\Documents\Codes\Spectral_analysis\Figures\';
figure('COlor','w','Position', [ 300 400 220 180],'Renderer', 'painters')
imagesc((-lwin:rwin)./200,freq2.freq,smooth2a(nanmean(aMAP(:,:,:),3),1,5)) ;  %Stim
colormap(jet)
set(gca,'Clim',[-0.004 .004])
colorbar
%set(gca,'Clim',[0.987 1.012])
axis xy
hold on,plot( (-lwin:rwin)./200 ,((zscore(nanmean(Atrace,2)).*10)+40),'k')
xlim([-0.3 0.3])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'CFC_all_Pow_stim' num2str(stim_freq) '.pdf'])



