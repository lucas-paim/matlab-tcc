% 1) Set-up workspace
% 2) EEGLAB section
% 3) FIELDTRIP section
% 4) Count the number of trials
% 5) Functions

clear;
close all;

%% 1) Set-up Workspace

restoredefaultpath

cd('C:\Users\lucas\source\matlab_eeg\dados_lucasp')
cwd = [cd filesep]; 

wpms = []; % pre-allocate a workspace variables in struct

wpms.DATAIN     = fullfile(cwd, 'Datain', filesep); 
wpms.DATAOUT    = [cwd 'Dataout' filesep];
wpms.FUNCTIONS  = fullfile(cwd, 'Functions', filesep);
wpms.software   = 'C:\Users\lucas\source\matlab_eeg\dados_lucasp\Functions\';

load([wpms.FUNCTIONS filesep 'chanlocs.mat'])
load([wpms.FUNCTIONS 'neighbour_template.mat'])
neighbours = neighboursCopy;

addpath([wpms.software 'eeglab_current\eeglab2024.2'])
addpath(genpath([wpms.software 'eeglab_current\eeglab2024.2\functions']));
addpath([wpms.software 'eeglab_current\eeglab2024.2\plugins\xdfimport1.19']);
addpath([wpms.software 'eeglab_current\eeglab2024.2\plugins\bva-io']);
%addpath([wpms.software 'eeglab_current\eeglab2024.2\plugins']); % add all plugins 
addpath(genpath([wpms.software 'eeglab_current\eeglab2024.2\plugins\PrepPipeline0.56.0']));
addpath(genpath([wpms.software 'eeglab_current\eeglab2024.2\plugins\firfilt']));

addpath([wpms.software 'fieldtrip-20241025']);
addpath([wpms.software 'fieldtrip-20241025\external\eeglab']);

% Structure and store participant (px) codes
% pxlist = dir([wpms.DATAIN 'px-*']); 
% test only on px-n
pxlist = dir([wpms.DATAIN 'px-08*']); 

mkdir(wpms.DATAOUT); % Create Output directory
for subout = 1:length(pxlist)
    mkdir([wpms.DATAOUT pxlist(subout).name]) % Create folder for each px
end

%% 2) EEGLAB Section
eeglab
close all

skipped = 1; % contador para participantes cujos dados não serão processados corretamente.

for px =1:length(pxlist)
    clearvars all_eeg_files EEG n_triggers triggers rejected
    fprintf(['\n Analysing participant: ' pxlist(px).name '\n\n']);
    all_eeg_files  = dir([wpms.DATAIN pxlist(px).name filesep 'EEG' filesep 'Raw' filesep '*eeg_php*.vhdr']);
    
    % another loop within initial px loop, loops through each resting state file for the current px.
    for file = 1:length(all_eeg_files)
    
        clearvars EEG end_of_window end_time n_triggers start_of_window triggers rejected filename
        filename = get_filename(all_eeg_files(file).name);
        
        %fprintf(['\n Analysing session: ' all_eeg_files(file).name '\n\n']);
        fprintf(['\n Analysing session: ' filename '\n\n']);
       
        % Loads the raw EEG file
        [EEG, ~] = pop_loadbv(all_eeg_files(file).folder, [all_eeg_files(file).name]); % Do you need to have opened and cleared eeglab?
        
        % store all the triggers/events in a variable
        for i = 1:length(EEG.event)
            triggers{i} = EEG.event(1, i).type; 
        end
        
         triggers = string(triggers);

         %triggers are markers in the EEG file that mark when something happened. In this case each marker means:
         % R5 - heat beginning to ramp up
         % R4 - heat reached 46 degrees (pain)
         % R6 - heat begining to ramp down
         % R15 - heat reached 32 degrees (neutral temp - no pain)
         % So at this point in the script you would need to check whether it have 5*each of the different triggers.

        % pause;
        % ask for input depending
        % if 1 entered, mark file as being processed and continue
        % if 0 entered, mark file as not being processed 
        prompt = 'Please indicate whether file has 4 codes indicating start and stop of resting states (1 = yes, 0 = no):'; 

        process = input(prompt); % prompt saved in this variable
        processed{px, file} = process;
        
        % n_triggers = length(triggers(triggers == 'Start/End'));
        n_triggers = length(triggers); % count how many triggers

        if process == 1

            % downsample to 500
            EEG = pop_resample(EEG, 500);

            EEG = pop_select(EEG, 'nochannel', {'GSR','HR','RESP'});
            %EEG = pop_reref(EEG, []); % common average % skip till after the removing channels
            EEG = pop_eegfiltnew(EEG, 2, 100, 826, 0, [], 1);
            close all;

            pop_eegplot( EEG, 1, 1, 1); %Check if the triggers are correct.

            if sum(count(triggers, "EC")) == 1
                % all good
            %elseif sum(count(triggers,"EC")) == 2
                %c = find( triggers == "EC", 2);
                %c = c(1);
                %EEG.event(1, c).type = 'EO';
                %c = find( triggers == "Start/End", 2);
                %c = c(1);
                %EEG.event(1, c).type = 'EC';
            
            else
                eo_count = sum(count(triggers, "EO Start/End"));
                se_count = sum(count(triggers, "Start/End"));
                if eo_count == 2
                    EEG.event(1, find( triggers == "Start/End", 1)).type = 'EC'; 
                elseif eo_count == 0 && se_count == 4
                    c = find( triggers == "Start/End", 3);
                    c = c(3);
                    EEG.event(1, c).type = 'EC';
                end
            end
            
            %if n_triggers == 2;
            %    EEG = pop_rmdat(EEG, {'Start'},[3 303], 0);
            %else
            %    EEG.event(1, (length(EEG.event)+1)).type = 'End';
            %    EEG.event(1, end).latency = length(EEG.times);
            %    pop_eegplot(EEG, 1, 1, 1);
            %    prompt = 'Please identify the end time';
            %    end_time = input(prompt);
            %    end_of_window = end_time - length(EEG.times)/500;
            %    start_of_window = end_of_window - 300;
            %    EEG = pop_rmdat(EEG, {'End'},[start_of_window end_of_window], 0);
            %end
            
            % This cuts the file based around the trigger we want from 3 seconds after to 183 seconds after the eyes closed trigger
            EEG = pop_rmdat(EEG, {'EC'},[3 183], 0);
            %EEG = pop_rmdat(EEG, {'EC'},[0 210], 0);
            %EEG = pop_rmdat(EEG, {'R  4'},[0 40], 0);
            %pop_eegplot( EEG, 1, 1, 1);

            close all 

            % rejecting bad channels/electrodes/sensors
            % opens up two figures to view the data
            figure; pop_spectopo(EEG, 1, [0 15000], 'EEG', 'freqrange', [2 100], 'electrodes','off');
            pop_eegplot(EEG, 1, 1, 1);
            prompt = 'Please identify sensor # to remove in matrix format (e.g., [12,14]):'; 

            rejected = input(prompt); % prompt saved in this variable
            rejected_channels{px, file} = rejected; % saves all rejected channels 

            % if statement 
            % if the number in rejected is not zero, then remove the bad channel storred in the variable 'rejected'
            if size(rejected,1)~=0
                EEG = pop_select(EEG, 'nochannel',rejected);  
            end

            % re-reference to the common average of all electrodes
            EEG = pop_reref(EEG, []);

            % pops up another figure to double check what the spectra look like
            figure; pop_spectopo(EEG, 1, [0 15000], 'EEG', 'freqrange', [2 100], 'electrodes','off');

            % Saves the data 
            save([wpms.DATAOUT pxlist(px).name filesep pxlist(px).name '_' filename '_preICA_maxclean.mat'], 'EEG');
            fprintf(['\n Finished analysing session: ' filename '\n\n']);
                
            
            %close all
        else
            processed_no{skipped, 1} = all_eeg_files(file).name;
            processed_no{skipped, 2} = process; % save whether 1 or 0 for sanity check
            skipped = skipped + 1;
            fprintf(['\n File: ' filename ' did not contain correct resting state recording triggers, so was skipped \n\n']);
        end
    end
    fprintf(['\n Finished analysing participant: ' pxlist(px).name '\n\n']);
end

save([wpms.DATAOUT 'rejected_channels.mat'], 'rejected_channels');
save([wpms.DATAOUT 'not_processed_list.mat'], 'processed_no');
save([wpms.DATAOUT 'processed_codes.mat'], 'processed');

%% 3) FIELDTRIP Section

clearvars all_eeg_files EEG n_triggers prompt rejected subjlist triggers px file
pxlist = dir([wpms.DATAOUT 'px-*']);
not_2min_count = 1; 

for px = length(pxlist)

    clearvars all_eeg_files cfg data data_rejected data_freq data_pruned EEG freq paf power promt n_rejected rejected temp temp2 X y z
    fprintf(['\n Analysing participant: ' pxlist(px).name '\n\n']);

    all_eeg_files  = dir([wpms.DATAOUT pxlist(px).name filesep '*_preICA_maxclean.mat']);

    % loop over each file for current px
    for file = 1:length(all_eeg_files)
        clearvars cfg data data_freq data_pruned EEG freq paf power promt rejected temp temp2 X y z n_rejected filename
        filename = get_filename_ICA(all_eeg_files(file).name);

        load([all_eeg_files(file).folder filesep all_eeg_files(file).name], 'EEG');
        fprintf(['\n Analysing participant and session: ' pxlist(px).name ' ' filename '\n\n']);

        % EEG is the data structure for eeglab, and here we convert EEG using the function 'eeglab2fieldtrip' to a different structre that we name data
        data = eeglab2fieldtrip(EEG, 'preprocessing');
        data.label = {EEG.chanlocs.labels};

        % cut the signal we have into epochs of 5 seconds
        cfg = []; % set up empty configuration 
        cfg.length = 5;
        cfg.overlap = 0;
        data=ft_redefinetrial(cfg,data);

        chanSel = 1:63; % number of channels we have
        cfg = [];
        cfg.channel = 'all';
        cfg.ylim = [0 20]; % auto set scale, but can be changed manually in window once it pops up
        cfg.blocksize = 6;
        cfg.viewmode = 'vertical';

        % opens the figure to browse data
        visualArtifacts = ft_databrowser(cfg,data);

        % Read out the trials that contain artifacts
        artifacts = visualArtifacts.artfctdef.visual.artifact;

        a_trials = []; % trials with artifacts will be numbered here
        for n=1:size(artifacts,1)
            a_trials = [a_trials;find(data.sampleinfo(:,1) <= artifacts(n,1) ...
                & data.sampleinfo(:,2) >= artifacts(n,1))]; % concatenate the trials that contain artifacts
        end

        % find out which trials are artefact free
        trials = 1:numel(data.trial); % list all the trials
        trials_visual = find(~ismember(trials,a_trials))';  % trials without visual artifacts

        data.trialinfo = trials'; 
        % the ' makes it long format (i.e. a column rather than a row)

        % add column with artifact detection index, 1 is clean, 0 is artifact
        noVisual = ismember(1:length(trials),trials_visual)';
        data.trialinfo = [data.trialinfo noVisual];

        % reject artifacts that were manually highlighted
        cfg = [];
        cfg.artfctdef.reject = 'complete'; % default to delete complete trial
        cfg.artfctdef.visual.artifact = visualArtifacts.artfctdef.visual.artifact;
        data_rejected = ft_rejectartifact(cfg, data);

        % store the number of epoched rejected
        n_rejected = length(data.trial) - length(data_rejected.trial);

        rest_num = 1;

        %n_rejected = length(data.trial) - length(data_rejected.trial);
        n_clean_epochs{px, rest_num} = length(data_rejected.trialinfo);

        if n_rejected > 12
            not_2mins{not_2min_count, 1} = all_eeg_files(file).name;
            not_2mins{not_2min_count, 2} = n_rejected; % save whether 1 or 0 for sanity check
            not_2min_count = not_2min_count + 1;
        end


        %Run ICA with 15 components to observe
        cfg = [];
        cfg.method = 'runica'; % the standard method
        cfg.runica.pca = 15; % pca - which is slightly different to ica
        % data struct is temporarily saved as a variable called 'X' here
        X = ft_componentanalysis(cfg,data_rejected);

        %Create plots to look at the components
        cfg.component = [1:15];
        cfg.layout = 'EEG1010.lay';
        ft_topoplotIC(cfg,X); % topoplots (i.e. heads)
        ft_databrowser(cfg,X); % browse the signal across time

        %Reject Visual Components
        prompt = 'Please Identify Component ro reject in matrix form, (i.e. [1,4]) :';
        rejected = input(prompt); % stores components to reject
        rejected_components{px, file} = rejected;

        cfg.component = rejected;
        data_postICA = ft_rejectcomponent(cfg,X); % reject components 

        %Fourier Transform with Hanning taper
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.taper = 'hanning';
        cfg.foi = 2:.20:50; %freq of interest 2-50Hz with 0.2Hz bins
        cfg.keeptrials = 'yes';
        data_freq = ft_freqanalysis(cfg, data_postICA);
        close all

        save([wpms.DATAOUT pxlist(px).name filesep pxlist(px).name '_' filename '_postICA_maxclean.mat'],'data','data_freq', 'data_postICA','data_rejected', 'n_rejected');% Nahian - Adjust based on what you th
        fprintf(['\n Finishined analysing participant and timepoint: ' pxlist(px).name ' ' filename '\n\n']);

EEG = fieldtrip2eeglab(data_postICA);
    end
    fprintf(['\n Finished analysing participant: ' pxlist(px).name '\n\n']);

    save([wpms.DATAOUT 'n_clean_epochs.mat'], 'n_clean_epochs');
    save([wpms.DATAOUT 'rejected_components.mat'], 'rejected_components');

end
save([wpms.DATAOUT 'n_clean_epochs.mat'], 'n_clean_epochs');
save([wpms.DATAOUT 'rejected_components.mat'], 'rejected_components');
fprintf('\n Finished ICA part of the script. \n\n');

%% 4) Count the number of trials

clearvars all_eeg_files EEG n_triggers prompt rejected subjlist triggers px file data_postICA file filename px n_clean_epochs n_clean_epochs_check not_2mins data_freq artifacts
clearvars not_2min_count not_2mins_check

pxlist = dir([wpms.DATAOUT 'px-*']); % list of all px in output folder
not_2min_count = 1;

for px = 1:length(pxlist) % loop through participant list
    all_eeg_files  = dir([wpms.DATAOUT pxlist(px).name filesep '*postICA_maxclean.mat']);
    
    % loop through each file for the current participant
    for file = 1:length(all_eeg_files)
        % load file saved after fourier transform in fieldtrip
        load([all_eeg_files(file).folder filesep all_eeg_files(file).name], 'data_freq');

        rest_num = get_rest_num(all_eeg_files(file).name);
        if rest_num == 11
            rest_num = 1;
        end
        filename = get_filename_ICA(all_eeg_files(file).name);
        fprintf(['\n Analysing participant and session: ' pxlist(px).name ' ' filename '\n\n']);

        %n_rejected = length(data.trial) - length(data_rejected.trial);
        n_clean_epochs_check{px, rest_num} = length(data_freq.trialinfo);

        if n_clean_epochs_check{px, rest_num} < 24
            not_2mins_check{not_2min_count, 1} = all_eeg_files(rest_num).name;
            not_2mins_check{not_2min_count, 2} = n_clean_epochs_check{px, rest_num}; % save whether 1 or 0 for sanity check
            not_2min_count = not_2min_count + 1;
        end
        fprintf(['\n Finishined analysing participant and timepoint: ' pxlist(px).name ' ' filename '\n\n']);
    end
    fprintf(['\n Finished analysing participant: ' pxlist(px).name '\n\n']);
    save([wpms.DATAOUT 'n_clean_epochs_check.mat'], 'n_clean_epochs_check');


end
save([wpms.DATAOUT 'not_2mins_check_02052022.mat'], 'not_2mins_check');

not_2mins_check;
%pause;

fprintf('\n Finished check of participants requiring manual restart to get more epochs. \n\n');
data = EEG.data;

    
%% 5) Functions

function filename = get_filename(full_filename)
    filename = split(full_filename, '_');
    filename = filename(length(filename)); 
    filename = split(filename, '.'); 
    filename = filename{1};
end

function filename = get_filename_ICA(full_filename)
    filename = split(full_filename, '_');
    filename = filename{2}; 
    %filename = split(filename, '.'); 
    %filename = filename{1};
end

function rest_num = get_rest_num(full_filename)
    rest_num = split(full_filename, '_');
    rest_num = split(rest_num(2), 'rest'); 
    rest_num = rest_num{2};
    rest_num = str2double(rest_num);
end
