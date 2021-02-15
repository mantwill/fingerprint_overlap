
clc
clearvars


% This script uses functional connectivity to predict behaviour from our preprint:
% https://www.biorxiv.org/content/10.1101/2021.02.07.429922v1

% Input connectivity matrices used for prediction have been vecotrised after preprocessing
% and need to be in the format of: atlas_subjectID_session (session being either rest1 or rest2)

% This script uses the Connectome-based modelling approach developed by
% Finn et al. 2015 (https://www.nature.com/articles/nn.4135) which we addapted from
% the openly available script found here: https://www.nitrc.org/frs/?group_id=51

% Reference:
% Finn, E. S., Shen, X., Scheinost, D., Rosenberg, M. D., Huang, J., Chun, M. M., ... & Constable, R. T. (2015). 
% Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity. 
% Nature neuroscience, 18(11), 1664-1671.

% Martin Gell and Maron Mantwill 15/02/2021



%%%% SETUP %%%%
% PARAMETERS
atlasname = 'atlasname'; % Name of atlas
n_node = 268; % Number of ROIs in atlas
thresholds = [0.05]; % Threshold
n_fold = 10; % Number of folds for CV
n_repeat = 100; % Number of repetitions for CV

% Directory Setup
inputdir = fullfile(pwd(),'fc/'); % folder with vecotrised fc matrices
beh_path = fullfile(pwd(),'behavioral_data/Behaviour.csv'); % individual behavioural scores
outputdir = fullfile(pwd(), 'results/prediction_results/');
functiondir = fullfile(pwd(), 'functions/'); % Function for reasembeling triangles to matrices
addpath(genpath(functiondir));
% ----------------------------------------



% load behavioural data
b = readtable(beh_path);

%% PREPARE REST 1 RS FMRI FOR PREDICTION DATA
% data has to be a matrix with dimensions ROI by ROI

% upper triangle indices for reasembeling vectors to matrices
FC_mat = zeros(n_node);
upper = find(triu(ones(length(FC_mat(1,:))),1));

cd(inputdir);
% Rest_1_mat
% Load subj IDs
sample_subs = table2array(b(:,1));
rest_1_mats = zeros(length(FC_mat(1,:)),length(FC_mat(1,:)),length(sample_subs));

for i = 1:length(sample_subs)
    sub_i = sample_subs(i);
    fc_vec = csvread(strjoin({atlasname, num2str(sub_i), 'rest1.csv'}, '_')); % read in functional connectivity 
    fc = reasemble_FC_mat(upper, fc_vec, length(FC_mat(1,:)), 1); % this step can be skipped if fuctional connectivity is already matrix
    rest_1_mats(:,:,i) = fc;
end
clear fc fc_vec file_i


% Behavioural prediction
% Loop over all thresholds
for threshold_i = 1:length(thresholds)
    
    %%%% SETUP %%%%
    subs = b(:,1);
    behaviour = b(:,2:end);             % define columns with behavioral values
    thresh = thresholds(threshold_i);   % threashold for feature selection in current loop itteration
    % ----------------------------------------
    
    % varaible to save all results into
    all_results = cell(size(behaviour,2), 10);
    
    %%%% CROSS VALIDATION %%%%
    % Loop through behavioural measures
    for beh_i = 1:size(behaviour,2)

        %%%% SETUP %%%%
        beh = [subs behaviour(:,beh_i)];      % table of behaviour to be predicted
        beh = standardizeMissing(beh,-99);
        beh = rmmissing(beh);                 % remove missing beh. scores
        beh_name = beh.Properties.VariableNames(2);
        beh_name = beh_name{1};               % name of behaviour to be predicted
        n_sub = size(beh, 1);                 % number of subjects
        % ----------------------------------------
        
        fprintf('Working on >> %s that is: %i of %i \n\n', beh_name, beh_i, size(behaviour,2));
        
        % Extract behaviour to be predicted from table
        curr_beh = table2array(beh(:,2));
        
        % Extract data to be used for prediction from all rest1 mats
        subsample_subs = ismember(sample_subs,table2array(beh(:,1)));
        all_mats = rest_1_mats(:,:,subsample_subs);
        
        % Check if subjects used
        if ~all(sample_subs(subsample_subs) == table2array(beh(:,1)))
            error('Subjects arent the same!')
        end
        
        % Variables for saving results - feature selection
        selected_features_pos = zeros(length(FC_mat));
        selected_features_neg = zeros(length(FC_mat));
        
        % Variables for saving results - accuracy
        predictions_pos = zeros(n_sub,n_repeat);
        predictions_neg = zeros(n_sub,n_repeat);
        
        % loop through M repeats
        for repeat = 1:n_repeat
            
            % Create cv indices
            indices = crossvalind('Kfold',n_sub,n_fold);
            
            fprintf('Working on >> rep %i of %i \n', repeat, n_repeat);
            
            % loop through K folds
            for fold = 1:n_fold
                
                %fprintf('test fold: %i \n', fold);
                
                % Identify train and test subjects
                train_ind = indices ~= fold;
                test_ind = ~train_ind;
                train_subs = all_mats(:,:,train_ind);
                test_subs =  all_mats(:,:,test_ind);
                
                % Identify train and test behaviour
                train_beh = curr_beh(train_ind);
                test_beh = curr_beh(test_ind);
                
                %%%% PREDICTION --- CPM %%%%
                
                % Variables to save positive and negative model
                behav_pred_pos = zeros(size(test_subs,3),1);
                behav_pred_neg = zeros(size(test_subs,3),1);
                
                train_vcts = reshape(train_subs,[],size(train_subs,3));
                
                % correlate all edges with behavior
                [r_mat,p_mat] = corr(train_vcts',train_beh);
                
                r_mat = reshape(r_mat,n_node,n_node);
                p_mat = reshape(p_mat,n_node,n_node);
                
                % set threshold and define masks
                pos_mask = zeros(n_node,n_node);
                neg_mask = zeros(n_node,n_node);
                
                pos_edges = find(r_mat > 0 & p_mat < thresh);
                neg_edges = find(r_mat < 0 & p_mat < thresh);
                                
                pos_mask(pos_edges) = 1;
                neg_mask(neg_edges) = 1;
                
                % Update concencus maps
                selected_features_pos = selected_features_pos + pos_mask;
                selected_features_neg = selected_features_neg + neg_mask;
                
                % get sum of all edges in TRAIN subs (divide by 2 to control for the
                % fact that matrices are symmetric). That is, sum up edges that pass
                % the threashold in each subject separately.
                train_sumpos = zeros(size(train_subs,3),1);
                train_sumneg = zeros(size(train_subs,3),1);
                
                for ss = 1:size(train_sumpos)
                    train_sumpos(ss) = sum(sum(train_subs(:,:,ss).*pos_mask))/2;
                    train_sumneg(ss) = sum(sum(train_subs(:,:,ss).*neg_mask))/2;
                end
                
                % build model on TRAIN subs
                fit_pos = polyfit(train_sumpos, train_beh,1);
                fit_neg = polyfit(train_sumneg, train_beh,1);
                                
                % run model on TEST subs
                test_sumpos_plot = zeros(size(test_subs,3),1);
                for test_sub_i = 1:size(test_subs,3)
                    test_sumpos = sum(sum(test_subs(:,:,test_sub_i).*pos_mask))/2;
                    test_sumneg = sum(sum(test_subs(:,:,test_sub_i).*neg_mask))/2;
                    
                    test_sumpos_plot(test_sub_i) = test_sumpos;
                    behav_pred_pos(test_sub_i) = fit_pos(1)*test_sumpos + fit_pos(2);
                    behav_pred_neg(test_sub_i) = fit_neg(1)*test_sumneg + fit_neg(2);
                end
                % collect results better
                predictions_pos(test_ind,repeat) = behav_pred_pos;
                predictions_neg(test_ind,repeat) = behav_pred_neg;
                
            end
        end
        
        % save average predictions of models
        mean_behav_pred_pos = mean(predictions_pos,2);
        mean_behav_pred_neg = mean(predictions_neg,2);
        
        % r value for analysis
        [c_pos,p_pos] = corr(mean_behav_pred_pos,curr_beh);
        [c_neg,p_neg] = corr(mean_behav_pred_neg,curr_beh);
        
        mse_pos = sum((mean_behav_pred_pos - curr_beh).^2)*(1/n_sub);
        mse_neg = sum((mean_behav_pred_neg - curr_beh).^2)*(1/n_sub);
        
        all_results(beh_i,:) = [{beh_name} c_pos p_pos mse_pos length(pos_edges) c_neg p_neg mse_neg length(neg_edges) n_sub];
        
        csvwrite(strjoin({outputdir,'average_predictions/', beh_name, '_', num2str(thresh), '_pos.csv'}, ''),[mean_behav_pred_pos curr_beh]);
        csvwrite(strjoin({outputdir,'average_predictions/', beh_name, '_', num2str(thresh), '_neg.csv'}, ''),[mean_behav_pred_neg curr_beh]);
        
        % Save results for each CV
        csvwrite(strjoin({outputdir,'selected_edges/', beh_name, '_', num2str(thresh), '_pos.csv'}, ''),selected_features_pos);
        csvwrite(strjoin({outputdir,'selected_edges/', beh_name, '_', num2str(thresh), '_neg.csv'}, ''),selected_features_neg);
        
    end
    
    % Convert cell to a table and use first row as variable names
    T = cell2table(all_results,'VariableNames',{'Behaviour' 'c_pos' 'p_pos' 'mse_pos' 'n_pos_edges' 'c_neg' 'p_neg' 'mse_neg' 'n_neg_edges' 'N'});
    
    % Write the table to a CSV file
    writetable(T,strjoin({outputdir,'overall_accuracies/', 'all_behaviours', '_', num2str(thresh), '_results.csv'},''));
    
end

