clc
clear

% This script calcualtes the edgewise overlap of edges used for prediction of 
% behaviour and discriminatory edges used in fingerprinitng from our preprint:
% https://www.biorxiv.org/content/10.1101/2021.02.07.429922v1 

% This script uses the function randomizer_bin_und from 
% Brain Connectivty Toolbox (https://sites.google.com/site/bctnet/). 
% Before running, please make sure to add the toolbox to your Matlab path.

% Martin Gell and Maron Mantwill 15/02/2021



%%%% SETUP %%%%
% Threshold for number of edges to pass as 'predictive'
% With 100*10 fold cv, maximum consensus will be 1000, i.e. edges that were selected in every fold of the cv.
% Here we use 800 which means only edges selected on 80% of cv folds are used 
consensus = 800;
n_perm = 100; % Number of permutations for degree preserving random matrices

% File name of DP edge threshold to use
filename_discriminatoryedges = 'example_DP_edges_0.01.txt';

% Directory Setup
discriminatory_dir = fullfile(pwd(), 'results/discriminatory_edges/');
prediction_dir = fullfile(pwd(), 'results/prediction_results/selected_edges/');
outputdir = fullfile(pwd(), 'results/overlap_results/');
functiondir = fullfile(pwd(), 'functions/');
addpath(functiondir);
% ----------------------------------------



% Load discrim edges
fing_edges = dlmread(strjoin({discriminatory_dir, filename_discriminatoryedges},''));

% Prediction weights
cd(prediction_dir);
behs = dir('*.csv');

% Save the number of overlaps with each random network here
overlaps = zeros(n_perm,1);

% Save results here
res = cell(length(behs), 7);

% start counter
i = 1;


% Loop over all prediction edges
for beh_i = behs'
        
    % Load in prediction weights for beh_i
    weights = dlmread(strjoin({beh_i.folder beh_i.name}, '/'));
    
    % Binarize prediction weights
    bin = weights;
    bin(bin < consensus)  = 0;
    bin(bin >= consensus) = 1;
    pred_edges = bin;
    n_pred_edges = (sum(sum(pred_edges))/2); % save n of bin edges
    
    % Save binarized weights as prediction edges
    beh_name = beh_i.name; % name of the behaviour in question
    dlmwrite(strjoin({outputdir  'binarized_matrices/binarized_prediction_edges_' num2str(consensus) '_' beh_name '.txt'},''), pred_edges, 'delimiter', ' ');

    % Perform permutation of overlap with degree preserving random networks
    for perm_i = 1:n_perm
        
        % Create a degree preserving random graph
        rand_edges = randomizer_bin_und(pred_edges,1);
        
        % Calculate overlap
        overlap = rand_edges + fing_edges;
        overlap(overlap == 1) = 0;
        overlap(overlap == 2) = 1;
        
        % Save n of overlap
        n_overlap = sum(overlap(overlap == 1))/2;
        overlaps(perm_i) = n_overlap;

        clear overlap
    end
    
    % plot histogram of all permutation overlaps
    figure(1);
    histogram(overlaps,length(unique(overlaps)));
    xticks(unique(overlaps));

    % Calculate empirical overlap wiht prediction edges
    overlap = pred_edges + fing_edges;
    overlap(overlap == 1) = 0;
    overlap(overlap == 2) = 1;
    n_overlap_emp = sum(overlap(overlap == 1))/2;
    
    % calculate p-value
    p = length(overlaps(overlaps >= n_overlap_emp))/n_perm;
    
    % Save
    n_fing_edges = sum(sum(fing_edges))/2;
    res(i,:) = [{beh_name} n_fing_edges n_pred_edges n_overlap_emp mean(overlaps) std(overlaps) p];
            
    % Save matrix of overlapping edges
    dlmwrite(strjoin({outputdir  'binarized_matrices/binarized_empirical_overlap_' num2str(consensus) '_' beh_name '.txt'},''), overlap, 'delimiter', ' ');

    % Update counter
    i = i + 1;
end


% Save results
T = cell2table(res, 'VariableNames', {'Behaviour' 'N_fing_edges' 'N_predictive_edges' 'Overlap' 'Mean_permuted_overlap' 'SD_permuted_overlap' 'p_permuted'});
writetable(T,strjoin({outputdir,  'permutation_results_allbehaviours_bin_' num2str(consensus) '.csv'},''));
