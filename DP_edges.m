
% This script calcualtes the discriminatory potential of edges used in our preprint:
% https://www.biorxiv.org/content/10.1101/2021.02.07.429922v1 

% It is closely adapted from Shen et al. mentioned below.

% Martin Gell and Maron Mantwill 15/02/2021


% Copyright 2015 Xilin Shen and Corey Horien

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% The individual functional connectome is unique and stable over months to years
% Corey Horien, Xilin Shen, Dustin Scheinost, R. Todd Constable, bioRxiv
% doi: https://doi.org/10.1101/238113 

%this code allows the calculation of which edges are helpful in ID (the
%differential power (DP) metric). For a given edge, the likelihood that
%within-subject similarity is higher than that between subjects, is
%calculated. This is done by comparing the product of edge values from
%time1 and time2 from the same subject to the product of time1 and time 2
%from unmatched subjects. Edges with high DP are helpful in ID.

clc;
clear;


%%%% SETUP %%%%
nodes = 268;
atlas = 'atlasname';
output_filename = 'DP_edges';
thresh = [0.001 0.005 0.01 0.02 0.05 0.1];

% Directory Setup
inputdir = fullfile(pwd(),'fc/');
outputdir = fullfile(pwd(), 'results/discriminatory_edges/');
functiondir = fullfile(pwd(), 'functions/');
addpath(genpath(inputdir),functiondir);
% ----------------------------------------


% upper triangle indices for reasembeling vectors to matrices 
FC_mat = zeros(nodes);
upper = find(triu(ones(length(FC_mat(1,:))),1));

% read out the subject ids from behavioural file for later loading of csv files with fc matrices
subid_path = fullfile(pwd(),'behavioral_data/Behaviour.csv');
id = readtable(subid_path);
sample_subs = table2array(id(:,1));

% Load fc matrices; final dimensions must be ROI by ROI
cd(inputdir);
rest_1_mats = zeros(length(FC_mat(1,:)),length(FC_mat(1,:)),length(sample_subs));

for i = 1:length(sample_subs)
    sub_i = sample_subs(i);
    fc_vec = csvread(strjoin({atlas, num2str(sub_i), 'rest1.csv'}, '_')); %read in functional connectivity 
    fc = reasemble_FC_mat(upper, fc_vec, length(FC_mat(1,:)), 1); % this step can be skipped if fuctional connectivity is already matrix
    rest_1_mats(:,:,i) = fc;
end
clear fc fc_vec file_i %only keeping rest_1_mats

sample_subs = table2array(id(:,1));
rest_2_mats = zeros(length(FC_mat(1,:)),length(FC_mat(1,:)),length(sample_subs));

for i = 1:length(sample_subs)
    sub_i = sample_subs(i);
    fc_vec = csvread(strjoin({atlas, num2str(sub_i), 'rest2.csv'}, '_'));
    fc = reasemble_FC_mat(upper, fc_vec, length(FC_mat(1,:)), 1);
    rest_2_mats(:,:,i) = fc;
end
clear fc fc_vec file_i %only keeping rest_2_mats


% From here on down, this is the original script by Finn et al., 2015 available at:
% https://www.nitrc.org/frs/?group_id=51

all_default_se1 = rest_1_mats;
all_default_se2 = rest_2_mats;

dim = size(all_default_se1);
no_sub = dim(3);
no_node = dim(1);

aa = ones(no_node, no_node);
aa_upp = triu(aa, 1);
upp_id = find(aa_upp);
upp_len = length(upp_id);


all_default_se1 = reshape(all_default_se1, no_node*no_node, dim(3));
all_default_se1 = all_default_se1(upp_id,:);
all_default_se2 = reshape(all_default_se2, no_node*no_node, dim(3));
all_default_se2 = all_default_se2(upp_id,:);

nn_default_s1 = zeros(upp_len, no_sub); % save the normalized connectivity matrices for condition1
nn_default_s2 = zeros(upp_len, no_sub); % save the normalized connectivity matrices for condition2

% normalize the connectivity matrices for each subject
for i=1:no_sub;
    c1 = all_default_se1(:,i);
    c1 = c1-mean(c1);
    nn_default_s1(:,i) = c1/norm(c1); 
    
    c2 = all_default_se2(:,i);
    c2 = c2-mean(c2);
    nn_default_s2(:,i) = c2/norm(c2);
end

edge_stat = zeros(1, upp_len); % save the differential statistics 


for edge = 1:upp_len;
    cur_edge_corr = nn_default_s1(edge,:)'*nn_default_s2(edge,:);  % this is of size no_sub by no_sub
    
    for i = 1:no_sub;
        cur_row = cur_edge_corr(i,:);
        cur_col = cur_edge_corr(:,i);
        cur_dig = cur_edge_corr(i,i);
        
        cur_p = (length(find(cur_row>cur_dig))+length(find(cur_col>cur_dig)))/(2*no_sub); 
        
        if (isnan(cur_row) == 1) 
          cur_p = nan;
        end
        
        if( cur_p==0)
            cur_p = 0.0000001; % because log(0) is not defined so we set a minimum here 
        
        end
        
        edge_stat(edge) = edge_stat(edge)+log(cur_p);     
    end
end

sort_stat = sort(-edge_stat(:), 'descend');

% convert the differential statistic vector to a matrix
edge_stat_matrix = zeros(no_node, no_node);
edge_stat_matrix(upp_id) = edge_stat;
edge_stat_matrix = edge_stat_matrix + edge_stat_matrix';

% threshold the differential statistic matrix and create binary mask for visualization

th_len = round(thresh*upp_len); 

for i=1:length(thresh);
    mask = (-edge_stat_matrix>=sort_stat(th_len(i)));
    dlmwrite([outputdir, output_filename, '_',  num2str(thresh(i)),'.txt'], mask, ' ');
    clear mask;
end