
% This function reasembles an fc matrix from a vector extracted 
% from the upper triangle of an fc matrix. 
%
% Input:
%   'FC_vec' is the vector extracted from the upper triangle. 
%   'upper' are the indeces of the original positions in the matrix
%           used to extract the upper triangle.
%   'n_nodes' is the number of nodes
%   'diagonal' is what integer should appear on the diagonal 
%
% Output:
%   'FC_M' is the full FC matrix
% 
% Example:
% FC_mat = csvread('C:\Users\Martin Gell\Desktop\Projects\FC\HCP\Measures\Static_mats\StaticMat_101006_rest1.csv');
% FC_vec = csvread('C:\Users\Martin Gell\Desktop\Projects\FC\HCP\Measures\Static\Static_101006_rest1.csv');
% upper = find(triu(ones(length(FC_mat(1,:))),1));
% FC_M = reasemble_FC_mat(upper, FC_vec, 268, Inf);
% all(all(FC_M == FC_mat))


function FC_M = reasemble_FC_mat(upper, FC_vec, n_nodes, diagonal)
FC_M_upper = zeros(n_nodes,n_nodes);
FC_M_upper(upper) = FC_vec;
FC_M_lower = FC_M_upper';
FC_M = FC_M_lower + FC_M_upper + diag(repelem(diagonal,n_nodes));
end

