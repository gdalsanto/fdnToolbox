clear; clc; close all;
addpath(genpath('../../fdnToolbox'))

init_path = "../output/test_fdn_rirs.mat";
optim_path = "../output/test_fdn_rirs_optim.mat";
%% Analyse reference RIR 

init_fdns = load(init_path);
optim_fdns = load(optim_path);

n_fdns = length(init_fdns.onset_time);

%% construct FDN

fs = init_fdns.fs(1);
D = zeros(1, 1);
irLen = 2*fs;
for i = 1:n_fdns
    [B, C, A, delays] = unpack_params(init_fdns, i);
    [residues, poles, direct, isConjugatePolePair, metaData] = dss2pr(delays, A, B, C, D);
    
end

%% functions 
function Y = skew(X)
    X = triu(X,1);
    Y = X - transpose(X);
end

function [B,C,A,delays] = unpack_params(struct, indx)
    B = double(struct.input_gains(indx, :))';
    C = double(struct.output_gains(indx, :));
    U = double(expm(skew(squeeze(struct.feedback_matrix(indx, :, :)))));
    g = double(struct.attenuation(indx));
    delays = double(struct.delays(indx, :));
    Gamma = diag(g.^delays);
    A = Gamma*U;    
end