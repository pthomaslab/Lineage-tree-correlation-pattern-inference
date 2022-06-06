%% bootstrapping_single.m
% finds correlation coefficients for choice and cell type
% Please read README.txt for full details

%% INPUT DATA
% Uses [celltype]_[choice]_pairs.mat from P2_pair_columns.m

prompt1 = 'Cell type ';
answer1 = input(prompt1,'s');

prompt2 = 'What did you choose to pair ';
answer2 = input(prompt2,'s');

load(strcat(answer1,'_',answer2,'_pairs','.mat'));

rng(10)

%% OUTPUT
% [celltype]_[choice]_corrs.txt which will be used in
%   bayesian_sampling.jl
% For format of this file see README.txt

%% SET-UP
% CREATE VECTOR WITH ALL VALID CELL CYCLE DURATIONS

data_all = data(:,col);
v = data(:,1) == data(:,2);
data_all(v,:) = [];

%% REMOVE ALL 0s
data_all(data_all==0) = [];

v = sum(MD_data==0,2)>0;
MD_data(v,:)=[]; 

v = sum(SS_data==0,2)>0;
SS_data(v,:)=[];

v = sum(CC_data==0,2)>0;
CC_data(v,:)=[]; 

v = sum(GG_data==0,2)>0;
GG_data(v,:)=[]; 

v = sum(GR_data==0,2)>0;
GR_data(v,:)=[]; 

v = sum(GR2_data==0,2)>0;
GR2_data(v,:)=[]; 

v = sum(GR3_data==0,2)>0;
GR3_data(v,:)=[]; 
v = sum(GR4_data==0,2)>0;
GR4_data(v,:)=[]; 
v = sum(GR5_data==0,2)>0;
GR5_data(v,:)=[]; 
v = sum(GR6_data==0,2)>0;
GR6_data(v,:)=[]; 

% repeats = 10000
repeats = 10000;
% CONFIDENCE INTERVAL LOCATIONS
% eg if 1000 repeats, 95% CI is at 25th and 975th value
L = round(repeats*0.025);
U = round(repeats*0.975);

for i = 1:repeats
    MD_sample = datasample(MD_data,length(MD_data));
    SS_sample = datasample(SS_data,length(SS_data));
    CC_sample = datasample(CC_data,length(CC_data));
    GG_sample = datasample(GG_data,length(GG_data));
    
    GR_sample = datasample(GR_data,length(GR_data));
    GR2_sample = datasample(GR2_data,length(GR2_data));
    GR3_sample = datasample(GR3_data,length(GR3_data));
    GR4_sample = datasample(GR4_data,length(GR4_data));
    GR5_sample = datasample(GR5_data,length(GR5_data));
    GR6_sample = datasample(GR6_data,length(GR6_data));
    
    data_sample = datasample(data_all,length(data_all));
    
    if length(MD_sample) > 2;
    MD_corr = corrcoef(MD_sample);
    else
    MD_corr = ones(2)*NaN;
    end
    if length(SS_sample) > 2;
    SS_corr = corrcoef(SS_sample);
    else
    SS_corr = ones(2)*NaN;
    end
    if length(CC_sample) > 2;
    CC_corr = corrcoef(CC_sample);
    else
    CC_corr = ones(2)*NaN;
    end
    if length(GG_sample) > 2;
    GG_corr = corrcoef(GG_sample);
    else
    GG_corr = ones(2)*NaN;
    end
    if length(GR_sample) > 2;
    GR_corr = corrcoef(GR_sample);
    else
    GR_corr = ones(2)*NaN;
    end
    if length(GR2_sample) > 2;
    GR2_corr = corrcoef(GR2_sample);
    else
    GR2_corr = ones(2)*NaN;
    end
    if length(GR3_sample) > 2;
    GR3_corr = corrcoef(GR3_sample);
    else
    GR3_corr = ones(2)*NaN;
    end
    if length(GR4_sample) > 2;
    GR4_corr = corrcoef(GR4_sample);
    else
    GR4_corr = ones(2)*NaN;
    end
    if length(GR5_sample) > 2;
    GR5_corr = corrcoef(GR5_sample);
    else
    GR5_corr = ones(2)*NaN;
    end
    if length(GR6_sample) > 2;
    GR6_corr = corrcoef(GR6_sample);
    else
    GR6_corr = ones(2)*NaN;
    end

    data_mean = mean(data_sample);
    data_cov = cov(data_sample);
    
    sample_vals(i,:) = [ data_mean data_cov ...
        MD_corr(1,2) SS_corr(1,2) CC_corr(1,2) GG_corr(1,2) ...
        GR_corr(1,2) GR2_corr(1,2) GR3_corr(1,2) GR4_corr(1,2) GR5_corr(1,2) GR6_corr(1,2)];
end

%% CALCULATE TRUE MEAN, COV, CORR

    if length(MD_data) > 2;
    mdc = corrcoef(MD_data);
    else
    mdc = ones(2)*NaN;
    end
    if length(SS_data) > 2;
    ssc = corrcoef(SS_data);
    else
    ssc = ones(2)*NaN;
    end
    if length(CC_data) > 2;
    ccc = corrcoef(CC_data);
    else
    ccc = ones(2)*NaN;
    end
    if length(GG_data) > 2;
    ggc = corrcoef(GG_data);
    else
    ggc = ones(2)*NaN;
    end
    if length(GR_data) > 2;
    grc = corrcoef(GR_data);
    else
    grc = ones(2)*NaN;
    end
    if length(GR2_data) > 2;
    gr2c = corrcoef(GR2_data);
    else
    gr2c = ones(2)*NaN;
    end
    if length(GR3_data) > 2;
    gr3c = corrcoef(GR3_data);
    else
    gr3c = ones(2)*NaN;
    end
    if length(GR4_data) > 2;
    gr4c = corrcoef(GR4_data);
    else
    gr4c = ones(2)*NaN;
    end
    if length(GR5_data) > 2;
    gr5c = corrcoef(GR5_data);
    else
    gr5c = ones(2)*NaN;
    end
    if length(GR6_data) > 2;
    gr6c = corrcoef(GR6_data);
    else
    gr6c = ones(2)*NaN;
    end

true_vals = [ mean(data_all) ...
    cov(data_all) ...
    mdc(1,2) ssc(1,2) ccc(1,2) ggc(1,2) ...
    grc(1,2) gr2c(1,2) gr3c(1,2) gr4c(1,2) gr5c(1,2) gr6c(1,2) ];

sz = size(sample_vals);

covbias = cov(data_all) - mean(sample_vals(:,2));
sample_vals(:,2) = sample_vals(:,2) + covbias;

for i = 1:sz(2);
    
    adj_vals(:,i) = sample_vals(:,i) - true_vals(i);
    
end

adj_vals = sort(adj_vals,1);
adj_vals(:,1) = adj_vals(:,1).*sqrt(length(data_all));

dataset = [true_vals ; ...
    adj_vals(L,:) ; ...
    adj_vals(U,:) ;
    [ std(data_all) std(sample_vals(:,2:end))] ] ;

dataset = dataset';

% CV calculation
cv = sqrt(sample_vals(:,2))./sample_vals(:,1);
cvmean = mean(cv);
cvstd = std(cv);
cvout = [cvmean ; cvstd];

% save(strcat(answer1,'_',answer2,'_bootstrapped'))
writematrix(dataset,strcat(answer1,'_',answer2,'_corrs.txt'),'Delimiter','tab')
% not used, coefficient of variation
% writematrix(cvout,strcat(answer1,'_',answer2,'_cv.txt'),'Delimiter','tab')
