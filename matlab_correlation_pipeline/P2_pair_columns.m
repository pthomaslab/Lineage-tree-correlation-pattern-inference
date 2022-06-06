%% P2_pair_columns.m
% Pairs columns according to cell ID pairs
% Please read README.txt for full details

%% INPUT DATA
% Uses [celltype]_paired_cell_ids.mat from P1_pair_cells.m

prompt1 = 'Type cell name (without _paired_cell_ids extension): ';
answer1 = input(prompt1,'s');

load(strcat(answer1,'_paired_cell_ids.mat'));

prompt2 = 'What is it that you are pairing? ';
% type something to go into output file name, eg 'IDT'
%   for interdivision time
answer2 = input(prompt2,'s');

prompt3 = strcat('Which column (number) do you want to pair? ');
answer3 = input(prompt3);

%% OUTPUT
% [celltype]_[choice]_pairs.mat which will be used in
%   bootstrapping_single.m

%%
col = answer3;
lgt = length(col);

js = [1:lgt;lgt+1:2*lgt];

%% CREATE MATRICES WITH CELL CYCLE DURATION

MD_size = size(MD_pairs);
MD_data = zeros(MD_size(1),2);
for i = 1:MD_size(1);
    
    for j = 1:2,
    loc = MD_pairs(i,j) == data(:,1);
    if sum(loc)~=0
        MD_data(i,js(j,:)) = data(loc,col);
    end
    end

end
MD_data = MD_data(all(MD_data(:,:)~=0,2),:);
MD_pairs = MD_pairs(all(MD_data(:,:)~=0,2),:);

SS_size = size(SS_pairs);
SS_data = zeros(SS_size(1),2);
for i = 1:SS_size(1);
    
    for j = 1:2;
    loc = SS_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    SS_data(i,js(j,:)) = data(loc,col);
    end
    end

end
SS_data = SS_data(all(SS_data(:,:)~=0,2),:);

CC_size = size(CC_pairs);
CC_data = zeros(CC_size(1),2);
for i = 1:CC_size(1);
    
    for j = 1:2;
    loc = CC_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    CC_data(i,js(j,:)) = data(loc,col);
    end
    end

end
CC_data = CC_data(all(CC_data(:,:)~=0,2),:);

GG_size = size(GG_pairs);
GG_data = zeros(GG_size(1),2);
for i = 1:GG_size(1);
    
    for j = 1:2;
    loc = GG_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    GG_data(i,js(j,:)) = data(loc,col);
    end
    end

end
GG_data = GG_data(all(GG_data(:,:)~=0,2),:);

GR_size = size(GR_pairs);
GR_data = zeros(GR_size(1),2);
for i = 1:GR_size(1);
    
    for j = 1:2;
    loc = GR_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    GR_data(i,js(j,:)) = data(loc,col);
    end
    end

end
GR_data = GR_data(all(GR_data(:,:)~=0,2),:);

GR2_size = size(GR2_pairs);
GR2_data = zeros(GR2_size(1),2);
for i = 1:GR2_size(1);
    
    for j = 1:2;
    loc = GR2_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    GR2_data(i,js(j,:)) = data(loc,col);
    end
    end

end
GR2_data = GR2_data(all(GR2_data(:,:)~=0,2),:);

GR3_size = size(GR3_pairs);
GR3_data = zeros(GR3_size(1),2);
for i = 1:GR3_size(1);
    
    for j = 1:2;
    loc = GR3_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    GR3_data(i,js(j,:)) = data(loc,col);
    end
    end

end
GR3_data = GR3_data(all(GR3_data(:,:)~=0,2),:);

GR4_size = size(GR4_pairs);
GR4_data = zeros(GR4_size(1),2);
for i = 1:GR4_size(1);
    
    for j = 1:2;
    loc = GR4_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    GR4_data(i,js(j,:)) = data(loc,col);
    end
    end

end
GR4_data = GR4_data(all(GR4_data(:,:)~=0,2),:);

GR5_size = size(GR5_pairs);
GR5_data = zeros(GR5_size(1),2);
for i = 1:GR5_size(1);
    
    for j = 1:2;
    loc = GR5_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    GR5_data(i,js(j,:)) = data(loc,col);
    end
    end

end
GR5_data = GR5_data(all(GR5_data(:,:)~=0,2),:);

GR6_size = size(GR6_pairs);
GR6_data = zeros(GR6_size(1),2);
for i = 1:GR6_size(1);
    
    for j = 1:2;
    loc = GR6_pairs(i,j) == data(:,1);
    if sum(loc)~=0
    GR6_data(i,js(j,:)) = data(loc,col);
    end
    end

end
GR6_data = GR6_data(all(GR6_data(:,:)~=0,2),:);

save(strcat(answer1,'_',answer2,'_pairs'),'MD_pairs','MD_data','SS_data','CC_data','GG_data',...
    'GR_data','GR2_data','GR3_data','GR4_data','GR5_data','GR6_data','data','col');

clear