%% P1_pair_cells.m
% Pairs cell IDs together
% Please read README.txt for full details

%% INPUT DATA
% 1 - CELL ID | 2 - PARENT ID
% The remaining columns can be anything you want that you want to pair

prompt1 = 'Type file name (without .txt extention) ';
answer1 = input(prompt1,'s');
data = table2array(readtable(strcat(answer1,'.txt')));

prompt2 = 'Type the celltype: ';
answer2 = input(prompt2,'s');

%% OUTPUT
% [celltype]_paired_cell_ids.mat which will be used in
%   P2_pair.columns.m

%% MOTHER DAUGHTER

MD_pairs = zeros(length(data),2);

for i = 1:length(data)
    
    cell_id = data(i,1);
    mother_id = data(i,2);
    
    mother_loc = data(:,1) == mother_id;
    gmother_id = data(mother_loc,2);
    
    if isempty(gmother_id) == 0 && ...
            gmother_id ~= mother_id && ...
            mother_id ~= cell_id 
        MD_pairs(i,:) = [cell_id, mother_id];
    end
end


%% SISTER SISTER

s = 1;
SS_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell1_id = data(i,1);
    mother1_id = data(i,2);
    
    for j = 1:length(data)
        cell2_id = data(j,1);
        mother2_id = data(j,2);
        
        if cell1_id ~= cell2_id && ...
                cell1_id ~= mother1_id && ...
                cell2_id ~= mother2_id && ...
                mother1_id == mother2_id 
                SS_pairs(s,:) = [ cell1_id, cell2_id ];
                s = s+1 ;
        end
    end
end

%% COUSIN COUSIN

c = 1;
CC_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell1_id = data(i,1);
    mother1_id = data(i,2);
    
    mother1_loc = data(:,1) == mother1_id;
    gmother1_id = data(mother1_loc,2);
    
    for j = 1:length(data)
        cell2_id = data(j,1);
        mother2_id = data(j,2);
        
        gmother_loc = data(:,1) == mother2_id;
        gmother2_id = data(gmother_loc,2);
        
        if  cell1_id ~= cell2_id && ...
                cell1_id ~= mother1_id && ...
                cell2_id ~= mother2_id && ...
                isempty(gmother1_id) == 0 && ...
                isempty(gmother2_id) == 0 && ...
                mother1_id ~= gmother1_id && ...
                mother2_id ~= gmother2_id && ...
                mother1_id ~= mother2_id && ...
                gmother1_id == gmother2_id
            CC_pairs(c,:) = [ cell1_id, cell2_id ];
                c = c+1;
        end
    end
end

%% GRANDMOTHER GRANDDAUGHTER

g = 1;
GG_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell_id = data(i,1);
    mother_id = data(i,2);
    
    mother_loc = data(:,1) == mother_id;
    gmother_id = data(mother_loc,2);
    
    if isempty(gmother_id) == 0
        gmother_loc = data(:,1) == gmother_id;
        grgmother_id = data(gmother_loc,2);
    else
        grgmother_id = [];
    end
    
    if cell_id ~= mother_id && ...
            isempty(gmother_id) == 0 && ...
            mother_id ~= gmother_id && ...
            isempty(grgmother_id) == 0 && ...
            gmother_id > grgmother_id
        % IF VALID GRANDMOTHER
        GG_pairs(g,:) = [ cell_id , gmother_id ];
        g = g+1;
    end
end


%% GREAT GRANDDAUGHTER ( GR )

r = 1;
GR_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell_id = data(i,1);
    mother_id = data(i,2);
    
    mother_loc = data(:,1)==mother_id;
    gmother_id = data(mother_loc,2);
        
    if isempty(gmother_id) == 0
        gmother_loc = data(:,1) == gmother_id;
        grgmother_id = data(gmother_loc,2);
    else
        grgmother_id = [];
    end
    
    if isempty(grgmother_id) == 0
        grgmother_loc = data(:,1)==grgmother_id;
        gr2gmother_id = data(grgmother_loc,2);
    else
        gr2gmother_id = [];
    end
    
    if cell_id ~= mother_id && ...
           isempty(gmother_id) == 0 && ...
           mother_id ~= gmother_id && ...
           isempty(grgmother_id) == 0 && ...
           gmother_id ~= grgmother_id && ...
           isempty(gr2gmother_id) == 0 && ...
           grgmother_id ~= gr2gmother_id
        % IF VALID
        GR_pairs(r,:) = [ cell_id , grgmother_id ] ;
        r = r+1;
    end
end

%% GREAT GREAT GRANDDAUGHTER ( GR2 )

r = 1;
GR2_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell_id = data(i,1);
    mother_id = data(i,2);
    
    mother_loc = data(:,1)==mother_id;
    gmother_id = data(mother_loc,2);
    
    if isempty(gmother_id) == 0
        gmother_loc = data(:,1) == gmother_id;
        grgmother_id = data(gmother_loc,2);
    else
        grgmother_id = [];
    end
    
    if isempty(grgmother_id) == 0
        grgmother_loc = data(:,1)==grgmother_id;
        gr2gmother_id = data(grgmother_loc,2);
    else
        gr2gmother_id = [];
    end
    
    if isempty(gr2gmother_id) == 0
        gr2gmother_loc = data(:,1)==gr2gmother_id;
        gr3gmother_id = data(gr2gmother_loc,2);
    else
        gr3gmother_id = [];
    end
    
    if cell_id ~= mother_id && ...
           isempty(gmother_id) == 0 && ...
           mother_id ~= gmother_id && ...
           isempty(grgmother_id) == 0  && ...
           gmother_id ~= grgmother_id && ...
           isempty(gr2gmother_id) == 0 && ...
           grgmother_id ~= gr2gmother_id && ...
           isempty(gr3gmother_id) == 0 && ...
           gr2gmother_id ~= gr3gmother_id
       % IF VALID
        GR2_pairs(r,:) = [ cell_id , gr2gmother_id ];
        r = r+1;
    end
end

%% GREAT GREAT GREAT GRANDDAUGHTER ( GR3 )

r = 1;
GR3_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell_id = data(i,1);
    mother_id = data(i,2);
    
    mother_loc = data(:,1)==mother_id;
    gmother_id = data(mother_loc,2);
    
    if isempty(gmother_id) == 0
        gmother_loc = data(:,1) == gmother_id;
        grgmother_id = data(gmother_loc,2);
    else
        grgmother_id = [];
    end
    
    if isempty(grgmother_id) == 0
        grgmother_loc = data(:,1)==grgmother_id;
        gr2gmother_id = data(grgmother_loc,2);
    else
        gr2gmother_id = [];
    end
    
    if isempty(gr2gmother_id) == 0
        gr2gmother_loc = data(:,1)==gr2gmother_id;
        gr3gmother_id = data(gr2gmother_loc,2);
    else
        gr3gmother_id = [];
    end
    
    if isempty(gr3gmother_id) == 0
        gr3gmother_loc = data(:,1)==gr3gmother_id;
        gr4gmother_id = data(gr3gmother_loc,2);
    else
        gr4gmother_id = [];
    end
 
    if cell_id ~= mother_id && ...
           isempty(gmother_id) == 0 && ...
           mother_id ~= gmother_id && ...
           isempty(grgmother_id) == 0  && ...
           gmother_id ~= grgmother_id && ...
           isempty(gr2gmother_id) == 0 && ...
           grgmother_id ~= gr2gmother_id && ...
           isempty(gr3gmother_id) == 0 && ...
           gr2gmother_id ~= gr3gmother_id && ...
           isempty(gr4gmother_id) == 0 && ...
           gr3gmother_id ~= gr4gmother_id
       % IF VALID
        GR3_pairs(r,:) = [ cell_id , gr3gmother_id ];
        r = r+1;
    end
end

%% GREAT GREAT GREAT GREAT GRANDDAUGHTER ( GR4 )

r = 1;
GR4_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell_id = data(i,1);
    mother_id = data(i,2);
    
    mother_loc = data(:,1)==mother_id;
    gmother_id = data(mother_loc,2);
    
    if isempty(gmother_id) == 0
        gmother_loc = data(:,1) == gmother_id;
        grgmother_id = data(gmother_loc,2);
    else
        grgmother_id = [];
    end
    
    if isempty(grgmother_id) == 0
        grgmother_loc = data(:,1)==grgmother_id;
        gr2gmother_id = data(grgmother_loc,2);
    else
        gr2gmother_id = [];
    end
    
    if isempty(gr2gmother_id) == 0
        gr2gmother_loc = data(:,1)==gr2gmother_id;
        gr3gmother_id = data(gr2gmother_loc,2);
    else
        gr3gmother_id = [];
    end
    
    if isempty(gr3gmother_id) == 0
        gr3gmother_loc = data(:,1)==gr3gmother_id;
        gr4gmother_id = data(gr3gmother_loc,2);
    else
        gr4gmother_id = [];
    end
    
    if isempty(gr4gmother_id) == 0
        gr4gmother_loc = data(:,1)==gr4gmother_id;
        gr5gmother_id = data(gr4gmother_loc,2);
    else
        gr5gmother_id = [];
    end
 
    if cell_id ~= mother_id && ...
           isempty(gmother_id) == 0 && ...
           mother_id ~= gmother_id && ...
           isempty(grgmother_id) == 0  && ...
           gmother_id ~= grgmother_id && ...
           isempty(gr2gmother_id) == 0 && ...
           grgmother_id ~= gr2gmother_id && ...
           isempty(gr3gmother_id) == 0 && ...
           gr2gmother_id ~= gr3gmother_id && ...
           isempty(gr4gmother_id) == 0 && ...
           gr3gmother_id ~= gr4gmother_id && ...
           isempty(gr5gmother_id) == 0 && ...
           gr4gmother_id ~= gr5gmother_id
       % IF VALID
        GR4_pairs(r,:) = [ cell_id , gr4gmother_id ];
        r = r+1;
    end
end

%% GREAT GREAT GREAT GREAT GREAT GRANDDAUGHTER ( GR5 )

r = 1;
GR5_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell_id = data(i,1);
    mother_id = data(i,2);
    
    mother_loc = data(:,1)==mother_id;
    gmother_id = data(mother_loc,2);
    
    if isempty(gmother_id) == 0
        gmother_loc = data(:,1) == gmother_id;
        grgmother_id = data(gmother_loc,2);
    else
        grgmother_id = [];
    end
    
    if isempty(grgmother_id) == 0
        grgmother_loc = data(:,1)==grgmother_id;
        gr2gmother_id = data(grgmother_loc,2);
    else
        gr2gmother_id = [];
    end
    
    if isempty(gr2gmother_id) == 0
        gr2gmother_loc = data(:,1)==gr2gmother_id;
        gr3gmother_id = data(gr2gmother_loc,2);
    else
        gr3gmother_id = [];
    end
    
    if isempty(gr3gmother_id) == 0
        gr3gmother_loc = data(:,1)==gr3gmother_id;
        gr4gmother_id = data(gr3gmother_loc,2);
    else
        gr4gmother_id = [];
    end
    
    if isempty(gr4gmother_id) == 0
        gr4gmother_loc = data(:,1)==gr4gmother_id;
        gr5gmother_id = data(gr4gmother_loc,2);
    else
        gr5gmother_id = [];
    end
    
    if isempty(gr5gmother_id) == 0
        gr5gmother_loc = data(:,1)==gr5gmother_id;
        gr6gmother_id = data(gr5gmother_loc,2);
    else
        gr6gmother_id = [];
    end
 
    if cell_id ~= mother_id && ...
           isempty(gmother_id) == 0 && ...
           mother_id ~= gmother_id && ...
           isempty(grgmother_id) == 0  && ...
           gmother_id ~= grgmother_id && ...
           isempty(gr2gmother_id) == 0 && ...
           grgmother_id ~= gr2gmother_id && ...
           isempty(gr3gmother_id) == 0 && ...
           gr2gmother_id ~= gr3gmother_id && ...
           isempty(gr4gmother_id) == 0 && ...
           gr3gmother_id ~= gr4gmother_id && ...
           isempty(gr5gmother_id) == 0 && ...
           gr4gmother_id ~= gr5gmother_id && ...
           isempty(gr6gmother_id) == 0 && ...
           gr5gmother_id ~= gr6gmother_id
       % IF VALID
        GR5_pairs(r,:) = [ cell_id , gr5gmother_id ];
        r = r+1;
    end
end

%% GREAT GREAT GREAT GREAT GREAT GREAT GRANDDAUGHTER ( GR6 )

r = 1;
GR6_pairs = zeros(length(data),2);

for i = 1:length(data)
    cell_id = data(i,1);
    mother_id = data(i,2);
    
    mother_loc = data(:,1)==mother_id;
    gmother_id = data(mother_loc,2);
    
    if isempty(gmother_id) == 0
        gmother_loc = data(:,1) == gmother_id;
        grgmother_id = data(gmother_loc,2);
    else
        grgmother_id = [];
    end
    
    if isempty(grgmother_id) == 0
        grgmother_loc = data(:,1)==grgmother_id;
        gr2gmother_id = data(grgmother_loc,2);
    else
        gr2gmother_id = [];
    end
    
    if isempty(gr2gmother_id) == 0
        gr2gmother_loc = data(:,1)==gr2gmother_id;
        gr3gmother_id = data(gr2gmother_loc,2);
    else
        gr3gmother_id = [];
    end
    
    if isempty(gr3gmother_id) == 0
        gr3gmother_loc = data(:,1)==gr3gmother_id;
        gr4gmother_id = data(gr3gmother_loc,2);
    else
        gr4gmother_id = [];
    end
    
    if isempty(gr4gmother_id) == 0
        gr4gmother_loc = data(:,1)==gr4gmother_id;
        gr5gmother_id = data(gr4gmother_loc,2);
    else
        gr5gmother_id = [];
    end
    
    if isempty(gr5gmother_id) == 0
        gr5gmother_loc = data(:,1)==gr5gmother_id;
        gr6gmother_id = data(gr5gmother_loc,2);
    else
        gr6gmother_id = [];
    end
    
    if isempty(gr6gmother_id) == 0
        gr6gmother_loc = data(:,1)==gr6gmother_id;
        gr7gmother_id = data(gr6gmother_loc,2);
    else
        gr7gmother_id = [];
    end
 
    if cell_id ~= mother_id && ...
           isempty(gmother_id) == 0 && ...
           mother_id ~= gmother_id && ...
           isempty(grgmother_id) == 0  && ...
           gmother_id ~= grgmother_id && ...
           isempty(gr2gmother_id) == 0 && ...
           grgmother_id ~= gr2gmother_id && ...
           isempty(gr3gmother_id) == 0 && ...
           gr2gmother_id ~= gr3gmother_id && ...
           isempty(gr4gmother_id) == 0 && ...
           gr3gmother_id ~= gr4gmother_id && ...
           isempty(gr5gmother_id) == 0 && ...
           gr4gmother_id ~= gr5gmother_id && ...
           isempty(gr6gmother_id) == 0 && ...
           gr5gmother_id ~= gr6gmother_id && ...
           isempty(gr7gmother_id) == 0 && ...
           gr6gmother_id ~= gr7gmother_id
       % IF VALID
        GR6_pairs(r,:) = [ cell_id , gr6gmother_id ];
        r = r+1;
    end
end

%% REMOVE ZEROS
v = MD_pairs(:,1)==0 | MD_pairs(:,2)==0;
MD_pairs(v,:) = [];

v = SS_pairs(:,1)==0 | SS_pairs(:,2)==0;
SS_pairs(v,:) = [];

v = CC_pairs(:,1)==0 | CC_pairs(:,2)==0;
CC_pairs(v,:) = [];

v = GG_pairs(:,1)==0 | GG_pairs(:,2)==0;
GG_pairs(v,:) = [];

v = GR_pairs(:,1)==0 | GR_pairs(:,2)==0;
GR_pairs(v,:) = [];

v = GR2_pairs(:,1)==0 | GR2_pairs(:,2)==0;
GR2_pairs(v,:) = [];

v = GR3_pairs(:,1)==0 | GR3_pairs(:,2)==0;
GR3_pairs(v,:) = [];

v = GR4_pairs(:,1)==0 | GR4_pairs(:,2)==0;
GR4_pairs(v,:) = [];

v = GR5_pairs(:,1)==0 | GR5_pairs(:,2)==0;
GR5_pairs(v,:) = [];

v = GR6_pairs(:,1)==0 | GR6_pairs(:,2)==0;
GR6_pairs(v,:) = [];

save(strcat(answer2,'_paired_cell_ids'),'MD_pairs','SS_pairs','CC_pairs','GG_pairs',...
    'GR_pairs','GR2_pairs','GR3_pairs','GR4_pairs','GR5_pairs','GR6_pairs','data');
clear
