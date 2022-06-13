clear
%% INITIAL SET UP
global d

% CHOOSE NUMBER OF COMPONENTS
d = 2;

%% VARIABLES

prompt1 = 'Simulation name: ';
answer1 = input(prompt1,'s');

load(strcat(answer1,'_sim_initial_conditions.mat'))

x11 = original(1);
x12 = original(2);
x21 = original(3);
x22 = original(4);

lambda1 = original(5);
lambda2 = original(6);

gamma12 = original(7);

delta11 = original(8);
delta12 = original(9);
delta22 = original(10);

numtrees = ceil(cellcap/63);

%% VARIABLES
global mu theta covmatz beta
% mu : mean of latent variables y
mu = zeros(d,1);

% theta : inheritance matrix
theta = [x11 x12 ; x21 x22];

% beta : mean of noise terms e
beta = (eye(d) - theta)*(mu);

varmat = [lambda1 ; lambda2 ]*[lambda1 lambda2];

gammat = [ 1 gamma12 ; gamma12 1];

delmat = [ delta11 delta12 ; delta12 delta22 ];

S1 = varmat.*gammat;
S2 = varmat.*delmat;

covmatz = [ S1 S2 ; S2 S1 ];

initialtracks = zeros(numtrees,3);
ngen = 5;
prob_vec = [0 0 1];
all_tracks = zeros(numtrees*63,5); % CHANGE THIS

%% GENERATES PARENT ID

for j = 1:numtrees

clear parents x_tracks initialx

parents = branch(ngen,prob_vec);
parents = parents';

%% GENERATE INITIAL CELLULAR QUANTITIES
initialtracks(j,:) = initial_conditions(randi([1 length(initial_conditions)],1),3:5);
initialx = initialtracks(j,1:2);

%% DATA SET UP

% IN x_tracks - 1: CELL ID, 2: PARENT ID, 3 ONWARDS: LATENT VARIABLES
x_tracks(:,1) = 1:length(parents);
x_tracks(:,2) = parents;
x_tracks(1,3:3+d-1) = initialx;

% SIMULATE TREE FROM ONE INITIAL CONDITION
for i = 1:max(parents);
    ic = x_tracks(i,3:3+d-1);
    ic = ic';
    ckO = ck(ic);
    
    d1 = find(x_tracks(:,2)==i,1,'first');
    d2 = find(x_tracks(:,2)==i,1,'last');
    
    % Adds latent variables to end of x_tracks
    x_tracks(d1,3:3+d-1) = ckO(1,:);
    x_tracks(d2,3:3+d-1) = ckO(2,:);
end

x_tracks(1,2) = 1;
x_tracks(:,3+d) = sum(x_tracks(:,3:3+d-1),2);

point1 = find(all_tracks==0,1,'first')
if point1==1;
    all_tracks(point1:point1+62,:) = x_tracks
else;
point = point1-1;
add = all_tracks(point,1);
x_tracks(:,1:2) = x_tracks(:,1:2) + add;

all_tracks(point1:point1+62,:) = x_tracks
end

end

data = [all_tracks(:,1:2), all_tracks(:,5), all_tracks(:,3:4)];
data(:,3) = data(:,3)+answer3;
data(:,4:5) = data(:,4:5) + answer3/2;
data = datasample(data,round(0.85*length(data)),'Replace',false)
save(strcat(answer1,'_sim_final_data'))

writematrix(data,strcat(answer1,'_sim_final_data.txt'),'Delimiter','tab')

%% FUNCTION
function ckO = ck(y)

global theta covmatz beta d
 
% e : correlated noise terms
e = mvnrnd([beta;beta],covmatz);

e1 = e(1:d)';
e2 = e(d+1:end)';

d1 = theta*y + e1;
d2 = theta*y + e2;

ckO = [d1' ; d2'];

end