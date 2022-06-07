clear
%% INITIAL SET UP
global d

% CHOOSE NUMBER OF COMPONENTS
d = 2;

%% VARIABLES

prompt0 = 'mlp .txt file name (excluding .txt extension) ';
answer0 = input(prompt0,'s');

prompt1 = 'Simulation name: ';
answer1 = input(prompt1,'s');

prompt2 = 'Number of cells to sim ';
cellcap = input(prompt2);

original = importdata(strcat(answer0,'.txt'));

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

%% FIRST PASS - GENERATE INITIAL CONDITIONS

%% GENERATES PARENT ID
ngen = 10;
prob_vec = [0 0 1];

parents = branch(ngen,prob_vec);
parents = parents';

%% GENERATE INITIAL CELLULAR QUANTITIES
initialx = normrnd(0,1,[1,d]);

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

initial_conditions = x_tracks(1000:end,:)

save(strcat(answer1,'_sim_initial_conditions'),'initial_conditions','cellcap')

clear

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