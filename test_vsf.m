%%
% author: Chayan Sarkar
% email: c.sarkar@tudelft.nl
% file description: main file to test VSF
%%
%-- global variables --%
global S P;             % variable contains sensed and predicted data
global n total;         % n=number of nodes, total=sample size
global Tp Op Rp p;      
global link parent;
global th corr_th;
global e_res;           % residual energy of the nodes
global T;               % sensing period in seconds

%%-- parameter definition --%
Tp = 40;                % length of training period
Op = 20;                % length of operational period
Rp = 10;                % length of revalidation period
p = 3;                  % order of autoregression funtion
th = 1.5;               % estimation error threshold to decide when to transmit
corr_th = 0.95;         % threshold to decide high correlation

%%-- dataset related variable --%
INTELLAB = 1;
GREENORB = 2;
TEMP=1; HUMI=2;

%-- read deployment-specific sensor data from the files
[n,total,S,link] = deployment(INTELLAB,TEMP);
T = 31;                 % sensing period in INTELLAB was 31 seconds 
% [n,total,S,link] = deployment(GREENORB,TEMP);
% T = 10*60;              % sensing period in GreenOrb was 10 minutes

%-- create collection tree based on the adjacency matrix (link)  --%
parent = collection_tree(n+1, 5*ones(n+1,1), link);

%-- initialize the matrix to store predicted values --%
P = zeros(total,n);
    
%-- residual energy of the nodes --%    
e_res = 5*1e6*ones*ones(n,1);

%-- if a fixed number of correlated groups to be formed without 
% considering the correlation threshold set cnum;
% else set cnum to 1
cnum = 1;

%-- call the main function of VSF --%
% Return values
% err -> error matrix, which is simply S-P
% tx -> number of packet transmitted from the source
% ttx -> total number of packets within the network (tx + packets from the
% forwarding nodes)
% rmse -> RMSE per node, mean is calucalted using n
% rmse2 -> RMSE only for the values when data is predicted
% eres -> residual energy of the nodes
[err,tx,ttx,rmse,rmse2,eres] = vsf_activity(cnum);

%% statistics and plot
disp([sum(tx),ttx]);
disp(mean(rmse));

errbin = 0:0.05:4;
N = hist(err,errbin);

ss = sum(N);
cdf = 0;
NN = zeros(length(errbin));
for i=1:length(errbin)
    cdf = cdf + N(i);
    NN(i) = cdf/ss;
end    
plot(errbin, NN, '*-');
