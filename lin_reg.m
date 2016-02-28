%%
% author: Chayan Sarkar
% email: c.sarkar@tudelft.nl
%%
function [beta, chi] = lin_reg(D, A, T)
trainingLen = round((3*T)/4);
testLen = T - trainingLen;

d_obs = zeros(trainingLen,1);
A_mat = zeros(trainingLen,2);
for i=1:trainingLen
    d_obs(i,1)   = D(i,1);
    A_mat(i,1:2) = [1 A(i,1)];
end
beta = (pinv(A_mat))*d_obs;


d_obs = zeros(testLen,1);
A_mat = zeros(testLen,2);
for i=1:testLen
    d_obs(i,1)   = D(i+trainingLen,1);
    A_mat(i,1:2) = [1 A(i+trainingLen,1)];
end
d_est = A_mat*beta;

err = d_obs - d_est;
if(var(D) == 0)
    chi = 1;
else
    chi = 1 - var(err)/var(D);
end
