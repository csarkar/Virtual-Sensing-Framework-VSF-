%%
% author: Chayan Sarkar
% email: c.sarkar@tudelft.nl
%%
function [alpha, mu] = autoreg(D, p)

ss =size(D);
Tp = ss(1);
n = ss(2);

alpha = zeros(p, n);
mu = zeros(n,1);

for k=1:n
    d_vec = zeros(Tp-p,1);
    D_mat = zeros(Tp-p,p);

    for i=1:Tp-p
        d_vec(i,1)   = D(i+p,k);
        D_mat(i,1:p) = D(i+p-1:-1:i,k)';
    end

    alpha(:,k) = (pinv(D_mat))*d_vec;
    mu(k) = (Tp*0.005)/(D(1:Tp,k)'*D(1:Tp,k));
end

% function alpha = autoreg(D, Tp, p)
% 
% d_vec = zeros(Tp-p,1);
% D_mat = zeros(Tp-p,p);
% 
% for i=1:Tp-p
%     d_vec(i,1)   = D(i+p,1);
%     D_mat(i,1:p) = D(i+p-1:-1:i,1)';
% end
% 
% alpha = (pinv(D_mat))*d_vec;