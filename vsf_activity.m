%%
% author: Chayan Sarkar
% email: c.sarkar@tudelft.nl
% file description: main file describe VSF method
%%
function [err,tx,ttx,rmse,rmse2,eres] = vsf_activity(cnum)
    
    global S P;
    global n total Tp Op Rp p;
    global th;
    global e_res;

    global alpha mu;
    global beta;
    global parent;
    global companion;
    global T e_active e_dormant e_send e_receive
    
    global sCache pCache;
    global TRAINING OPERATION REVALID;
    global ZEROVALUE;
    

    
    %% setting system parameters
    %-- constant values
    TRAINING = 1;
    OPERATION = 2; 
    REVALID = 3;
    ZEROVALUE = -9999.0;

    %-- energy measurement parameters (approximate estimate)
    %-- the energy values are calculated for sky-node from its datasheet;
    cpu_active = 1.8;
    cpu_sleep = (0.001*5.1);
    radio_idle = (0.001*21);
    radio_sleep = (0.001*1);
    
    %-- every sensor node is assumed to be sense and transmit every 31sec;
    %-- the duty cycling is assumed to be 1 percent of the sensing period.
    %-- the energy values are in uJ.
    e_sense = 0.33;
    e_active = (0.95*T)*(cpu_sleep+radio_sleep)*3 +...
        (0.05*T)*(cpu_active+radio_idle)*3 + e_sense;
    e_dormant = T*(cpu_sleep+radio_sleep)*3;
    e_send = 4.74;
    e_receive = 5.23;
    eres = e_res;
    
    
    %% initial parameter setup
    tx = ones(n,1);
    ttx = 2*n;
    
    % two small cache is maintained at the sensor nodes;
    % the current sensed data is predicted based on a certain 
    % number of past values, which are stored in the cache.
    sCache = ZEROVALUE*ones(p,n);
    pCache = zeros(p,n);
    
    for j=1:n
        sCache(1,j) = S(1,j);
        P(1,j) = S(1,j);    
        pCache(1:p,j) = P(1,j);
    end
    
    alpha = (1/p)*ones(p,n);
    mu = 0.005*ones(n,1);
    
    companion(1:n) = 0;
    period = TRAINING;
    errTh = th/2;
    
    mode = zeros(total,1);
    %% VSF starts
    i = 2;
    prev_i = 1;
    while(i<=total)
        if(period==TRAINING || period==REVALID)
            mode(i) = period;
        elseif(companion(12)==0)
            mode(i) = period;
        else
            mode(i) = 0;
        end
        
        %-- first predict for all the active nodes %
        for j=1:n
            if(companion(j) == 0)
                %-- update sCache
                sCache(2:p,j) = sCache(1:p-1,j);
                sCache(1,j) = S(i,j);

                % predict the data
                P(i,j) = pCache(1:p,j)'*alpha(1:p,j);
                e = S(i,j) - P(i,j);
                %-- if prediction error crosses error threshold, 
                %-- use the sensed value
                if(abs(e) > errTh)
                    P(i,j) = S(i,j);
                    for k=1:p
                        if(sCache(k,j) ~= ZEROVALUE)
                            pCache(k,j) = sCache(k,j);
                        end
                    end
                    alpha(1:p,j) = alpha(1:p,j) + pCache(1:p,j)*e*mu(j);
                    tx(j) = tx(j)+1;
                    ttx = ttx+1;
                    
                    e_res(j) = e_res(j) - (e_active + e_send);
                    jj = parent(j);
                    while(jj ~= 0)
                        ttx = ttx+1;
                        e_res(jj) = e_res(jj) - (e_receive + e_send);
                        jj = parent(jj);
                    end
                else
                    pCache(2:p,j) = pCache(1:p-1,j);
                    pCache(1,j) = P(i,j);
                    e_res(j) = e_res(j) - e_active;
                end
            end
        end

        %-- then predict for dormant nodes %
        for j=1:n
            if(companion(j) ~= 0)
                c = companion(j);
                cData = P(i,c);
                sCache(2:p,j) = sCache(1:p-1,j);
                sCache(1,j) = ZEROVALUE;

                P(i,j) = [1; cData]'*beta(:,j);
                e_res(j) = e_res(j) - e_dormant;
            end
        end
        i = i+1;
        
        if(period == TRAINING)
            if((i-prev_i) >= Tp)
                [alpha(:,1:n),mu(1:n,1)] = autoreg(P(i-Tp:i-1,1:n), p);
                [period,errTh] = decideNextMode(i,cnum);
                prev_i = i;
            end
        elseif(period == OPERATION)
            if((i-prev_i) >= Op)
               companion(1:n) = 0;
               period = REVALID;
               errTh = th/2;
               prev_i = i;
            end
        elseif(period == REVALID)
            if((i-prev_i) >= Rp)
                [period,errTh] = decideNextMode(i,cnum);
                prev_i = i;                
            end
        end
    end
    
    
    E = S - P;
    err = zeros(n*total,1);
    rmse = zeros(n,1);
    rmse2 = zeros(n,1);
    for j=1:n
        e=0;
        for i=1:total
            e = e+E(i,j)*E(i,j);
        end
        rmse(j) = sqrt(e/total);
        rmse2(j) = sqrt(e/(total-tx(j)));
        eres(j) = eres(j) - e_res(j);
        err((j-1)*total+1:j*total,1) = abs(E(1:total,j));
    end
    
end
    
    
%% this function decides the next operating mode of the sensor nodes    
function [period,errTh] = decideNextMode(i,cnum)
    global P;
    global beta delta companion;
    global link parent;
    global n Tp;
    global th corr_th e_res;
    global TRAINING OPERATION REVALID;
    
    %parent = collection_tree(n+1, [5*1e5; e_res], link);
    [beta,delta,companion,flag] = active_selection...
                    (P(i-Tp:i-1,:),Tp,n,corr_th,e_res,parent,cnum);
                
    %% in the dynamic mode selection, an operational period is 
    % not guranteed to resume after an Revalid or a training period;
    % rather another revalidation or a new training period might start.
                
    if(flag == 0)
        period = OPERATION;
        errTh = th;
    elseif(flag == 1)
        companion(1:n) = 0;
        period = REVALID;
        errTh = th/2;
    elseif(flag >= 2)
        companion(1:n) = 0;
        period = TRAINING;
        errTh = th/2;
    end
end
    
    
%% active nodes
function [pData,tx,ttx] = AVS(sData,j,th,tx,ttx)
    global alpha p mu;
    global parent;
    global e_res e_send e_receive e_active;
    global sCache pCache;
    global ZEROVALUE;

    % update sCache
    sCache(2:p,j) = sCache(1:p-1,j);
    sCache(1,j) = sData;
    
    % predict the data
    pData = pCache(1:p,j)'*alpha(1:p,j);
    e = sData - pData;
    % if prediction error crosses error threshold, 
    % use the sensed value
    if(abs(e) > th)
        pData = sData;
        for i=1:p
            if(sCache(i,j) ~= ZEROVALUE)
                pCache(i,j) = sCache(i,j);
            end
        end
        alpha(1:p,j) = alpha(1:p,j) + pCache(1:p,j)*e*mu(j);
        tx(j) = tx(j)+1;
        ttx = ttx+1;
        e_res(j) = e_res(j) - (e_active + e_send);
        jj = parent(j);
        while(jj ~= 0)
            ttx = ttx+1;
            e_res(jj) = e_res(jj) - (e_receive + e_send);
            jj = parent(jj);
        end
    else
        pCache(2:p,j) = pCache(1:p-1,j);
        pCache(1,j) = pData;
        e_res(j) = e_res(j) - e_active;
    end
end
    
    

      