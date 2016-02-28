%%
% author: Chayan Sarkar
% email: c.sarkar@tudelft.nl
%%
%-- S: matrix conatining sensed data from the nodes.
%-- n: total number of nodes in the network;
%-- 1 sink, n-1 sening nodes.
%-- th: correlation threshold; used for creating constellation;
%-- in a constellation, correlation among the nodes is atleast 'th'.
%-- cnum: number of constellation that need to be created;

%-- constellation formation stops if either of the conditions 
%-- ('th' or 'cnum') are reached; if constellations are formed 
%-- based on only one condition, the other parameter is set to 
%-- minimum (th=0, cnum=1).

function corr_comp = get_corr_companion(S, n, th, cnum)

    %-- matrix showing potential companion of a node;
    %-- if (i,j)th element is 1, then jth node can be a
    %-- companion of the ith element; the matrix is symmetric.
    %-- node 1 is the sink node and it can be 
    %-- correlated only with itself.
    corr_comp = zeros(n,n);
    corr_comp(1,1) = 1;
    
    %-- correlation among the sensor nodes is found.
    corrmat = abs(corrcoef(S(:,2:n)));

    %-- find the constellations of sensor nodes (2nd to nth node);
    %-- any node within a constellation can be companion of
    %-- any other node of the same constellation.
    %-- constellations are found using hierarchical clustering,
    %-- where correlation among the nodes is used as 
    %-- the distance measure.
    corr_comp(2:n,2:n) = hierarchical_constellation(n-1,...
        corrmat, cnum, th);

end



% n: how many variables are there
% corrmat: correlation matrix
% cnum: number of constellation required
function cmat = hierarchical_constellation(n, corrmat, cnum, th)

    count = n;
    rmat = corrmat;
    v(1:n,1) = 1:n;
    l = ones(count,1);

    
    if(cnum<=0)
        cnum = 1;
    end
    
    %%  the clusters are being merged until the required 
    %-- number of clusters are left
    while(count > cnum)
        
        %-- select the two constellations that can be merged;
        %-- the selection is done based on the distance between 
        %-- two constellation; the constellation pair with minimal 
        %-- distance (maximum correlation among their nodes) among
        %-- the remaining constellations are selected.
        rr = 0;
        for i=1:count-1
            for j=i+1:count
                r = rmat(i,j);
                if(r > rr && i~=j)
                    rr = r;
                    ii = i;
                    jj = j;
                end
            end
        end
        
        %-- if the distance (correlation) between the selected 
        %-- constellation pairs is more than the defined (correlation 
        %-- is less than the threshold, no more merging is possible
        %-- and the clustering process is exited.
        if(rr<th)
            break;
        end

        
        %-- cut the rows and columns from the old correlation matrix, 
        %-- except the two rows that are being merged
        rrmat = zeros(count-1, count-1);
        vv = zeros(count-1,n);
        ll = zeros(count-1,1);
        k = 1;
        for i=1:count
            if(i ~= ii && i~= jj)
                rrmat(k,1:count-2) = [rmat(i,1:ii-1) ...
                    rmat(i,ii+1:jj-1) rmat(i,jj+1:count)];
                ll(k) = l(i);
                vv(k,1:ll(k)) = v(i,1:l(i));
                k = k+1;
            end
        end

        %-- merge the two clusters
        ll(count-1) = l(ii) + l(jj);
        vv(count-1,1:ll(count-1)) = [v(ii,1:l(ii)) v(jj,1:l(jj))];
        %-- also calculate the correlation b/t the new (merged) 
        %-- and the old clusters
        for i=1:count-2
            rrmat(i,k) = getR(vv(i,:)',ll(i),vv(count-1,:)',...
                ll(count-1),corrmat);
            rrmat(k,i) = rrmat(i,k);
        end
        rrmat(count-1,count-1) = rr;
        
        %-- replacing the old clusters with new clusters.
        rmat = rrmat;
        v = vv;
        l = ll;
        count = count-1;
    end
    
    %%
    cmat = zeros(n,n);
    for k=1:count 
        for i=1:l(k)
            ii = v(k,i);
            for j=1:l(k)
                jj = v(k,j);
                cmat(ii,jj) = 1;
                cmat(jj,ii) = 1;
            end
        end
    end
end


%%  find the distance (relative correlation) between 
%-- two constellations; v1, v2 are the vectors with the 
%-- constelllation members and l1, l2 are their respective lengths.
%-- corrmat contains the original correlation among the nodes
%-- based on their sensed data.
function min = getR(v1, l1, v2, l2, corrmat)
    
    min = 1;
    for i=1:l1
        for j=1:l2
            if(corrmat(v1(i), v2(j)) < min)
                min = corrmat(v1(i), v2(j));
            end
        end
    end
end
