%%
% author: Chayan Sarkar
% email: c.sarkar@tudelft.nl
%%
function [beta,chi,companion,flag] = ...
    active_selection(S,Tp,n,th,e_res,parent,cnum)

    %-- find a collection tree among the nodes incorporating 
    %-- residual energy of the nodes in the link weight.
    %parent = collection_tree(n+1, 5*ones(n+1,1), link);

    %-- find the constellations of nodes, 
    %-- i.e., potential companions of a node.
    corr_comp = get_corr_companion(S, n, th, cnum);
    
    %-- finally choosing active set of nodes 
    %-- and assign companion for the dormant nodes.
    companion = active_node_set(n, e_res, corr_comp, parent);

    
    flag = 0;
    beta = zeros(2,n);
    chi = ones(n,1);
    for i=1:n
        if(companion(i) ~= 0)
            [beta(:,i), chi(i)] = lin_reg(S(1:Tp,i),...
                S(1:Tp,companion(i)),Tp);
            if(chi(i) <= 0)
                flag = flag + 1;
            end
        end
    end
end


function companion = active_node_set(n, e_res, corr_comp, parent)
    state = zeros(n,1);
    state(1) = 1;
    companion = zeros(n,1);
    companion(1) = 0;
    cover = 1;

    while(cover<n)
        e_max = 0;
        i_max = 0;
        for i=2:n
            if(e_max < e_res(i) && state(i)==0)
                e_max = e_res(i);
                i_max = i;
            end
        end
        [cover,state,companion] = mark_as_active(i_max,cover,...
            state,parent,companion,corr_comp,n);
    end
end


%%
function [cover,state,companion] = mark_as_active(i,cover,...
    state,parent,companion,corr_comp,n)

    if(state(i)==0)
        cover = cover+1;
    end

    state(i) = 1;
    companion(i) = 0;

    if(parent(i) ~= 0 && state(parent(i)) ~= 1)
        [cover,state,companion] = mark_as_active(parent(i),cover,...
            state,parent,companion,corr_comp,n);
    end

    for k=1:n
        if(state(k)==0 && corr_comp(i,k)==1)
            state(k) = -1;
            cover = cover+1;
            companion(k) = i;
        end
    end
end
