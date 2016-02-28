%%
% author: Chayan Sarkar
% email: c.sarkar@tudelft.nl
%%
function parent = collection_tree(n, e_res, link)

    max_res = e_res(1);
    relay_cost = 2*max_res*ones(n,1);
    state = zeros(n,1);
    parent = -1*ones(n,1);


    parent(1) = 1;
    relay_cost(1) = 0;
    k = 1;
    cover = 1;

    while(cover < n)
        % update distance for the remaining nodes based on the newly select node
        for i=2:n
            if(state(i)==0 && link(i,k)>0 && i~=k)
                relay = 2*max_res - (e_res(i)+e_res(k));
                relay = relay + relay_cost(k);

                if(relay<relay_cost(i))
                  relay_cost(i) = relay;
                  parent(i) = k;
                end
            end
        end

        min_cost = 2*max_res;
        % find the node with minimum cost (among unmarked nodes)
        for i=2:n
            if(relay_cost(i)<min_cost && state(i)==0)
                min_cost = relay_cost(i);
                k = i;
            end
        end

        state(k) = 1;
        cover = cover+1;
    end
    
    
    pp = zeros(n-1,1);
    for i=2:n
        pp(i-1) = parent(i) - 1;
    end
    parent = pp;
end
    
%     max_res = e_res(1);
% relay_cost = 2*max_res*ones(n,1);
% state = zeros(n,1);
% parent = zeros(n,1);
% 
% 
% parent(1) = 1;
% relay_cost(1) = 0;
% k = 1;
% cover = 1;
% 
% while(cover < n)
%     % update distance for the remaining nodes based on the newly select node
%     for i=2:n
%         if(state(i)==0 && link(i,k)>0 && i~=k)
%             relay = 2*max_res - (e_res(i)+e_res(k));
%             relay = relay + relay_cost(k);
%              
%             if(relay<relay_cost(i))
%               relay_cost(i) = relay;
%               parent(i) = k;
%             end
%         end
%     end
%     
%     min_cost = 2*max_res;
%     % find the node with minimum cost (among unmarked nodes)
%     for i=2:n
%         if(relay_cost(i)<min_cost && state(i)==0)
%             min_cost = relay_cost(i);
%             k = i;
%         end
%     end
% 
%     state(k) = 1;
%     cover = cover+1;
%     
% end
% 
% pp = zeros(n-1,1);
% for i=1:n-1
%     pp(i) = parent(i+1) - 1;
% end
% parent = pp;