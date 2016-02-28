%%
% author: Chayan Sarkar
% email: c.sarkar@tudelft.nl
%%
function [n, total, S, link, loc] = deployment(deploy,datatype)
    
    if(deploy == 1)
        total = 5000;
        nodes = [1,2,3,4,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,...
           22,23,24,25,27,28,29,30,31,32,33,34,35,36,37,38,...
           39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54]';
        n = length(nodes);
        if(datatype == 1)
            path = 'IntelLab/temp';
        else
            path = 'IntelLab/humi';
        end
        S = zeros(total,n);
        for i=1:n
            name = strcat(path, int2str(nodes(i)));
            Sen = load(name);
            S(1:total,i) = Sen(1:total,1);
        end
        link = load('IntelLab/link_intellab');
        loc = get_location('IntelLab/locations',nodes,n);
    else
        total = 432;
        nodes = load('GreenOrb/nodes');
        if(datatype == 1)
            path = 'GreenOrb/temp';
        else
            path = 'GreenOrb/humi';
        end
        n = length(nodes);
        S = zeros(total,n);
        for i=1:n
            name = strcat(path, int2str(nodes(i)));
            Sen = load(name);
            S(1:total,i) = Sen(1:total,1);
        end
        link = load('GreenOrb/link_greenorb');
        loc = get_location('GreenOrb/locations',nodes,n);
    end
end

function loc = get_location(fname,nodes,n)
    ll = load(fname);
    loc = zeros(n,2);
    for i=1:n
        for j=1:length(ll)
            if(ll(j,1) == nodes(i))
                loc(i,1) = ll(j,2);
                loc(i,2) = ll(j,3);
                break;
            end
        end
    end
end
