%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consttructs adjacency matrix from partial using an AND rule
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [adjFinal]  = construct_adj_from_partial(file_dir,offset,penalty)

    fprintf('nargin = %d\n',nargin);
    if nargin < 3
        file_dir = 'adrn_out/partial';
        offset = 'out_perm_adrn_';
        penalty = 10;
    end
    nRes = 348;
    
    adjFinal = ones(nRes,nRes);
    
    for i=1:nRes
        fn = strcat(offset,num2str(penalty),'_',num2str(i),'.mat');
        fpath = fullfile(file_dir,fn);
        load(fpath);
        adjFinal(i,:) = full(adj)';
    end
    adjFinal = adjFinal .* adjFinal';
  
end