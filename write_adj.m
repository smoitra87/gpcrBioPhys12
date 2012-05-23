%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write an adjacency matrix as an outpout text file
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = write_adj(fname,adj_mat)
    [idx,idy] = find(triu(adj_mat));
  
    [s,s_idx]= sort(idx);
    idx = idx(s_idx);
    idy = idy(s_idx);
    fout = fopen(fname,'w');
    for i = 1:length(idx)
        fprintf(fout,'%d %d\n',idx(i),idy(i));
     
    end
    fclose(fout);
    
    % Find overlap with old matrix
    
end

