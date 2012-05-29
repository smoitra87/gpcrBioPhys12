%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compares adjacency matrices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = compare_adj(adjf1,adjf2)

load(adjf1)
adj_orig = full(adjFinal);

load(adjf2)
adj_new = full(adjFinal);

% common edges
[bla,adjf1,bla] = fileparts(adjf1);
[bla,adjf2,bla] = fileparts(adjf2);

offset = length('_perm_and');

commonf = strcat(['common_' adjf1(1:end-offset) '_'...
    adjf2(1:end-offset) '.dat']);
fprintf('Common file is %s\n',commonf);

common = adj_orig.*adj_new;
fprintf('Percentage common edges: %f\n',length(find(common))/...
    length(find(adj_new)));

write_adj(commonf,common);

% different edges
difff1 = strcat(['diff_' adjf1(1:end-offset) '_'...
    adjf2(1:end-offset) '.dat']);
fprintf('Diff1 file is %s\n',difff1);

diff1 = adj_orig - adj_new;
diff1(diff1<0) =0;
write_adj(difff1,diff1);

difff2 = strcat(['diff_' adjf2(1:end-offset) '_'...
    adjf1(1:end-offset) '.dat']);
fprintf('Diff2 file is %s\n',difff2);

diff2 = adj_new - adj_orig;
diff2(diff2<0)=0;
write_adj(difff2,diff2);


% % reachability analysis
% uc = triu(common);
% unew = triu(adj_new);
% 
% for ii = 1:size(unew,1)
%    disp(nnz(uc(ii,:))/nnz(unew(ii,:)))
% end
% 
% end



% hop-analysis


