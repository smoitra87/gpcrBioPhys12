%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program loads all the partial files as specified by hardcoded
% strings, constructs the adjacency matrix for the same and finds the
% percent overlap of edges with full alignment edges
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DO_PART = 2;
nRes =348;

if DO_PART == 1
    
    
%--------------------------------------------------------------------------
% Analyze each of the indivisual permutation graphs


% Permut graphs
perm_graph_list = { 'GPCR_ranga_perm_and.mat','chemo_perm_and.mat', ...
    'adrn_perm_and.mat' };
for i = 1:length(perm_graph_list)
    fname = perm_graph_list{i};
    load(fname);
    figure;
    spy(adjFinal);
    title(fname,'Interpreter','None');
    saveas(gcf,strcat(fname,'.png'),'png');
    
    % calculate graph density
    density = nnz(adjFinal) / (nRes*(nRes-1));
    fprintf('Density of %s : %f\n',fname,density);
end


end

if DO_PART == 2

%--------------------------------------------------------------------------
% Compare two or more graphs with each other

adrn_dir = 'adrn_out/partial';
chemo_dir = 'chemo_out/partial';
adrn_penalty = 10;
chemo_penalty = 20;
gpcr_rangaf = 'GPCR_ranga_perm_and.mat';

adjFinal = construct_adj_from_partial(adrn_dir,...
    'out_perm_adrn_',adrn_penalty);
save('adrn_partial_adj.mat','adjFinal');

figure;
spy(adjFinal);
title(strcat('Topology adrenergic penalty',num2str(adrn_penalty)));
saveas(gcf,strcat('adrn_top_',num2str(adrn_penalty),'.png'),'png');

% calculate graph density
density = nnz(adjFinal) / (nRes*(nRes-1));
fprintf('Density of %s : %f\n','adrn',density);

adjFinal = construct_adj_from_partial(chemo_dir,...
    'out_perm_chemo_',chemo_penalty);
save('chemo_partial_adj.mat','adjFinal');

figure;
spy(adjFinal);
title(strcat('Topology Chemokine penalty',num2str(chemo_penalty)));
saveas(gcf,strcat('chemo_top_',num2str(chemo_penalty),'.png'),'png');

% calculate graph density
density = nnz(adjFinal) / (nRes*(nRes-1));
fprintf('Density of %s : %f\n','chemo',density);

adjfiles = {'GPCR_ranga_perm_and.mat','chemo_partial_adj.mat',...
    'adrn_partial_adj.mat'};

adjfiles = perm_graph_list;

C = combnk(adjfiles,2);
for ii = 1:size(C,1) 
      compare_adj(C{ii,1},C{ii,2});
end

end