%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Makes a call to compare_adjm.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

adjfiles = {'GPCR_ranga_perm_and.mat','chemo_perm_and.mat',...
    'adrn_perm_and.mat'};

C = combnk(adjfiles,2);
for ii = 1:size(C,1) 
      compare_adj(C{ii,1},C{ii,2})
end
