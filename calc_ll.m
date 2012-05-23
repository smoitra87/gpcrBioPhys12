%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates likelihoods of seqeunces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mutf = 'mut_oznur.fasta';
mut_seqs = fastaread(mutf);
modelf = 'GPCR_ranga_perm_and_fix_1_1.mat';
do_display = 1;
mutfoldf = 'mut_oznur.fold';

mut_msa = cell(length(mut_seqs),1);
testNdx= 1:length(mut_seqs);

for seqid = 1:length(mut_seqs) 
    mseq = mut_seqs(seqid).Sequence;
    mut_msa{seqid} = mseq;

end
nmsa = converttonumericmsa(mut_msa);

% Load the computed model
load(modelf,'nodeWeights','edgeWeights',...
    'edgeStruct')
foldperc = load(mutfoldf);

% Calculate likelihoods

edgeWeights = squeeze(edgeWeights);

score_edge = zeros(1,size(edgeWeights,2));
for i=1:size(edgeWeights,2)
    score_edge(i) = sum(edgeWeights(:,i).^2);
end

[sorted_score_edge indices] = sort(score_edge,'descend');

naa = 21;

nodeWeights = squeeze(nodeWeights);
naasq = naa*naa;

% For the gap characters
if (size(edgeWeights,1) ~= 1)
    edgeWeights(naasq,:) = 0;
else
    edgeWeights(1,naasq) = 0;
end
nodeWeights(naa,:) = 0;
score_plot = zeros(length(testNdx),1);
score_plot_node = zeros(length(testNdx),1);
score_plot_edge = zeros(length(testNdx),1);

for i = testNdx
    score_edge = 0;
    score_node = 0;
    b = nmsa(i,:);
        %b = b+1 ?
    
    for j = 1:size(edgeStruct.edgeEnds,1)
        a = edgeStruct.edgeEnds(j,:);
        if (size(edgeWeights,1) ~= 1)
            ew = reshape(edgeWeights(:,j),naa,naa);
        else
            ew = reshape(edgeWeights,naa,naa);
        end
        score_edge = score_edge + (ew(b(a(1)),b(a(2)))) ;
    end
     for j = 1:size(nodeWeights,2)
         score_node = score_node + (nodeWeights(b(j),j));
     end
    score_plot(i) = score_edge + score_node;
    score_plot_node(i) = score_node;
    score_plot_edge(i) = score_edge;
end

[s,sidx] = sort(score_plot);
[R,P] = corrcoef(score_plot,foldperc);

if(do_display)
    figure;
    scatter(score_plot,foldperc)
    saveas(gcf,'scatter_oz.png','png')
    fprintf('Corr coeff between gremlin and foldperc %f\n',R(2,1))
%     plot(score_plot);
%     plot(s);
%     hold on;
%     plot(foldperc(sidx),'g');
    %plot(score_plot_edge,'g');
    %plot(score_plot_node,'r');
end;