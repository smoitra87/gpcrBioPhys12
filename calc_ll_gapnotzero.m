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

DO_SECTION = 1;


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

% % For the gap characters
% if (size(edgeWeights,1) ~= 1)
%     edgeWeights(naasq,:) = 0;
% else
%     edgeWeights(1,naasq) = 0;
% end
% nodeWeights(naa,:) = 0;

score_plot = zeros(length(testNdx),1);
score_plot_node = zeros(length(testNdx),1);
score_plot_edge = zeros(length(testNdx),1);


if DO_SECTION == 1

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
        xlabel('Likelihood');
        ylabel('Folding perc');
        title('Scatter plot Likelihood vs Folding perc');           
        saveas(gcf,'scatter_oz.png','png')
        fprintf('Corr coeff between gremlin and foldperc %f\n',R(2,1))
    %     plot(score_plot);
    %     plot(s);
    %     hold on;
    %     plot(foldperc(sidx),'g');
        %plot(score_plot_edge,'g');
        %plot(score_plot_node,'r');
    end;


    % Compare the mutant score with natural bovine rhodopsin
    natf = '1F88.fasta';
    natseq = fastaread(natf);
    natmsa = {natseq(1).Sequence};
    nnat = converttonumericmsa(natmsa);
    score_edge = 0;
    score_node = 0;
    for j = 1:size(edgeStruct.edgeEnds,1)
            
            b = nnat;
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
    score_nat = score_node+ score_edge;

    figure;
    plot(s,'ro-');
    xlabel('Mutant Seq id');
    ylabel('Likelihood score');
    hold on;
    line([0,length(mut_seqs)],[score_nat,score_nat]);
    title('Plot Mutant Seq likelihood vs natural');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some post processing 

    % Find the mutants with the lowest score

    % 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now for HGMD data

if DO_SECTION == 1

    mutf = 'hgmd.fasta';
    mut_seqs = fastaread(mutf);
    trainf = 'GPCR_ranga.fasta';
    train_seqs = fastaread(trainf);
    %train_seqs = train_seqs(randperm(length(mut_seqs)));


    train_msa = cell(length(train_seqs),1);
    mut_msa = cell(length(mut_seqs),1);

    testNdx= 1:length(mut_seqs);
    testtrdx = 1:length(train_msa);

    for seqid = 1:length(mut_seqs) 
        mseq = mut_seqs(seqid).Sequence;
        mut_msa{seqid} = mseq;
    end
    nmsa = converttonumericmsa(mut_msa);


    for seqid = 1:length(train_seqs) 
        tseq = train_seqs(seqid).Sequence;
        train_msa{seqid} = tseq;
    end
    ntmsa = converttonumericmsa(train_msa);


    % Calculate likelihoods

    score_plot = zeros(length(testNdx),1);
    score_plot_node = zeros(length(testNdx),1);
    score_plot_edge = zeros(length(testNdx),1);

    score_train_plot = zeros(length(testtrdx),1);
    score_train_plot_node = zeros(length(testtrdx),1);
    score_train_plot_edge = zeros(length(testtrdx),1);


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

    % For training data

    e_gap = zeros(testtrdx,1) ;
    n_gap = zeros(testtrdx,1) ;

    for i = testtrdx
        score_edge = 0;
        score_node = 0;
        b = ntmsa(i,:);
            %b = b+1 ?
        e_gapcount = 0;
        n_gapcount = 0;

        for j = 1:size(edgeStruct.edgeEnds,1)
            a = edgeStruct.edgeEnds(j,:);
            if (size(edgeWeights,1) ~= 1)
                ew = reshape(edgeWeights(:,j),naa,naa);
            else
                ew = reshape(edgeWeights,naa,naa);
            end
            score_edge = score_edge + (ew(b(a(1)),b(a(2)))) ;
            if b(a(1)) == 21 || b(a(2)) == 21
                e_gapcount = e_gapcount + 1;
            end

        end
         for j = 1:size(nodeWeights,2)
             score_node = score_node + (nodeWeights(b(j),j));
             if b(j) == 21 
                  n_gapcount = n_gapcount + 1;
             end
         end
    %     score_edge = score_edge * (e_gapcount + size(edgeStruct.edgeEnds,1))/ ...
    %         size(edgeStruct.edgeEnds,1) ;
    %     score_node = score_node * (n_gapcount + size(nodeWeights,2))/ ...
    %         size(nodeWeights,2) ;
        e_gap(i) = e_gapcount;
        n_gap(i) = n_gapcount;

        score_train_plot(i) = score_edge + score_node;
        score_train_plot_node(i) = score_node;
        score_train_plot_edge(i) = score_edge;
    end

    [s,sidx] = sort(score_plot);

    if(do_display)
        figure;
        %saveas(gcf,'scatter_oz.png','png')
        %fprintf('Corr coeff between gremlin and foldperc %f\n',R(2,1))
        plot(sort(score_plot));
        hold on;
        plot(sort(score_train_plot),'g');
        xlabel('sequence id');
        ylabel('lieklihood gremlin');
        title('Plottign likelihood vs sequence id');
        legend('Sorted Mutant seq','Sorted Blast Seq');
        %plot(score_plot_edge,'g');
        %plot(score_plot_node,'r');
    end;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some post processing 

% Plot mutant score vs natural bovine score
figure;
plot(sort(score_plot),'ro-');
title('Comparing HGMD mutants with nat bov rhodopsin')
xlabel('Mutant Seq id');
ylabel('Likelihood score');
hold on;
line([0,length(mut_seqs)],[score_nat,score_nat]);


% Find the mutants with the lowest score and highest score
[s,idx] = sort(score_plot);
fprintf('\n--------Lowest scoring mutants---------\n');
for i = 1:10 
   fprintf('Mutant seq %s\n',mut_seqs(idx(i)).Header)
end
[s,idx] = sort(score_plot,'descend');
fprintf('\n--------Highest scoring mutants---------\n');
for i = 1:10 
   fprintf('Mutant seq %s\n',mut_seqs(idx(i)).Header)
end

% Plot the natscores vs blast score on horiz line
figure;
plot(sort(score_train_plot),'ro-');
title('Comparing Blast with nat bov rhodopsin')
xlabel('Blast Seq id');
ylabel('Likelihood score');
hold on;
line([0,length(train_seqs)],[score_nat,score_nat]);


% Find the blast seqs with the lowest score and highest score
[s,idx] = sort(score_train_plot);
fprintf('\n--------Lowest scoring Blast seqs---------\n');
for i = 1:10 
   fprintf('Mutant seq %s\n',train_seqs(idx(i)).Header)
end
[s,idx] = sort(score_train_plot,'descend');
fprintf('\n--------Highest scoring Blast seqs---------\n');
for i = 1:10 
   fprintf('Mutant seq %s\n',train_seqs(idx(i)).Header)
end





