function [ret] = auc_comparisons(numeric_type,protein_name,fname,display)

load(fname,'adjFinal','Xedge', 'edgeWeights', 'Xnode', 'nodeWeights','edgeStruct','infoStruct');
edgeWeights = squeeze(edgeWeights);

score_edge = zeros(1,size(edgeWeights,2));
if (numeric_type == 2)
    for i=1:size(edgeWeights,2)
        score_edge(i) = sum(edgeWeights(:,i).^2);
    end
end
if (numeric_type == 3)
    for i=1:size(edgeWeights,2)
        score_edge(i) = max(edgeWeights(:,i));
        if (score_edge(i) < 0)
            score_edge(i) = 0;
        end
    end
end
pname = {'1b1k','1befa','1bhp','1bk7a','1bm8','1btn','1bx7','1c25','1cc7a','1cdza','1ctf','1f94a'};
cbk_edges = [16 112 26 74 8 2 2 12 2 7 58 17];

[sorted_score_edge indices] = sort(score_edge,'descend');
indices(1:3)
for i=1:length(pname)
    if(strcmp(pname{i},protein_name) == 1)
        save_index = i;
    end
end

number_of_edges = cbk_edges(save_index);

if (number_of_edges >= size(edgeWeights,2))
    fprintf('CBK has more edges, answers same as before');
    return;
end

fprintf('Considering %d edges only\n',number_of_edges);

trainf = strcat(protein_name,'_train.tp');
trainf = strcat('./more_folds/',trainf);
testf = strcat(protein_name,'_test.tp');
testf = strcat('./more_folds/',testf);
addpath(genpath(pwd));

naa = 21;

nodeWeights = squeeze(nodeWeights);
naasq = naa*naa;

Ba = dlmread(trainf);
train_sequences = Ba(:,1:end-1);

Aa = dlmread(testf);
test_sequences = Aa(:,1:end-1);


testNdx = 1:(size(Ba,1) + size(Aa,1));
%testNdx
if (size(edgeWeights,1) ~= 1)
    edgeWeights(naasq,:) = 0;
else
    edgeWeights(1,naasq) = 0;
end
nodeWeights(naa,:) = 0;
score_plot = zeros(length(testNdx),1);

size_train = size(train_sequences,1);

for i = testNdx
    score_edge = 0;
    score_node = 0;
    if (i > size_train)
        b = test_sequences(i-size_train,:);
        b = b+1;
    else
        b = train_sequences(i,:);
        b = b+1;
    end
    for j = 1:number_of_edges
        a = edgeStruct.edgeEnds(indices(j),:);
        if (size(edgeWeights,1) ~= 1)
            ew = reshape(edgeWeights(:,indices(j)),naa,naa);
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
if(display)
    figure;
    plot(score_plot);
    hold on
    plot(score_plot_edge,'g');
    plot(score_plot_node,'r');
end;
number_of_train_seq = {339,28,31,120,176,139,55,144,624,66,195,144};
number_of_pos_seq = {509,42,48,179,266,208,82,217,936,99,292,215};
    [x,auc1,y,z] = roc(score_plot(number_of_train_seq{save_index}:number_of_pos_seq{save_index}),score_plot(number_of_pos_seq{save_index}+1:end));
    [x,auc2,y,z] = roc(score_plot_edge(number_of_train_seq{save_index}:number_of_pos_seq{save_index}),score_plot_edge(number_of_pos_seq{save_index}+1:end));
    [x,auc3,y,z] = roc(score_plot_node(number_of_train_seq{save_index}:number_of_pos_seq{save_index}),score_plot_node(number_of_pos_seq{save_index}+1:end));
    fprintf('Actual AUC with CBK number of edges: %d\n',auc1);
    fprintf('Edge only AUC with CBK number of edges: %d\n',auc2);
    fprintf('Node only AUC with CBK number of edges: %d\n',auc3);
    
    
ret(1) = auc1;
ret(2) = auc2;
ret(3) = auc3;

















clearvars -except ret numeric_type protein_name fname display;



load(fname,'adjFinal','Xedge', 'edgeWeights', 'Xnode', 'nodeWeights','edgeStruct','infoStruct');
edgeWeights = squeeze(edgeWeights);

pname = {'1b1k','1befa','1bhp','1bk7a','1bm8','1btn','1bx7','1c25','1cc7a','1cdza','1ctf','1f94a'};
cbk_edges = [16 112 26 74 8 2 2 12 2 7 58 17];

for i=1:length(pname)
    if(strcmp(pname{i},protein_name) == 1)
        save_index = i;
    end
end

number_of_edges = cbk_edges(save_index);
cbk_edge_fname = strcat('./cbk_results/',protein_name);
cbk_edge_fname = strcat(cbk_edge_fname,'_edges.out');
actual_edges = dlmread(cbk_edge_fname);

fprintf('Considering %d edges only\n',number_of_edges);

trainf = strcat(protein_name,'_train.tp');
trainf = strcat('./more_folds/',trainf);
testf = strcat(protein_name,'_test.tp');
testf = strcat('./more_folds/',testf);
addpath(genpath(pwd));

naa = 21;


nodeWeights = squeeze(nodeWeights);
naasq = naa*naa;

Ba = dlmread(trainf);
train_sequences = Ba(:,1:end-1);

Aa = dlmread(testf);
test_sequences = Aa(:,1:end-1);


testNdx = 1:(size(Ba,1) + size(Aa,1));
%testNdx
if (size(edgeWeights,1) ~= 1)
    edgeWeights(naasq,:) = 0;
else
    edgeWeights(1,naasq) = 0;
end
nodeWeights(naa,:) = 0;
score_plot = zeros(length(testNdx),1);

size_train = size(train_sequences,1);
if (size(edgeWeights,1) ~= 1)
    numEdges = size(edgeWeights,2);
else
    numEdges = 1;
end
ew = zeros(number_of_edges,naa,naa);
    for j = 1:number_of_edges
        for k = 1:numEdges
            temp = edgeStruct.edgeEnds(k,:);
            if ((temp(:,1) == actual_edges(j,1))&&(temp(:,2) == actual_edges(j,2)))
                if (size(edgeWeights,1) ~= 1)
                    ew(j,:,:) = reshape(edgeWeights(:,k),naa,naa);
                else
                    ew(j,:,:) = reshape(edgeWeights,naa,naa);
                end
            end
        end
    end


    
for i = testNdx
    score_edge = 0;
    score_node = 0;
    if (i > size_train)
        b = test_sequences(i-size_train,:);
        b = b+1;
    else
        b = train_sequences(i,:);
        b = b+1;
    end
    for j = 1:number_of_edges
        if (sum(sum((ew(j,:,:).^2))) == 0)
            if (i == 1)
                fprintf('This edge was not found: %d %d\n',actual_edges(j,1),actual_edges(j,2));
            end
        else
            score_edge = score_edge + (ew(j,b(actual_edges(j,1)),b(actual_edges(j,2))));
        end
    end
     for j = 1:size(nodeWeights,2)
         score_node = score_node + (nodeWeights(b(j),j));
     end
    score_plot(i) = score_edge + score_node;
    score_plot_node(i) = score_node;
    score_plot_edge(i) = score_edge;
end
if(display)
    figure;
    plot(score_plot);
    hold on
    plot(score_plot_edge,'g');
    plot(score_plot_node,'r');
end;
number_of_train_seq = {339,28,31,120,176,139,55,144,624,66,195,144};
number_of_pos_seq = {509,42,48,179,266,208,82,217,936,99,292,215};
    [x,auc1,y,z] = roc(score_plot(number_of_train_seq{save_index}:number_of_pos_seq{save_index}),score_plot(number_of_pos_seq{save_index}+1:end));
    [x,auc2,y,z] = roc(score_plot_edge(number_of_train_seq{save_index}:number_of_pos_seq{save_index}),score_plot_edge(number_of_pos_seq{save_index}+1:end));
    [x,auc3,y,z] = roc(score_plot_node(number_of_train_seq{save_index}:number_of_pos_seq{save_index}),score_plot_node(number_of_pos_seq{save_index}+1:end));
    fprintf('Actual AUC with CBK edges: %d\n',auc1);
    fprintf('Edge only AUC with CBK edges: %d\n',auc2);
    fprintf('Node only AUC with CBK edges: %d\n',auc3);

ret(4) = auc1;
ret(5) = auc2;
ret(6) = auc3;











