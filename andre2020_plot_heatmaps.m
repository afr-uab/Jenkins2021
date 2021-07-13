% Jenkins et al. 2021
%     "Lung dendritic cells migrate to the spleen to prime long-lived 
%     memory CD8+ T cell precursors after influenza virus infection"
%
% Generates gene expression heat maps (Figures 5a and 5d)
%
% Alex Rosenberg
% University of Alabama at Birmingham
% 8/6/2020


% load data
andre2020_load_data  

% load gene lists
L(1).name = 'Characteristic Memory Genes';
L(1).genes = readcell('andre_gene_list_CHARACTERISTIC_MEMORY_GENES.txt');
L(2).name = 'Inhibitory Receptor Genes';
L(2).genes = readcell('andre_gene_list_INHIBITORY_RECEPTOR_GENES.txt');
L(3).name = 'Teff Differentiation Genes';
% L(3).genes = readcell('andre_gene_list_TEFF_DIFFERENTIATION_GENES.txt');
L(3).genes = readcell('andre_gene_list_TEFF_DIFFERENTIATION_GENES_v2.txt');

% get indices of gene lists
for i = 1:3
    for k = 1:length(L(i).genes)
        L(i).idx(k, 1) = find(strcmp(L(i).genes{k}, A2.feat.gene));
    end
end

% set up data structure for heat map 1
j = [find(strcmp(A2.samp.group, 'LA')); find(strcmp(A2.samp.group, 'SA'))];
F(1).gene = [L(1).genes; L(2).genes];
F(1).set = [...
    repmat({'Memory'}, length(L(1).genes), 1);...
    repmat({'Inhibitory Receptor'}, length(L(2).genes), 1)]; 
F(1).idx = [L(1).idx; L(2).idx];
F(1).cpm = A2.cpm(F(1).idx, j);
F(1).zscore = zscore(F(1).cpm, [], 2);
F(1).samp = cellstr([char(A2.samp.group(j)) char(A2.samp.rep(j))]);

% set up data structure for heat map 2
F(2).gene = L(3).genes;
F(2).set = repmat({'Teff'}, length(L(3).genes), 1);
F(2).idx = L(3).idx;
F(2).cpm = A2.cpm(F(2).idx, :);
F(2).zscore = zscore(F(2).cpm, [], 2);
F(2).samp = cellstr([char(A2.samp.group) char(A2.samp.rep)]);

% two heat maps
for ip = 1:2

    % size of data set
    [nf, ns] = size(F(ip).zscore);

    % plot configureation
    pix = 25;
    dw = 70;
    pw = pix * ns;
    ph = pix * nf;
    lw = 60;
    pw2 = 300;
    right = 10;
    sp1 = 25;
    ch = 25;
    sp2 = 50;
    top = 50;
    fw = dw + pw + lw;
    fh = sp1 + ch + sp2 + ph + top;
    posd = [0, (sp1 + ch + sp2) / fh, dw / fw, ph / fh];
    pos  = [dw / fw, (sp1 + ch + sp2) / fh, pw / fw, ph / fh];
    posl = [(dw + pw) / fw, (sp1 + ch + sp2) / fh, lw / fw, ph / fh];
    post = [dw / fw, (sp1 + ch + sp2 + ph) / fh, pw / fw, top / fh];
    if ns == 6, w = 1; else, w = .6667; end
    poss = [dw / fw, sp1 / fh, pw * w / fw, ch / fh];
    pref = struct(...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'xtick', [],...
        'ytick', [],...
        'box', 'off');
    crange = [-1 1] * 1.5;
    erange = [-2 12];
    figure('position', [20 20 fw fh]);

    % draw dendrogram and get sort order
    subplot('position', posd);
    z = linkage(pdist(F(ip).zscore, 'euclidean'), 'average');
    [hd, ~, isort] = dendrogram(z, 0, 'orientation', 'left');
    set(hd, 'color', [1 1 1] * .3, 'linewidth', 1);
    set(gca, pref, 'xlim', [0 max(xlim)], 'ylim', [.5 nf + .5], 'ydir', 'reverse');
    
    % gene names
    subplot('position', posl)
    gsort = F(ip).gene(isort);
    for i = 1:nf
        text(.1, i, gsort{i}, 'fontname', 'arial', 'fontsize', 16);
    end
    set(gca, pref, 'xlim', [0 1], 'ylim', [.5 nf + .5], 'ydir', 'reverse');
    
    % sample names
    subplot('position', post);
    for i = 1:ns
        text(i, .1, F(ip).samp{i}, 'fontname', 'arial', 'fontsize', 16, 'rotation', 90);
    end
    set(gca, pref, 'xlim', [.5 ns + .5], 'ylim', [0 1]);
    
    % draw heat map
    subplot('position', pos, 'nextplot', 'add');
    imagesc(F(ip).zscore(isort, :), crange);
    colormap(usa(256));
    set(gca, pref, 'xlim', [.5 ns + .5], 'ylim', [.5 nf + .5], 'ydir', 'reverse');
    plot(xlim, [1 1]' * (.5:1:(nf + .5)), '-', 'color', [1 1 1] * .3, 'linewidth', .5);
    plot([1 1]' * (.5:1:(ns + .5)), ylim, '-', 'color', [1 1 1] * .3, 'linewidth', .5);

    % draw color bar
    subplot('position', poss);
    imagesc(-100:100);
    colormap(usa);
    set(gca,...
        'ytick', [],...
        'ydir', 'normal',...
        'ticklength', [0 0],...
        'linewidth', .5,...
        'xtick', [min(xlim) mean(xlim) max(xlim)],...
        'xticklabel', {num2str(min(crange)), '0', num2str(max(crange))},...
        'fontname', 'arial',...
        'fontsize', 16);
    text(mean(xlim), 2, 'z-Score', 'fontname', 'arial', 'fontsize', 16, 'horizontalalignment', 'center');
   
    % clean up figure
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');
    
end
