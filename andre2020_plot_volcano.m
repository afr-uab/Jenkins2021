% Jenkins et al. 2021
%     "Lung dendritic cells migrate to the spleen to prime long-lived 
%     memory CD8+ T cell precursors after influenza virus infection"
%
% Generates volcano plot (Figure 5c)
%
% Alex Rosenberg
% University of Alabama at Birmingham
% 8/6/2020


% load data
andre2020_load_data  

% load gene lists
L(1).name = 'Characteristic Memory Genes';
L(1).genes = readcell('andre_gene_list_CHARACTERISTIC_MEMORY_GENES_v2.txt');
L(2).name = 'Inhibitory Receptor Genes';
L(2).genes = readcell('andre_gene_list_INHIBITORY_RECEPTOR_GENES.txt');

% comparisons
D(1).ylab = 'log_2 Fold Change (Activated Spleen over Activated MLN)';
D(2).ylab = 'log_2 Fold Change (Activated Spleen over Spleen Naive)';
D(3).ylab = 'log_2 Fold Change (Activated MLN over Spleen Naive)';

% x and y data
D(1).x = A2.feat.SA_over_LA_paired_logfc;
D(2).x = A2.feat.SA_over_SN_logfc;
D(3).x = A2.feat.LA_over_SN_logfc;
D(1).y = -log10(A2.feat.SA_over_LA_paired_pval);
D(2).y = -log10(A2.feat.SA_over_SN_pval);
D(3).y = -log10(A2.feat.LA_over_SN_pval);
D(1).q = A2.feat.SA_over_LA_paired_fdr;
D(2).q = A2.feat.SA_over_SN_fdr;
D(3).q = A2.feat.LA_over_SN_fdr;

% limits
D(1).xylim = [10 9];
D(2).xylim = [12 13];
D(3).xylim = [12 13]; 

% -log10 p-value corresponding to fdr < .05
D(1).ythr = min(D(1).y(D(1).q < .05));
D(2).ythr = min(D(2).y(D(2).q < .05));
D(3).ythr = min(D(3).y(D(3).q < .05));

% x threshold
D(1).xthr = log2(1.5);
D(2).xthr = 1;
D(3).xthr = 1;

% called out genes
i1 = nan(length(L(1).genes), 1);
for i = 1:length(i1), i1(i) = find(strcmp(L(1).genes{i}, A2.feat.gene)); end
i2 = nan(length(L(2).genes), 1);
for i = 1:length(i2), i2(i) = find(strcmp(L(2).genes{i}, A2.feat.gene)); end

% figure layout
pw = 650;
ph = 500;
left = 65;
bot = 65;
right = 10;
top = 10;
fw = left + pw + right;
fh = bot + ph + top;
pos = [left / fw, bot / fh, pw / fw, ph / fh];
font = struct('fontweight', 'bold', 'fontname', 'arial', 'fontsize', 18);

% plot for each comparison
for k = 1:3

    % draw axes and plot
    figure('position', [20 20 fw fh]);
    subplot('position', pos, 'nextplot', 'add');
    n = length(D(k).x);
    if k == 1
        plot(D(k).x, D(k).y, '.', 'color', [1 1 1] * .55, 'markersize', 8);
    else
        iup = intersect(find(D(k).x >  D(k).xthr), find(D(k).q < .05));
        idn = intersect(find(D(k).x < -D(k).xthr), find(D(k).q < .05));
        ino = setdiff((1:n), [iup; idn]);
        plot(D(k).x(ino), D(k).y(ino), '.', 'color', [1 1 1] * .55, 'markersize', 8);
        plot(D(k).x(iup), D(k).y(iup), '.', 'color', 'r', 'markersize', 8);
        plot(D(k).x(idn), D(k).y(idn), '.', 'color', 'b', 'markersize', 8); 
        text(D(k).xthr *  1.2, .95 * D(k).xylim(2), num2str(length(iup)), font, 'color', 'r');
        text(D(k).xthr * -1.2, .95 * D(k).xylim(2), num2str(length(idn)), font, 'color', 'b', 'horizontalalignment', 'right');
    end
    plot([1 1] * D(k).xthr, [0 D(k).xylim(2)], 'k-');
    plot([1 1] * -1 * D(k).xthr, [0 D(k).xylim(2)], 'k-');
    plot(D(k).xylim(1) * [-1 1], [1 1] * D(k).ythr, 'k--');
    if k == 1
        plot(D(k).x(i1), D(k).y(i1), 'ro', 'linewidth', 2, 'markersize', 8);
        plot(D(k).x(i2), D(k).y(i2), 'bo', 'linewidth', 2, 'markersize', 8);
        text(D(k).x(i1) + .25, D(k).y(i1), A2.feat.gene(i1), 'color', 'r', font);
        text(D(k).x(i2) - .25, D(k).y(i2), A2.feat.gene(i2), 'color', 'b', font, 'horizontalalignment', 'right');
    end
    set(gca,...
        'linewidth', 1,...
        'ticklength', [0 0],...
        'xlim', D(k).xylim(1) * [-1 1],...
        'ylim', [0 D(k).xylim(2)],...
        'fontname', 'arial',...
        'fontsize', 14);
    hx = xlabel(D(k).ylab,...
        'fontname', 'arial',...
        'fontsize', 20);
    hx.Position(2) = 0 - (.06667 * diff(ylim));
    hy = ylabel('-log_1_0 p-Value',...
        'fontname', 'arial',...
        'fontsize', 20);
    hy.Position(1) = min(xlim) - (.045 * diff(xlim));
    if k == 1
        text(-9.5, 8.85, L(1).name, 'fontname', 'arial', 'fontsize', 20, 'color', 'r');
        text(-9.5, 8.4, L(2).name, 'fontname', 'arial', 'fontsize', 20, 'color', 'b');
    end
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');         

end