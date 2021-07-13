% Jenkins et al. 2021
%     "Lung dendritic cells migrate to the spleen to prime long-lived 
%     memory CD8+ T cell precursors after influenza virus infection"
%
% Generates PCA plots (Figure 5b)
%
% Alex Rosenberg
% University of Alabama at Birmingham
% 8/6/2020


% load data
andre2020_load_data  

% perform PCA on CPM data
[coeff, score, ~, ~, pe] = pca(A2.cpm');
x = score(:, 1);
y = score(:, 2);

% sample groups, symbols, colors
G = {...
    'LA', 's', [0.0 0.7 0.0], 12, find(strcmp(A2.samp.group, 'LA'));...
    'SA', 'v', [0.7 0.0 0.7], 12, find(strcmp(A2.samp.group, 'SA'));...
    'SN', 'o', [0.5 0.5 0.5], 12, find(strcmp(A2.samp.group, 'SN'))};

% plot setup
pw = 300;
ph = 300;
left = 60;
bot = 60;
right = 15;
top = 10;
fw = left + pw + right;
fh = bot + ph + top;
r = 150;
pos = [left / fw, bot / fh, pw / fw, ph / fh];
figure('position', [20 20 fw fh]);

% draw pca plot
subplot('position', pos, 'nextplot', 'add');
plot([0 0], [-r r], '-', 'color', [1 1 1] * .3, 'linewidth', .5);
plot([-r r], [0 0], '-', 'color', [1 1 1] * .3, 'linewidth', .5);
h = nan(3, 1);
for i = 1:3
    h(i) = plot(x(G{i, 5}), y(G{i, 5}), G{i, 2},...
        'color', G{i, 3},...
        'markerfacecolor', G{i, 3},...
        'markersize', G{i, 4});
end
set(gca,...
    'xlim', [-r r],...
    'ylim', [-r r],...
    'ticklength', [0 0],...
    'linewidth', 1,...
    'box', 'on',...
    'fontname', 'arial',...
    'fontsize', 12);
set(gcf, 'inverthardcopy', 'off', 'color', 'w');  
hx = xlabel(['PC1 (' sprintf('%.1f', pe(1)) '%)'], 'fontname', 'arial', 'fontsize', 16);
hy = ylabel(['PC2 (' sprintf('%.1f', pe(2)) '%)'], 'fontname', 'arial', 'fontsize', 16);
hx.Position(2) = min(ylim) - (.11 * diff(ylim));
hy.Position(1) = min(xlim) - (.12 * diff(xlim));
legend(h, G(:, 1), 'location', 'northeast', 'fontname', 'arial', 'fontsize', 16, 'box', 'off');
set(gcf, 'inverthardcopy', 'off', 'color', 'w');         
