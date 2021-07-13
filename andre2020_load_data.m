% Jenkins et al. 2021
%     "Lung dendritic cells migrate to the spleen to prime long-lived 
%     memory CD8+ T cell precursors after influenza virus infection"
%
% Loads data from edgeR output into data structure to be used for other
% figure-generating scripts.
%
% Alex Rosenberg
% University of Alabama at Birmingham
% 8/6/2020


% path to analysis output *** NOTE - CHANGE TO REFLECT YOUR PATH TO DATA
pt = '/Users/afr/home/projects/BALLESTEROS-TATO/2018_07_RNAseq_spleen_MLN_naive/andre_2020_reprocessed/edger/';

% edger output files - col1 = file name, col2 = field name prefix
sfile = {...
    'andre2020_DEG_LA_over_SN.txt',                  'LA_over_SN';...
    'andre2020_DEG_SA_over_SN.txt',                  'SA_over_SN';...
    'andre2020_DEG_SA_over_LA.txt',                  'SA_over_LA';...
    'andre2020_DEG_SA_over_LA_activated_paired.txt', 'SA_over_LA_paired'};

% edger cpm file
cpmfile = 'andre2020_log2cpm.txt';

% edger rpkm file
rpkmfile = 'andre2020_log2rpkm.txt';

% read in stats files
cname = {'row', 'gene', 'length', 'logFC', 'logCPM', 'F', 'pValue', 'FDR'};
S = struct;
for i = 1:4
    S.(sfile{i, 2}) = readtable([pt sfile{i, 1}], 'NumHeaderLines', 1);
    S.(sfile{i, 2}).Properties.VariableNames = cname;
    if i > 1
        if ~isequal(S.(sfile{i, 2}).gene, S.(sfile{1, 2}).gene)
            error('gene order mismatch');
        end
    end
end

% read in expression files as cell arrays
e1 = readcell([pt cpmfile]);
e2 = readcell([pt rpkmfile]);

% get header rows
h1 = e1(1, 1:(end - 1))';
h2 = e2(1, 1:(end - 1))';
if ~isequal(h1, h2), error('cpm, rpkm col headings different'); end

% convert sample names into sample table
h1 = regexprep(h1, [h1{1}(1:max(strfind(h1{1}, '/'))) '|.count'], '');
h1 = regexprep(h1,...
    {'MLNode_activ', 'Spleen_activ', 'Spleen_naive'},...
    {'LA', 'SA', 'SN'});
temp = cell(length(h1), 3);
for i = 1:length(h1), temp(i, :) = strsplit(h1{i}, '_'); end
samp = cell2table(temp, 'VariableNames', {'group', 'rep', 'pool'});

% gene names
g1 = e1(2:end, 1);
g2 = e2(2:end, 1);
if ~isequal(g1, g2), error('cpm, rpkm gene names different'); end

% make feature table with gene names, stats
feat = table;
feat.gene = g1;
feat.length = S.(sfile{1, 2}).length;
for i = 1:4
    k = sfile{i, 2};
    feat.([k '_logfc']) = S.(k).logFC;
    feat.([k '_pval']) = S.(k).pValue;
    feat.([k '_fdr']) = S.(k).FDR;
end

% SA over LA score
feat.SA_over_LA_paired_score = ...
    -log10(feat.SA_over_LA_paired_pval) .* ...
    sign(feat.SA_over_LA_paired_logfc);

% master dataset container
A2.samp = samp;
A2.feat = feat;
A2.rpkm = cell2mat(e2(2:end, 2:end));
A2.cpm  = cell2mat(e1(2:end, 2:end));

% clean up
clearvars -except A2







