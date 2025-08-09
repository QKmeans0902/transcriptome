clear;clc;

% take AAL atlas for example
num_regions = 116

% where your analysis was done
base_dir = 'E:\Data\Analysis\AHBA';

% load gene expression data (from abagen toolbox)
gene_exp_path = [base_dir, '\AAL116_AHBA.csv'];
gene_exp = csvread(gene_exp_path, 1, 1);

% get brain network topological measure (take degree for example)
degree_path = [base_dir, '\degree.csv'];
degree = csvread(degree_path, 1, 1);

% get group label and other covariates
var_path = [base_dir, '\var.csv'];
var = csvread(var_path, 1, 1);
group = categorical(var(:,1)); % 1 for controls / 2 for patients
age = var(:,2);
sex = var(:,3);

for i = 1:num_regions
    region_degree = degree(:,i)
    tbl = table(region_degree, group, age, sex, 'Variable',{'measure','group','sex','age'}) 
    lm = fitlm(tbl,'measure~ 1+group+sex+age');
    stat = anova(lm);
    p_val(i,:) = lm.Coefficients{2,4};
    t_val(i,:) = lm.Coefficients{2,3};
end

stat_path = [base_dir, '\degree_tval.txt'];
save(stat_path, 't_val', '-ascii');

% exclude regions without AHBA expression values
ex_reg = find(mean(gene_exp, 2) == 0);
gene_exp(ex_reg, :) = [];
t_val(ex_reg, :) = [];

% z-score:
X = zscore(gene_exp);
Y = zscore(t_val);

% perform full PLS and plot variance in Y explained by top 10 components
[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X, Y);
dim = 10;
plot(1:dim, 100*PCTVAR(2, 1:dim), '-o', 'LineWidth', 1.5, 'Color', [140/255,0,0]);
set(gca, 'Fontsize',14)
xlabel('PLS components', 'FontSize', 14);
ylabel('Percent Variance Explained', 'FontSize', 14);
grid on

% decide number of components you want to keep and calculate again to keep data tidy
dim = 2;
[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X, Y, dim);
Rsq = cumsum(100*PCTVAR(2, 1:dim));

% align PLS components with desired direction%
[R1, P1] = corr(XS(:, 1), t_val);
if R1(1, 1) < 0
    stats.W(:, 1) = -1 * stats.W(:, 1);
    XS(:, 1) = -1 * XS(:, 1);
end
[R2, p2] = corr(XS(:, 2), t_val);
if R2(1, 1) < 0
    stats.W(:, 2) = -1 * stats.W(:, 2);
    XS(:, 2)= -1 * XS(:, 2);
end

% correlation between PLS components and with t-val:
for z = 1:dim
    [r(z, :), ~] = corr(XS(:, z), t_val);
end

% spin test for PLS regression results (surrogate maps were obtained from BrainSMASH)
surrogate_maps = load([base_dir, '\degree_tval_surrogate_maps.txt']);

for j = 1:size(surrogate_maps, 1)
    surrogate_map = surrogate_maps(j,:)';
    Y_spin = zscore(surrogate_map);
    [XL_spin, YL_spin, XS_spin, YS_spin, BETA_spin, PCTVAR_spin, MSE_spin, stats_spin] = plsregress(X, Y_spin, dim);
    PCTVAR_spin_all(j,:) = PCTVAR_spin(2,:);
    Rsq_spin_all(j, :) = cumsum(100*PCTVAR_spin(2, 1:dim));
    for z = 1:dim
        [r_spin_all(j, z), ~] = corr(XS(:,z), surrogate_map);
        [r_spin_all(j, z), ~] = corr(XS(:,z), surrogate_map);
    end
end

for z = 1:dim
    p_single(z, :) = length(find(PCTVAR_spin_all(:, z) >= PCTVAR(2, z)))/size(surrogate_maps, 1);
    p_cumsum(z, :) = length(find(Rsq_spin_all(:, z) >= Rsq(:, z)))/size(surrogate_maps, 1);
    p_corr(z, :) = length(find(r_spin_all(:, z) >= r(:, z)))/size(surrogate_maps, 1);
end

% get gene names
gene_exp_info = readtable(gene_exp_path);
gene_name = gene_exp_info.Properties.VariableNames';
gene_name(1) = []; % remove the first VaribleName that is not the gene name
gene_index = 1:length(gene_name);

% calculate gene weights
bootnum = 1000;

[PLS1w, x1] = sort(stats.W(:, 1), 'descend');
PLS1ids = gene_name(x1);
geneindex1 = gene_index(x1);
csvwrite('PLS1_ROIscores.csv', XS(:,1));
PLS1_score = XS(:,1);

[PLS2w, x2] = sort(stats.W(:, 2),'descend');
PLS2ids = gene_name(x2);
geneindex2 = geneindex(x2);
csvwrite('PLS2_ROIscores.csv', XS(:, 2));
PLS2_score = XS(:, 2);

% define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

% start bootstrap
for i = 1:bootnum
    myresample = randsample(size(X, 1), size(X, 1));
    res(i, :) = myresample;
    Xr = X(myresample, :);
    Yr = Y(myresample, :);
    [XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(Xr, Yr, dim);
      
    temp = stats.W(:, 1); % extract PLS1 weights
    newW = temp(x1); % order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w, newW) < 0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW = -1 * newW;
    end
    PLS1weights = [PLS1weights, newW]; % store (ordered) weights from this bootstrap run

    temp = stats.W(:, 2); % extract PLS2 weights
    newW = temp(x2);
    if corr(PLS2w, newW) < 0
        newW = -1 * newW;
    end
     PLS2weights = [PLS2weights, newW];
end

% get standard deviation of weights from bootstrap runs
PLS1sw = std(PLS1weights');
PLS2sw = std(PLS2weights');

% get bootstrap weights
temp1 = PLS1w./PLS1sw';
temp2 = PLS2w./PLS2sw';

% order bootstrap weights (Z) and names of regions
[Z1, ind1] = sort(temp1, 'descend');
[Z2, ind2] = sort(temp2, 'descend');
PLS1 = PLS1ids(ind1);
PLS2 = PLS2ids(ind2);
geneindex1 = geneindex1(ind1);
geneindex2 = geneindex2(ind2);

% print out results
fid1 = fopen([base_dir, 'PLS1_geneWeights.csv','w');
for i = 1:length(gene_name)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
end
fclose(fid1)

fid2 = fopen([base_dir, 'PLS1_geneWeights.csv','w');
for i = 1:length(gene_name)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i}, geneindex2(i), Z2(i));
end
fclose(fid2)
