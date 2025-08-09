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
group = categorical(var(:,1)); % 1:controls 2:patients
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
X=zscore(gene_exp);
Y=zscore(t_val);

%perform full PLS and plot variance in Y explained by top 10 components
[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X, Y);
dim = 10;
plot(1:dim, 100*PCTVAR(2,1:dim), '-o', 'LineWidth', 1.5, 'Color', [140/255,0,0]);
set(gca, 'Fontsize',14)
xlabel('PLS components', 'FontSize', 14);
ylabel('Percent Variance Explained', 'FontSize', 14);
grid on

% number of components you want to keep
dim=2;
[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats]=plsregress(X, Y, dim); % no need to do this but it keeps outputs tidy

% plot correlation of PLS component 1/2 with t-val:
[r1, p1]=corr(XS(:, 1), t_val);
[r2, p2]=corr(XS(:, 2), t_val);

% spin test
rep=10000;
% entropy_tval_surrogate_10000.txt are sourced from get_surrogate_map.ipynb
tval_surrogate_maps = load([root_path, '\4_Transcriptome\entropy_tval_surrogate_10000.txt']);

for j=1:rep
         tval_surrogate_map = tval_surrogate_maps(j,:)';
        [r_perm(j),p_perm(j)] = corr(XS(:,1), tval_surrogate_map);
end

P_corr = length(find(r_perm>=r))/rep;

for z = 1:size(t_val, 2)
    figure
    plot(XS(:, 1), t_val(:, z), 'r.')
    xlabel('XS scores for PLS component 1', 'FontSize',14);
    ylabel('t-statistic','FontSize',14);
    grid on
end

% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:

tval_surrogate_maps = load([root_path, '\4_Transcriptome\entropy_tval_surrogate_10000.txt']);
% tval_surrogate_maps = load([root_path, '\4_Transcriptome\entropy_tval_NoLhRh_surrogate_10000.txt']);
for dim=1:3
[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X, Y, dim);
temp = cumsum(100*PCTVAR(2, 1:dim));
Rsquared = temp(dim);
    for j=1:rep
         tval_surrogate_map = tval_surrogate_maps(j,:)';
        Yp = zscore(tval_surrogate_map);
        [XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X, Yp, dim);
        temp = cumsum(100*PCTVAR(2, 1:dim));
        Rsq(j) = temp(dim);
    end
R(dim) = Rsquared;
P(dim) = length(find(Rsq>=Rsquared))/rep;
end

figure
plot(1:dim, P,'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on

% regional PLS 1 (or PLS 2, PLS 3) used for brain mapping figure
PLS1_enorm = ones(n_reg, 1);
PLS1_enorm(ex_regID)=0;
PLS1_enorm(logical(PLS1_enorm)) = XS(:, 1); 
save([root_path, '\4_Transcriptome\entropy_PLS1.txt'],'PLS1_enorm','-ascii');

% keep data tidy
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,dim);

% get gene names
gene_exp_info = readtable(gene_exp_path);
gene_name = gene_exp_info.Properties.VariableNames';
gene_name(1) = []; % remove the first VaribleName that is not the gene name
gene_index = 1:length(gene_name);

% number of bootstrap iterations:
bootnum=1000;

% store regions' IDs and weights in descending order of weight for both components:
[R1,P1]=corr(XS(:,1), t_val);

% align PLS components with desired direction for interpretability 
if R1(1, 1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:, 1) = -1 * stats.W(:, 1);
    XS(:, 1) = -1 * XS(:, 1);
end

[PLS1w, x1] = sort(stats.W(:, 1), 'descend');
PLS1ids = gene_name(x1);
geneindex1 = gene_index(x1);

% print out results
csvwrite('PLS1_ROIscores.csv', XS(:,1));

% define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];

% start bootstrap
for i=1:bootnum
    myresample = randsample(size(X,1), size(X,1));
    res(i,:)=myresample; % store resampling out of interest
    Xr = X(myresample, :); % define X for resampled subjects
    Yr = Y(myresample, :); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp = stats.W(:,1); % extract PLS1 weights
    newW = temp(x1); % order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w, newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW]; % store (ordered) weights from this bootstrap run
end

% get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');

% get bootstrap weights
temp1=PLS1w./PLS1sw';

% order bootstrap weights (Z) and names of regions
[Z1, ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
pval_pos = 1-normcdf(Z1(Z1>0));
pval_neg = normcdf(Z1(Z1<0));
pval = [pval_pos;pval_neg];
fdr = [mafdr(pval_pos, 'BHFDR', 'true'); mafdr(pval_neg, 'BHFDR', 'true')];

% print out results
% later use first column of these csv files for pasting into GOrilla (for
% bootstrapped ordered list of genes) 
fid1 = fopen('PLS1_geneWeights_temp.csv','w');
for i=1:length(gene_name)
  fprintf(fid1,'%s, %d, %f, %f, %f\n', PLS1{i}, geneindex1(i), Z1(i), pval(i),fdr(i));
end
fclose(fid1)
