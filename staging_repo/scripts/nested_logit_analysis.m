function results = nested_logit_analysis()
%NESTED_LOGIT_ANALYSIS  Level-2 aggregate demand for sports-nutrition products.
%
%   Complete econometric pipeline for the applied-economics paper
%   "Price Elasticity and Brand Substitution in Sports Nutrition: A Nested
%   Logit Analysis of Protein Powders and Protein Bars".
%
%   1. Ingest the product catalogue and the weekly sales panel (plain CSV
%      reader, no toolboxes required - runs in MATLAB and GNU Octave).
%   2. Define the nests: protein POWDERS vs protein BARS, and define the
%      choice markets as flavour profiles (Mocha / Birthday Cake /
%      Salted Caramel are the focal markets).
%   3. Estimate the BLP Level-2 (aggregate) nested-logit model by Berry
%      inversion / two-step OLS:
%          ln(s_jt/s_0t) = x_jt*b + sigma*ln(s_j|g,t) + xi_jt
%   4. Compute own- and cross-price elasticities and print the elasticity
%      matrix (product level and brand level), also written to CSV.
%   5. Publication-quality figures: market shares, price/share scatter,
%      elasticity heat maps, substitution across price points, fit.
%   6. Validation block: (a) held-out promotional weeks of the weekly panel,
%      (b) held-out SKU subset of the catalogue, plus a benchmark model.
%
%   OUTPUT  results struct; CSV files in ../output/; figures in ../output/figures/
%   RUN     >> nested_logit_analysis

% =========================================================== 1. CONFIGURATION
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
cfg.catalogFile  = fullfile(root, 'data', 'sports_nutrition_supplements.csv');
cfg.panelFile    = fullfile(root, 'data', 'supplement_sales_weekly.csv');
cfg.outDir       = fullfile(root, 'output');
cfg.figDir       = fullfile(root, 'output', 'figures');
cfg.focalMarkets = {'Mocha / Coffee', 'Birthday Cake', 'Salted Caramel'};
cfg.promoThreshold   = 0.13;   % week is promotional if mean discount >= 13%
cfg.holdoutFraction  = 0.20;   % share of catalogue SKUs held out
if ~(exist(cfg.outDir, 'dir')), mkdir(cfg.outDir); end
if ~(exist(cfg.figDir, 'dir')), mkdir(cfg.figDir); end

fprintf('==================================================================\n');
fprintf(' NESTED LOGIT (BLP Level-2) - sports nutrition demand analysis\n');
fprintf('==================================================================\n\n');

% ============================================ 2. INGEST THE DATASET (CATALOGUE)
fprintf('[1] Ingesting product catalogue: %s\n', cfg.catalogFile);
C = read_csv(cfg.catalogFile);
hdr = C(1, :);
col = @(name) find(strcmp(hdr, name), 1);

nRaw = size(C, 1) - 1;
cat.category = C(2:end, col('Category'));
cat.name     = C(2:end, col('Name'));
cat.brand    = C(2:end, col('Brand'));
cat.flavor   = C(2:end, col('Flavor'));
cat.reviews  = to_number(C(2:end, col('Reviews')));
cat.price    = to_number(C(2:end, col('Price')));
cat.pps      = to_number(C(2:end, col('PricePerServing')));
% parent product identifier (brand | product name); built element-wise so that
% no whitespace is stripped from the separator
cat.prod = cell(nRaw, 1);
for k = 1:nRaw
    cat.prod{k} = [strtrim(cat.brand{k}), ' | ', strtrim(cat.name{k})];
end

ix = grp2id(cat.prod);
cat.nsku = zeros(nRaw, 1);
for k = 1:nRaw, cat.nsku(k) = sum(ix == ix(k)); end
% a parent product's cumulative review mass is split over its flavour/size
% SKUs so that a 20-flavour product does not receive 20x its review mass
cat.mass  = cat.reviews ./ cat.nsku;
cat.share = cat.mass / sum(cat.mass);
fprintf('      %d SKUs, %d distinct products, %d brands\n', nRaw, ...
        numel(unique(cat.prod)), numel(unique(cat.brand)));

% ---------------------------------------------------------- nests (2)
powderCats = {'WHEY PROTEIN', 'WHEY PROTEIN ISOLATE', 'MICELLAR CASEIN PROTEIN', ...
              'PLANT PROTEIN', 'LOW CARB PROTEIN'};
barCats = {'PROTEIN BARS'};
% The source catalogue lists several bar-form products (Quest Bars, Protein
% Crisp Bars, ...) under powder categories, so product-name form rules are
% applied inside the powder categories (must match reference_analysis.py).
cat.nest = repmat({''}, nRaw, 1);
for k = 1:nRaw
    nm = lower(strtrim(cat.name{k}));
    % substring form rules (Octave's regexp lacks \b word boundaries)
    isBarForm = ~isempty(strfind(nm, 'bar')) || ~isempty(strfind(nm, 'wafer')) || ...
                (~isempty(strfind(nm, 'cookie')) && isempty(strfind(nm, 'shake')));
    inPowderCat = any(strcmp(cat.category{k}, powderCats));
    if any(strcmp(cat.category{k}, barCats)) || (inPowderCat && isBarForm)
        cat.nest{k} = 'bar';
    elseif inPowderCat
        cat.nest{k} = 'powder';
    end
end
inside = ~cellfun(@isempty, cat.nest) & cat.mass > 0;
fprintf('      nests: %d powder SKUs vs %d bar SKUs\n', ...
        sum(strcmp(cat.nest(inside), 'powder')), sum(strcmp(cat.nest(inside), 'bar')));

% --------------------------------------- choice markets = flavour profiles
cat.market = cell(nRaw, 1);
for k = 1:nRaw, cat.market{k} = flavour_market(cat.flavor{k}); end

idx = find(inside);
D.market = cat.market(idx);  D.nest = cat.nest(idx);
D.brand   = cat.brand(idx);   D.prod = cat.prod(idx);
D.flavor  = cat.flavor(idx);
D.pps     = cat.pps(idx);     D.price = cat.price(idx);
D.share   = cat.share(idx);   D.n = numel(idx);
mkts = unique(D.market);

D.s0 = zeros(D.n, 1); D.sg = zeros(D.n, 1); D.sjg = zeros(D.n, 1);
for m = 1:D.n
    inM = strcmp(D.market, D.market{m});
    D.s0(m) = 1 - sum(D.share(inM));            % outside good = other flavours
    inG = inM & strcmp(D.nest, D.nest{m});
    D.sg(m)  = sum(D.share(inG));
    D.sjg(m) = D.share(m) / D.sg(m);
end
D.y = log(D.share ./ D.s0);
D.lnsjg = log(D.sjg);
fprintf('      %d estimation observations in %d flavour-profile markets\n', ...
        D.n, numel(mkts));
fprintf('      outside-good share ranges %.3f - %.3f\n\n', min(D.s0), max(D.s0));

% ================================================= 3. ESTIMATE THE NESTED LOGIT
fprintf('[2] Estimating nested logit (Berry inversion, two-step OLS)\n');
fprintf('    ln(s_jt/s_0t) = b*price_per_serving + sigma*ln(s_j|g,t) + market_FE + xi\n\n');
mnames = dummy_names(D.market);
X = [ones(D.n, 1), dummies(D.market), D.pps, D.lnsjg];
names = [{'const'}, mnames, {'pps', 'lns_jg'}];
[theta, V, stats] = ols_cluster(X, D.y, grp2id(D.prod));

iP = find(strcmp(names, 'pps'));
iS = find(strcmp(names, 'lns_jg'));
bPrice = theta(iP); sigma = theta(iS);
seP = sqrt(V(iP, iP)); seS = sqrt(V(iS, iS));
tP = bPrice/seP; tS = sigma/seS;
fprintf('    price per serving   b  = %8.4f (se %.4f, t %6.2f, p = %.4f)\n', ...
        bPrice, seP, tP, 2*(1 - 0.5*(1 + erf(abs(tP)/sqrt(2)))));
fprintf('    nest parameter   sigma  = %8.4f (se %.4f, t %6.2f, p = %.3g)\n', ...
        sigma, seS, tS, 2*(1 - 0.5*(1 + erf(abs(tS)/sqrt(2)))));
fprintf('    R-squared = %.3f | n = %d | clusters = %d\n', stats.R2, D.n, stats.G);

% robustness check with log price per serving
Xl = [ones(D.n, 1), dummies(D.market), log(D.pps), D.lnsjg];
nl_ = [{'const'}, mnames, {'lnpps', 'lns_jg'}];
[thL, VL] = ols_cluster(Xl, D.y, grp2id(D.prod));
bLog = thL(strcmp(nl_, 'lnpps')); sLog = thL(strcmp(nl_, 'lns_jg'));
fprintf('    robustness (log price): b = %.4f (se %.4f), sigma = %.4f\n\n', ...
        bLog, sqrt(VL(end-1, end-1)), sLog);

% market fixed effects (mean utility level of each flavour market)
D.fe = zeros(D.n, 1);
for m = 1:numel(mkts)
    ic = find(strcmp(names, ['mkt_' mkts{m}]), 1);
    if ~isempty(ic), lvl = theta(1) + theta(ic); else, lvl = theta(1); end
    D.fe(strcmp(D.market, mkts{m})) = lvl;
end

% ============================================ 4. ELASTICITIES & THE MATRIX
fprintf('[3] Own- and cross-price elasticities\n');
own = zeros(D.n, 1);
for m = 1:numel(mkts)
    im = find(strcmp(D.market, mkts{m}));
    [~, o] = elasticity_matrix(bPrice, sigma, D.pps(im), D.nest(im), ...
                               D.share(im), D.sjg(im));
    own(im) = o;
end
D.own = own;
fprintf('    mean own-price elasticity : %7.3f\n', mean(own));
fprintf('      ... protein powders     : %7.3f\n', mean(own(strcmp(D.nest, 'powder'))));
fprintf('      ... protein bars        : %7.3f\n', mean(own(strcmp(D.nest, 'bar'))));
fprintf('      range                   : %7.3f to %7.3f\n\n', min(own), max(own));

for f = 1:numel(cfg.focalMarkets)
    fm = cfg.focalMarkets{f};
    im = find(strcmp(D.market, fm));
    if isempty(im), continue; end
    [E, ~] = elasticity_matrix(bPrice, sigma, D.pps(im), D.nest(im), ...
                               D.share(im), D.sjg(im));
    labels = strcat(strtrim(D.brand(im)), ' - ', strtrim(D.flavor(im)));
    fprintf('    ELASTICITY MATRIX - %s market (%d products, %d brands)\n', ...
            fm, numel(im), numel(unique(D.brand(im))));
    fprintf('    diagonal = own-price elasticity, off-diagonal = cross-price effect\n');
    print_matrix(E, labels, 7);
    write_matrix(fullfile(cfg.outDir, ['elasticity_matrix_' tag(fm) '.csv']), E, labels);
    [B, Mb] = brand_matrix(D.brand(im), D.share(im), E);
    fprintf('    BRAND-LEVEL substitution matrix - %s\n', fm);
    print_matrix(Mb, B, 7);
    write_matrix(fullfile(cfg.outDir, ['brand_matrix_' tag(fm) '.csv']), Mb, B);
    fprintf('\n');
end
write_own_table(fullfile(cfg.outDir, 'own_elasticities.csv'), D, own);

% ================================================================ 5. FIGURES
fprintf('[4] Writing figures to %s\n', cfg.figDir);
fig1_shares(cfg, D, mkts);
fig2_price_share(cfg, D, bPrice, sigma);
fig3_heatmap(cfg, D, bPrice, sigma);
fig4_price_sweep(cfg, D, bPrice, sigma);
fig5_fit(cfg, D, bPrice, sigma);

% ============================================================== 6. VALIDATION
fprintf('\n[5] VALIDATION\n');
resHold  = validate_catalog_holdout(cfg, D, bPrice, sigma);
resPromo = validate_promo_holdout(cfg);

% ================================================================ 7. SUMMARY
results = struct();
results.config = cfg;
results.beta_price = bPrice;  results.se_price = seP;
results.sigma = sigma;        results.se_sigma = seS;
results.beta_logpps = bLog;   results.sigma_logpps = sLog;
results.R2 = stats.R2;        results.n = D.n;
results.own_mean   = mean(own);
results.own_powder = mean(own(strcmp(D.nest, 'powder')));
results.own_bar    = mean(own(strcmp(D.nest, 'bar')));
results.holdout = resHold;
results.promo   = resPromo;
write_summary(fullfile(cfg.outDir, 'estimation_summary.csv'), results);
fprintf('\n[6] Done. Summary written to output/estimation_summary.csv\n');
fprintf('==================================================================\n');
end

% =========================================================================
%  CHOICE MARKETS - map a flavour string onto its flavour profile
% =========================================================================
function m = flavour_market(name)
rules = { ...
  'Birthday Cake',    'birthday'; ...
  'Salted Caramel',   'salted caramel|caramel'; ...
  'Mocha / Coffee',   'mocha|coffee|cappuccino|espresso|cafe|latte'; ...
  'Cookies & Cream',  'cookie|oreo'; ...
  'Chocolate',        'chocolate|cocoa|fudge|brownie|double rich|milk chocolate'; ...
  'Vanilla',          'vanilla|creme|cream$'; ...
  'Strawberry/Berry', 'strawberr|berry|raspberry'; ...
  'Peanut Butter',    'peanut|pb&j'; ...
  'Banana',           'banana'; ...
  'Mint',             '\bmint\b'; ...
  'Unflavored',       'unflavored|unflavoured|flavorless'};
s = lower(strtrim(name));
m = 'Other';
for k = 1:size(rules, 1)
    if ~isempty(regexp(s, rules{k, 2}, 'once')), m = rules{k, 1}; return; end
end
end

% =========================================================================
%  NESTED-LOGIT ELASTICITY MATRIX for one market
%    price enters utility linearly, so d delta_m / d p_m = b; the column
%    factor pps(m) converts d ln s_j / d p_m into the standard elasticity
%    d ln s_j / d ln p_m:
%    own          : b*pps(m)*( (1 - sigma*s_jg)/(1-sigma) - s_jg*s_g )
%    cross, same  : -b*pps(m)* s_mg*( sigma/(1-sigma) + s_g )
%    cross, other : -b*pps(m)* s_m
% =========================================================================
function [E, own] = elasticity_matrix(b, sigma, pps, nests, share, sjg)
n = numel(pps);
E = zeros(n, n);
gs = zeros(n, 1);
for j = 1:n, gs(j) = sum(share(strcmp(nests, nests{j}))); end
oms = 1 - sigma;
for j = 1:n
    for m = 1:n
        if j == m
            d = (1 - sigma*sjg(j))/oms - sjg(j)*gs(j);
        elseif strcmp(nests{j}, nests{m})
            d = -sjg(m)*(sigma/oms + gs(m));
        else
            d = -share(m);
        end
        E(j, m) = b * pps(m) * d;
    end
end
own = diag(E);
end

% =========================================================================
%  brand-level elasticity matrix: elasticity of brand a's aggregate share
%  w.r.t. a uniform +/-1% price change in every SKU of brand c
%    Mb(a,c) = sum_{i in a} (s_i/S_a) * sum_{j in c} E(i,j)
% =========================================================================
function [B, Mb] = brand_matrix(brands, share, E)
B = unique(brands);
nb = numel(B);
Mb = zeros(nb, nb);
for a = 1:nb
    ia = find(strcmp(brands, B{a}));
    wa = share(ia) / sum(share(ia));
    for c = 1:nb
        ic = find(strcmp(brands, B{c}));
        row = sum(E(ia, ic), 2);   % effect on SKU i of brand-c price change
        Mb(a, c) = sum(wa .* row);
    end
end
end

% =========================================================================
%  closed-form nested-logit share mapping for one market
% =========================================================================
function [s, s0] = nl_shares(delta, nests, sigma)
n = numel(delta);
delta = delta(:);
oms = 1 - sigma;
u = unique(nests);
v = delta / oms;
lnD = zeros(numel(u), 1);
for k = 1:numel(u)
    vk = v(strcmp(nests, u{k}));
    mx = max(vk);
    lnD(k) = mx + log(sum(exp(vk - mx)));
end
terms = [0; oms * lnD];
mx = max(terms);
lnDen = mx + log(sum(exp(terms - mx)));
s = zeros(n, 1);
for k = 1:numel(u)
    ik = find(strcmp(nests, u{k}));
    for q = 1:numel(ik)
        s(ik(q)) = exp(v(ik(q)) - sigma*lnD(k) - lnDen);
    end
end
s0 = exp(-lnDen);
end

% =========================================================================
%  OLS with cluster-robust standard errors
% =========================================================================
function [theta, V, stats] = ols_cluster(X, y, groups)
[N, K] = size(X);
XtX = X' * X;
theta = XtX \ (X' * y(:));
u = y(:) - X*theta;
S = zeros(K, K);
ug = unique(groups);
G = numel(ug);
for g = 1:G
    m = groups == ug(g);
    a = X(m, :)' * u(m);
    S = S + a * a';
end
bread = XtX \ eye(K);
V = bread * S * bread * (G/(G-1)) * ((N-1)/(N-K));
stats = struct('R2', 1 - sum(u.^2)/sum((y(:) - mean(y)).^2), ...
               'N', N, 'K', K, 'G', G);
end

% =========================================================================
%  VALIDATION A - held-out SKUs of the catalogue
% =========================================================================
function out = validate_catalog_holdout(cfg, D, bPrice, sigma)
fprintf('    A. Catalogue hold-out: %.0f%% of SKUs held out at random\n', ...
        100*cfg.holdoutFraction);
% deterministic, reproducible split (no random generator needed)
key = zeros(D.n, 1);
for k = 1:D.n
    v = double(char(D.prod{k}));
    key(k) = mod(sum(v .* (1:numel(v))), 5);
end
hold_ = key < round(5*cfg.holdoutFraction);
train = ~hold_;
mk = unique(D.market);
if sum(train) < 10 || sum(hold_) < 5
    out = struct('n', 0, 'rmspe', NaN, 'mape', NaN, 'corr', NaN, ...
                 'naive_rmspe', NaN, 'beta_train', NaN, 'sigma_train', NaN);
    return;
end

% --- estimate on the training SKUs, same specification as the baseline
mkTr = unique(D.market(train));
Xt = [ones(sum(train), 1), dummies(D.market(train)), D.pps(train), D.lnsjg(train)];
th = Xt \ D.y(train);
resid = D.y(train) - Xt*th;
trM = D.market(train);
fe = zeros(numel(mk), 1);
for k = 1:numel(mk)
    kk = find(strcmp(mkTr, mk{k}), 1);   % column of this market dummy
    lvl = th(1);
    if ~isempty(kk) && kk > 1, lvl = th(1) + th(kk); end
    m = strcmp(trM, mk{k});
    if any(m), fe(k) = lvl + mean(resid(m)); else, fe(k) = lvl; end
end

holdIdx = find(hold_);
dh = zeros(numel(holdIdx), 1);
for k = 1:numel(holdIdx)
    j = holdIdx(k);
    kk = find(strcmp(mk, D.market{j}), 1);
    if isempty(kk), kk = 1; end
    dh(k) = fe(kk) + th(end-1)*D.pps(j);
end
mkH = D.market(holdIdx); nestsH = D.nest(holdIdx);
ums = unique(mkH);
sPred = zeros(numel(holdIdx), 1); sAct = zeros(numel(holdIdx), 1);
sNai  = zeros(numel(holdIdx), 1);
for k = 1:numel(ums)        g = find(strcmp(mkH, ums{k}));
        % strict out-of-sample: share mapping uses sigma from training SKUs only
        [sp, ~] = nl_shares(dh(g), nestsH(g), th(end));
    sPred(g) = sp / sum(sp);
    sa = D.share(holdIdx(g));  sAct(g) = sa / sum(sa);
    trm = D.share(train & strcmp(D.market, ums{k}));
    if isempty(trm), sn = ones(numel(g), 1); else, sn = ones(numel(g), 1) * mean(trm); end
    sNai(g) = sn / sum(sn);
end
e  = sPred - sAct;
en = sNai  - sAct;
cc = corrcoef(sPred, sAct);
cn = corrcoef(sNai, sAct);
out = struct('n', numel(sPred), ...
             'rmspe', sqrt(mean(e.^2)), ...
             'mape', mean(abs(e) ./ max(sAct, 1e-12)), ...
             'corr', cc(1, 2), ...
             'naive_corr', cn(1, 2), ...
             'naive_rmspe', sqrt(mean(en.^2)), ...
             'beta_train', th(end-1), 'sigma_train', th(end));
fprintf('      RMSPE = %.4f (benchmark %.4f) | MAPE = %.1f%% | corr = %.3f (benchmark %.3f)\n', ...
        out.rmspe, out.naive_rmspe, 100*out.mape, out.corr, out.naive_corr);
fprintf('      n held out = %d | training-sample parameters: b = %.4f, sigma = %.4f\n', ...
        out.n, out.beta_train, out.sigma_train);

fig = figure('Visible', 'off', 'Position', [100 100 900 600]);
plot(sAct, sPred, 'o', 'MarkerSize', 5); hold on;
lims = [0, max([sAct; sPred]) * 1.05];
plot(lims, lims, 'r-'); grid on;
xlabel('Actual share (held-out SKUs)', 'FontSize', 11);
ylabel('Predicted share', 'FontSize', 11);
title(sprintf('Validation A: held-out SKUs (RMSPE %.4f, corr %.2f)', ...
      out.rmspe, out.corr), 'FontSize', 12);
print(fig, fullfile(cfg.figDir, 'fig6a_holdout_skus.png'), '-dpng', '-r300');
close(fig);
end

% =========================================================================
%  VALIDATION B - held-out PROMOTIONAL weeks of the weekly sales panel
% =========================================================================
function out = validate_promo_holdout(cfg)
fprintf('    B. Weekly panel: hold out promotional weeks (mean discount >= %.0f%%)\n', ...
        100*cfg.promoThreshold);
P = read_csv(cfg.panelFile);
h = P(1, :); c = @(n) find(strcmp(h, n), 1);
nObs = size(P, 1) - 1;
p.date     = P(2:end, c('Date'));
p.prod     = P(2:end, c('Product Name'));
p.units    = str2double(P(2:end, c('Units Sold')));
p.price    = str2double(P(2:end, c('Price')));
p.discount = str2double(P(2:end, c('Discount')));

powderPanel = {'Whey Protein', 'BCAA', 'Creatine', 'Pre-Workout', ...
               'Electrolyte Powder', 'Collagen Peptides'};
p.nest = repmat({'solidose'}, nObs, 1);
for k = 1:nObs
    if any(strcmp(p.prod{k}, powderPanel)), p.nest{k} = 'powder'; end
end
wk = unique(p.date);  nw = numel(wk);
meanDisc = zeros(nw, 1);
promoWeek = zeros(nObs, 1);
isPromo = false(nw, 1);
for k = 1:nw
    m = strcmp(p.date, wk{k});
    meanDisc(k) = mean(p.discount(m));
    isPromo(k) = meanDisc(k) >= cfg.promoThreshold;
    promoWeek(m) = isPromo(k);
end
fprintf('      %d weeks total, %d held out as promotional (%.0f%%)\n', ...
        nw, sum(isPromo), 100*mean(isPromo));

s0a = 0.5;                       % assumed outside-good share (see paper)
p.share = zeros(nObs, 1);
for k = 1:nw
    m = strcmp(p.date, wk{k});
    p.share(m) = (1 - s0a) * p.units(m) / sum(p.units(m));
end
p.sg = zeros(nObs, 1); p.sjg = zeros(nObs, 1);
for k = 1:nObs
    m = strcmp(p.date, p.date{k}) & strcmp(p.nest, p.nest{k});
    p.sg(k)  = sum(p.share(m));
    p.sjg(k) = p.share(k) / p.sg(k);
end
p.y = log(p.share / s0a);
p.lnsjg = log(p.sjg);
p.lnp = log(p.price);
tr = promoWeek == 0;

% --- (i) pooled NL with week fixed effects on regular weeks (coefficients)
Xw = [ones(sum(tr), 1), dummies(p.date(tr)), p.lnp(tr), p.lnsjg(tr)];
[thw, Vw] = ols_cluster(Xw, p.y(tr), grp2id(p.prod(tr)));
bP = thw(end - 1); sS = thw(end);
seB = sqrt(Vw(end - 1, end - 1));
fprintf('      price coefficient = %.4f (se %.4f, p = %.2f) on regular weeks\n', ...
        bP, seB, 2*(1 - 0.5*(1 + erf(abs(bP/seB)/sqrt(2)))));
fprintf('      sigma = %.4f (se %.4f)\n', sS, sqrt(Vw(end, end)));
cc0 = corrcoef(p.price, p.units);
fprintf('      raw correlation of price and units in the panel: %.3f\n', cc0(1, 2));

% --- (ii) prediction model with product fixed effects
Xp = [ones(sum(tr), 1), dummies(p.prod(tr)), p.lnp(tr)];
thp = Xp \ p.y(tr);
pbase = unique(p.prod(tr));
basePrice = zeros(numel(pbase), 1);
baseShare = zeros(numel(pbase), 1);
for k = 1:numel(pbase)
    m = tr & strcmp(p.prod, pbase{k});
    basePrice(k) = median(p.price(m));
    baseShare(k) = mean(p.share(m));
end

% --- predict every product in every held-out promotional week
promIdx = find(promoWeek == 1);
sPred = zeros(numel(promIdx), 1); sAct = zeros(numel(promIdx), 1);
sBase = zeros(numel(promIdx), 1); sNai = zeros(numel(promIdx), 1);
promWeeks = wk(isPromo);
for k = 1:numel(promWeeks)
    g = find(strcmp(p.date, promWeeks{k}));
    pos = zeros(numel(g), 1);
    for q = 1:numel(g), pos(q) = find(promIdx == g(q), 1); end
    % predicted shares at the observed (promotional) prices
    dh = thp(end) * p.lnp(g);
    lb = zeros(numel(g), 1);
    nb = zeros(numel(g), 1);
    for q = 1:numel(g)
        cc = find(strcmp(pbase, p.prod{g(q)}), 1);
        if isempty(cc)
            dh(q) = dh(q) + thp(1);
            lb(q) = thp(1) + thp(end) * log(p.price(g(q)));
            nb(q) = mean(p.share(tr));
        else
            dh(q) = dh(q) + thp(1);
            if cc > 1, dh(q) = dh(q) + thp(cc); end
            lb(q) = thp(1) + thp(end) * log(basePrice(cc));
            if cc > 1, lb(q) = lb(q) + thp(cc); end
            nb(q) = baseShare(cc);
        end
    end
    [sp, ~] = nl_shares(dh, p.nest(g), sS);
    [sb, ~] = nl_shares(lb, p.nest(g), sS);
    sa = p.share(g);
    sPred(pos) = sp / sum(sp);
    sAct(pos)  = sa  / sum(sa);
    sBase(pos) = sb  / sum(sb);
    sNai(pos)  = nb  / sum(nb);
end
e  = sPred - sAct;
en = sNai  - sAct;
cc = corrcoef(sPred, sAct);
liftA = sAct ./ sNai - 1;      % actual lift vs regular-period benchmark
liftP = sPred ./ sBase - 1;    % model-predicted promotional lift
out = struct('n', numel(sPred), ...
             'rmspe', sqrt(mean(e.^2)), ...
             'mape', mean(abs(e) ./ max(sAct, 1e-12)), ...
             'corr', cc(1, 2), ...
             'naive_rmspe', sqrt(mean(en.^2)), ...
             'price_coef', bP, 'sigma', sS, ...
             'corr_price_units', cc0(1, 2), ...
             'lift_actual', mean(liftA), 'lift_pred', mean(liftP), ...
             'lift_rmse', sqrt(mean((liftA - liftP).^2)));
fprintf('      hold-out RMSPE = %.4f (benchmark %.4f), MAPE = %.1f%%, corr = %.3f\n', ...
        out.rmspe, out.naive_rmspe, 100*out.mape, out.corr);
fprintf('      promotional lift: actual %+.3f, predicted %+.3f (RMSE %.3f)\n', ...
        out.lift_actual, out.lift_pred, out.lift_rmse);

fig = figure('Visible', 'off', 'Position', [100 100 900 600]);
plot(sAct, sPred, 'o', 'MarkerSize', 5); hold on;
lims = [0, max([sAct; sPred]) * 1.05];
plot(lims, lims, 'r-'); grid on;
xlabel('Actual share (promotional weeks)', 'FontSize', 11);
ylabel('Predicted share', 'FontSize', 11);
title(sprintf('Validation B: held-out promotional weeks (RMSPE %.4f)', out.rmspe), ...
      'FontSize', 12);
print(fig, fullfile(cfg.figDir, 'fig6b_holdout_promotions.png'), '-dpng', '-r300');
close(fig);

fig = figure('Visible', 'off', 'Position', [100 100 900 600]);
edges = linspace(min([liftA; liftP]), max([liftA; liftP]), 30);
bar_hist(liftA, edges); hold on;
bar_hist(liftP, edges);
legend({'actual promotional lift', 'model-predicted lift'}, 'Location', 'northeast');
xlabel('Share lift in promotional weeks', 'FontSize', 11);
ylabel('Frequency', 'FontSize', 11);
title('Validation B: promotional lift, actual vs predicted', 'FontSize', 12);
grid on;
print(fig, fullfile(cfg.figDir, 'fig7_promo_lift.png'), '-dpng', '-r300');
close(fig);
end

% =========================================================================
%  FIGURES
% =========================================================================
function fig1_shares(cfg, D, mkts)
fig = figure('Visible', 'off', 'Position', [100 100 1000 560]);
pw = zeros(numel(mkts), 1); br = pw;
for k = 1:numel(mkts)
    m = strcmp(D.market, mkts{k});
    pw(k) = sum(D.share(m & strcmp(D.nest, 'powder')));
    br(k) = sum(D.share(m & strcmp(D.nest, 'bar')));
end
[~, order] = sort(pw + br, 'descend');
bar([pw(order), br(order)], 'stacked'); grid on;
set(gca, 'XTick', 1:numel(mkts), 'XTickLabel', mkts(order));
legend({'protein powders', 'protein bars'}, 'Location', 'northeast');
ylabel('Share of catalogue review mass', 'FontSize', 11);
title('Market shares by flavour profile and nest', 'FontSize', 13);
print(fig, fullfile(cfg.figDir, 'fig1_market_shares.png'), '-dpng', '-r300');
close(fig);
end

function fig2_price_share(cfg, D, b, sigma)
fig = figure('Visible', 'off', 'Position', [100 100 1000 560]);
ip = strcmp(D.nest, 'powder'); ib = strcmp(D.nest, 'bar');
scatter(D.pps(ip), D.share(ip), 36, [0 0.45 0.74], 'filled'); hold on;
scatter(D.pps(ib), D.share(ib), 72, [0.85 0.33 0.10], 'd', 'filled');
set(gca, 'YScale', 'log'); grid on;
xlabel('Price per serving (USD)', 'FontSize', 11);
ylabel('Market share (log scale)', 'FontSize', 11);
legend({'protein powders', 'protein bars'}, 'Location', 'northeast');
title(sprintf('Price per serving vs market share  (b = %.3f, sigma = %.3f)', b, sigma), ...
      'FontSize', 13);
print(fig, fullfile(cfg.figDir, 'fig2_price_share.png'), '-dpng', '-r300');
close(fig);
end

function fig3_heatmap(cfg, D, b, sigma)
for f = 1:numel(cfg.focalMarkets)
    fm = cfg.focalMarkets{f};
    im = find(strcmp(D.market, fm));
    if isempty(im), continue; end
    [E, ~] = elasticity_matrix(b, sigma, D.pps(im), D.nest(im), ...
                               D.share(im), D.sjg(im));
    labels = strcat(strtrim(D.brand(im)), ' - ', strtrim(D.flavor(im)));
    fig = figure('Visible', 'off', 'Position', [100 100 1100 880]);
    imagesc(E); colorbar; colormap(redblue());
    lim = max(abs(E(:))); if lim == 0, lim = 1; end
    caxis([-lim lim]);
    set(gca, 'XTick', 1:numel(im), 'XTickLabel', labels, 'XTickLabelRotation', 90, ...
             'YTick', 1:numel(im), 'YTickLabel', labels, 'FontSize', 7);
    title(sprintf('Own- and cross-price elasticities: %s market', fm), 'FontSize', 13);
    xlabel('price change of product m  \rightarrow', 'FontSize', 11);
    ylabel('response of product j', 'FontSize', 11);
    print(fig, fullfile(cfg.figDir, ['fig3_elasticity_heatmap_' tag(fm) '.png']), ...
          '-dpng', '-r300');
    close(fig);
end
end

function fig4_price_sweep(cfg, D, b, sigma)
% substitution patterns across price points: sweep the price of the market
% leader and trace the implied shares of the leader and its closest rivals
fm = 'Mocha / Coffee';
im = find(strcmp(D.market, fm));
if isempty(im), fm = 'Birthday Cake'; im = find(strcmp(D.market, fm)); end
if isempty(im), return; end
[~, ord] = sort(D.share(im), 'descend');
im = im(ord(1:min(6, numel(ord))));
focal = 1;
sweep = linspace(0.6, 1.4, 41);
sh = zeros(numel(sweep), numel(im));
for k = 1:numel(sweep)
    delta = D.fe(im) + b * D.pps(im);
    delta(focal) = D.fe(im(focal)) + b * D.pps(im(focal)) * sweep(k);
    [s, ~] = nl_shares(delta, D.nest(im), sigma);
    sh(k, :) = s(:)';
end
fig = figure('Visible', 'off', 'Position', [100 100 1000 600]);
labels = cell(numel(im), 1);
for j = 1:numel(im)
    labels{j} = [strtrim(D.brand{im(j)}), ' - ', strtrim(D.flavor{im(j)})];
    plot(100*(sweep - 1), 100*sh(:, j), 'LineWidth', 1.6); hold on;
end
grid on;
xlabel(sprintf('price change of %s (%%)', labels{focal}), 'FontSize', 11);
ylabel('market share (%%)', 'FontSize', 11);
legend(labels, 'Location', 'northeast', 'FontSize', 8);
title(sprintf('Substitution across price points: %s market', fm), 'FontSize', 13);
print(fig, fullfile(cfg.figDir, 'fig4_price_sweep.png'), '-dpng', '-r300');
close(fig);
end

function fig5_fit(cfg, D, b, sigma)
delta = D.fe + b * D.pps;
mk = unique(D.market);
sHat = zeros(D.n, 1); sAct = zeros(D.n, 1);
for k = 1:numel(mk)
    g = find(strcmp(D.market, mk{k}));
    [s, ~] = nl_shares(delta(g), D.nest(g), sigma);
    sHat(g) = s / sum(s);
    sa = D.share(g);
    sAct(g) = sa / sum(sa);
end
fig = figure('Visible', 'off', 'Position', [100 100 900 600]);
plot(sAct, sHat, 'o', 'MarkerSize', 5); hold on;
lims = [0, max([sAct; sHat]) * 1.05];
plot(lims, lims, 'r-'); grid on;
xlabel('Actual share (normalised within market)', 'FontSize', 11);
ylabel('Fitted share', 'FontSize', 11);
title('Nested-logit fit: actual vs fitted market shares', 'FontSize', 12);
print(fig, fullfile(cfg.figDir, 'fig5_fit.png'), '-dpng', '-r300');
close(fig);
end

% =========================================================================
%  SMALL UTILITIES (toolbox-free, MATLAB & Octave compatible)
% =========================================================================
function C = read_csv(fname)
fid = fopen(fname, 'r');
if fid < 0, error('Cannot open %s', fname); end
C = {};
tline = fgetl(fid);
while ischar(tline)
    parts = strsplit(tline, ',', 'CollapseDelimiters', false);
    C = [C; parts]; %#ok<AGROW>
    tline = fgetl(fid);
end
fclose(fid);
end

function x = to_number(strCells)
% strip currency symbols, thousands separators and non-ASCII blanks
n = numel(strCells);
x = zeros(n, 1);
for k = 1:n
    s = regexprep(strCells{k}, '[^0-9.\-+]', '');
    x(k) = str2double(s);
end
end

function D = dummies(names)
u = unique(names);
D = zeros(numel(names), numel(u) - 1);
for k = 2:numel(u), D(:, k-1) = double(strcmp(names, u{k})); end
end

function nm = dummy_names(names)
u = unique(names);
nm = cell(1, numel(u) - 1);
for k = 2:numel(u), nm{k-1} = ['mkt_' u{k}]; end
end

function id = grp2id(g)
u = unique(g);
id = zeros(numel(g), 1);
for k = 1:numel(g), id(k) = find(strcmp(u, g{k}), 1); end
end

function t = tag(s)
t = regexprep(s, '[^a-zA-Z0-9]+', '_');
end

function print_matrix(E, labels, prec)
n = size(E, 1);
w = 11;
fprintf('    %*s', w, '');
for j = 1:n, fprintf('%*s', w, abbrev(labels{j}, w)); end
fprintf('\n');
for i = 1:n
    fprintf('    %*s', w, abbrev(labels{i}, w));
    for j = 1:n, fprintf('%*.*f', w, prec - 1, E(i, j)); end
    fprintf('\n');
end
end

function s = abbrev(s, maxlen)
s = regexprep(s, '\s+', ' ');
if numel(s) > maxlen, s = [s(1:maxlen-1) '~']; end
end

function write_matrix(fname, E, labels)
fid = fopen(fname, 'w');
fprintf(fid, 'product');
for j = 1:numel(labels), fprintf(fid, ',%s', labels{j}); end
fprintf(fid, '\n');
for i = 1:numel(labels)
    fprintf(fid, '%s', labels{i});
    for j = 1:numel(labels), fprintf(fid, ',%.6f', E(i, j)); end
    fprintf(fid, '\n');
end
fclose(fid);
end

function write_own_table(fname, D, own)
fid = fopen(fname, 'w');
fprintf(fid, 'market,brand,product,flavor,nest,price_per_serving,share,own_elasticity\n');
for k = 1:D.n
    fprintf(fid, '%s,%s,%s,%s,%s,%.4f,%.8f,%.4f\n', D.market{k}, D.brand{k}, ...
            D.prod{k}, D.flavor{k}, D.nest{k}, D.pps(k), D.share(k), own(k));
end
fclose(fid);
end

function write_summary(fname, r)
fid = fopen(fname, 'w');
fprintf(fid, 'item,value\n');
fprintf(fid, 'n_observations,%d\n', r.n);
fprintf(fid, 'price_coefficient,%.6f\n', r.beta_price);
fprintf(fid, 'price_std_error,%.6f\n', r.se_price);
fprintf(fid, 'price_coefficient_logspec,%.6f\n', r.beta_logpps);
fprintf(fid, 'sigma,%.6f\n', r.sigma);
fprintf(fid, 'sigma_std_error,%.6f\n', r.se_sigma);
fprintf(fid, 'sigma_logspec,%.6f\n', r.sigma_logpps);
fprintf(fid, 'r_squared,%.6f\n', r.R2);
fprintf(fid, 'mean_own_elasticity,%.6f\n', r.own_mean);
fprintf(fid, 'mean_own_elasticity_powders,%.6f\n', r.own_powder);
fprintf(fid, 'mean_own_elasticity_bars,%.6f\n', r.own_bar);
fprintf(fid, 'holdout_sku_rmspe,%.6f\n', r.holdout.rmspe);
fprintf(fid, 'holdout_sku_benchmark_rmspe,%.6f\n', r.holdout.naive_rmspe);
fprintf(fid, 'holdout_sku_corr,%.6f\n', r.holdout.corr);
fprintf(fid, 'holdout_sku_benchmark_corr,%.6f\n', r.holdout.naive_corr);
fprintf(fid, 'holdout_sku_n,%d\n', r.holdout.n);
fprintf(fid, 'holdout_sku_beta,%.6f\n', r.holdout.beta_train);
fprintf(fid, 'holdout_sku_sigma,%.6f\n', r.holdout.sigma_train);
fprintf(fid, 'promo_holdout_rmspe,%.6f\n', r.promo.rmspe);
fprintf(fid, 'promo_holdout_benchmark_rmspe,%.6f\n', r.promo.naive_rmspe);
fprintf(fid, 'promo_holdout_corr,%.6f\n', r.promo.corr);
fprintf(fid, 'panel_price_coefficient,%.6f\n', r.promo.price_coef);
fprintf(fid, 'panel_sigma,%.6f\n', r.promo.sigma);
fprintf(fid, 'panel_corr_price_units,%.6f\n', r.promo.corr_price_units);
fprintf(fid, 'promo_lift_actual,%.6f\n', r.promo.lift_actual);
fprintf(fid, 'promo_lift_predicted,%.6f\n', r.promo.lift_pred);
fprintf(fid, 'promo_lift_rmse,%.6f\n', r.promo.lift_rmse);
fclose(fid);
end

function bar_hist(x, edges)
n = numel(edges) - 1;
cnt = zeros(n, 1);
for k = 1:n
    cnt(k) = sum(x >= edges(k) & x < edges(k+1));
end
centers = (edges(1:end-1) + edges(2:end)) / 2;
bar(centers, cnt);
end

function cmap = redblue()
n = 64;
r = [linspace(0, 1, n/2), ones(1, n/2)]';
g = [linspace(0, 1, n/2), linspace(1, 0, n/2)]';
b = [ones(1, n/2), linspace(1, 0, n/2)]';
cmap = [r, g, b];
end
