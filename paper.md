Price Elasticity and Format-Flavor Substitution in Sports Nutrition: A Nested-Logit Analysis of Differentiated Consumer Goods

**Aastha Boora — BFIN515: Economic Analysis — October 2026**

GitHub Repository: https://github.com/AasthaBoora/Price-Elasticity-Substitution-in-Sports-Nutrition-A-Nested-Logit-Analysis

**Abstract.** I estimate own-price elasticities and cross-brand substitution across two differentiated formats — whey powders and protein bars — using a Berry (1994) nested logit estimated by Berry inversion with cluster-robust OLS on 861 product–flavour SKUs (157 products, 11 flavour markets), validated against a 274-week promotional panel. The price coefficient is $-0.3633$ ($t = -4.11$) and $\sigma = 0.8311$ ($t = 17.54$), with $R^2 = 0.806$. Every SKU is elastic (mean $-3.45$): bars average $-5.22$ versus powders $-3.09$, and substitution is almost entirely within format. A three-tier validation ladder shows strong internal fit, hold-out performance below a naive benchmark, and no ability to predict promotional lift from a panel with a zero price coefficient. I derive asymmetric format-pricing, flavour-portfolio, and promotional-calendar recommendations, and delineate where static choice models may and may not guide pricing.

**JEL classification:** L81, M31, D12. **Keywords:** nested logit, price elasticity, demand estimation, sports nutrition, promotion evaluation.

## 1. Business Question

How should category directors and brand pricing managers set list prices and promotional depth across supplement formats and flavour profiles? A brand such as Optimum Nutrition, Quest or BSN sells both a *format* — powder tubs versus protein bars — and a *flavour profile* within each format. Each cell of that grid has its own demand sensitivity, and pricing one cell moves the others through substitution.

Mispricing cuts both ways. A deeply funded bar promotion may only transfer units from the brand's own powder tubs — cannibalizing margin instead of growing the category. Uniform pricing forfeits genuinely inelastic niche flavours: the loyalty that lets a Birthday Cake powder carry a premium. Both errors share one root — pricing without estimates of own-price elasticity and of the cross-elasticity structure that says who takes whose volume.

Brand managers need these estimates to set list prices, design price-pack architecture, and choose which SKUs fund promotions. Retail category buyers need them to allocate shelf space, sequence promotional calendars, and negotiate vendor funding: knowing that a temporary price cut on a bar moves volume almost exclusively between bars — and essentially never into powders — changes how promotions are bundled. I deliver both sets of estimates, with an honest account of where they can and cannot be trusted.

## 2. Data

Demand models can be supported at Level 1 (category series), Level 2 (product-level shares, prices, characteristics, ownership), or transaction microdata (Dong, 2026). My data sit squarely at Level 2 — enough to estimate differentiated demand, but without cost shifters or time-varying consumer states, a boundary that matters below.

**Estimation sample.** I use a public retail catalogue (Kaggle, *Sports nutrition supplements*; 2,719 SKUs, 108 brands). I retain the five protein-powder categories and protein bars with at least one review; assign each SKU to a bar or powder nest by category and product form, reassigning 75 misclassified rows; map flavour strings to flavour markets; and split each parent product's review mass across its flavour SKUs as the share proxy. The result is $N = 861$ SKUs across 157 products and 11 flavour markets (715 powder, 146 bar), with outside-good shares of $0.884$–$0.997$. Price per serving averages $\$1.64$ — $\$1.45$ for powders versus $\$2.58$ for bars, a 77% format premium.

**Validation panel.** A second Kaggle panel (*Supplement Sales Data*, Apache 2.0) covers 16 products over 274 weeks (4,384 product-weeks, 2020–2025) with units, prices, and discounts. Weeks averaging at least 13% discount — 1,680 promotional product-weeks — form the held-out regime in Section 5. Neither dataset contains Myprotein or Barebells (verified by search), so those brands enter only as market context.

**Limitations.** Review mass proxies sales; the catalogue is cross-sectional, projecting all time variation onto cross-sectional variation; and zero-review products are dropped, conditioning on observed demand. These boundaries motivate the fixed-effects design and bound Sections 5–6. (Sample exhibits: E1, E2, ET1.)

## 3. Method

I estimate a Berry (1994) Level-2 aggregate nested logit by two-step OLS with standard errors clustered on the 157 parent products. Consumers choose among products nested in powders or bars versus an outside good, with conditional within-nest shares following the nested-logit GEV distribution. Inverting the share system (Berry, 1994) gives my estimating equation:

$$
\ln(s_{jt} / s_{0t}) = \mathbf{x}_{jt}\boldsymbol{\beta} - \alpha \cdot \text{PPS}_{jt} + \sigma \ln(s_{j\vert g,t}) + \xi_{jt},
$$

which I estimate by constructing both share terms from observed shares and running OLS with flavour-profile fixed effects.

The nesting parameter $\sigma \in [0,1)$ governs within-nest substitution. As $\sigma \to 0$ the model collapses to multinomial logit, whose IIA restriction implies — counterfactually here — that powder loyalists and bar shoppers divert identically. As $\sigma \to 1$, goods within a nest become near-perfect substitutes and substitution mass concentrates inside the nest, eliminating the cross-format substitution that IIA would impose.

Closed-form own, intra-nest, and inter-nest elasticities follow by differentiating the share system; I verified them against five-point numerical differentiation (maximum deviation below $10^{-9}$) and aggregate them to brand level as share-weighted sums.

Two devices address unobserved quality $\xi_{jt}$. Flavour-profile fixed effects absorb market-level taste shocks, so $\alpha$ is identified from price dispersion within a flavour market, net of nest membership. Product-level clustering allows arbitrary correlation across a parent product's variants. Without cost instruments, $\alpha$ may absorb correlated unobserved quality, so I treat all magnitudes as upper bounds.

## 4. Results

Table 1 reports the main estimates. The price-per-serving coefficient is $-0.3633$ (SE $= 0.0884$, $t = -4.11$, $p < 0.001$): a one-dollar increase in price per serving cuts $\ln(s_{jt}/s_{0t})$ by 0.36 within a flavour market. The nesting parameter is $\sigma = 0.8311$ (SE $= 0.0474$, $t = 17.54$) — decisively above zero, so IIA is rejected and substitution is overwhelmingly within format. Fit is strong for cross-sectional demand work ($R^2 = 0.806$). A log-price robustness spec gives $\beta = -0.7407$ ($p < 0.001$) and $\sigma = 0.8223$, so the conclusions do not depend on price functional form.

**Table 1 — Nested-logit estimates (Berry inversion, two-step OLS, clustered SEs)**

| | Main spec (PPS, \$) | Log-price spec |
|---|---|---|
| Price coefficient | $-0.3633$ (0.0884) | $-0.7407$ (0.1645) |
| $t$-statistic | $-4.11$ | $-4.50$ |
| Nesting parameter $\sigma$ | 0.8311 (0.0474) | 0.8223 |
| $R^2$ / Adjusted | 0.806 / 0.803 | — |
| Observations | 861 | 861 |

Table 2 converts the estimates into elasticities over all 861 SKUs. Mean own-price elasticity is $-3.453$ (median $-2.940$; range $-11.25$ to $-1.06$); every SKU is elastic and 81% exceed $\lvert E \rvert = 2$. The decisive split: **powders average $-3.091$ versus bars $-5.221$** — bar demand is far more price-sensitive, consistent with smaller purchase size, lower switching costs, and a denser competitive set.

**Table 2 — Format elasticity summary (own-price elasticities, $N = 861$)**

| Statistic | All SKUs | Protein powders | Protein bars |
|---|---|---|---|
| Mean | $-3.453$ | $-3.091$ | $-5.221$ |
| Median | $-2.940$ | — | — |
| Range | $-11.25$ / $-1.06$ | — | — |
| Share with $\lvert E \rvert > 1$ | 100% | 100% | 100% |

The three focal markets sharpen the reading (Exhibit ET2). Mocha/Coffee (26 SKUs) averages $-3.19$ with intra-nest cross-elasticity of $0.099$ against inter-nest of $0.00013$ — roughly 760 times more intense within format. Salted Caramel (23 SKUs) is the most contested at $-4.10$, with intra-nest $0.317$ and brand-level cross-effects up to $3.54$; Birthday Cake averages $-3.24$. Substitution is asymmetric: Dymatize's share responds to an Optimum Nutrition price change with elasticity $0.913$, while Optimum responds to Dymatize with only $0.092$ — volume flows toward the dominant incumbent (matrices: Exhibits ET3, E3a–E3c).

A counterfactual sweep of the Mocha leader's price across $\pm 40\%$ (Exhibit E4) shows the same structure dynamically: the leader's own curve is steep ($-1.96$) and diversion flows almost entirely to other powders while rival bar lines stay flat. In-sample fit (Exhibit E5) clusters tightly on the 45-degree line for large markets, with dispersion concentrated in small markets — foreshadowing Section 5.

## 5. Validation

I evaluate the model on Dong's (2026) multi-tiered validation ladder, which separates internal fit from increasingly demanding out-of-sample tests.

**Level A — internal fit.** $R^2 = 0.806$ and in-sample share RMSPE $= 0.032$ confirm that the inversion, nests, and fixed effects are correctly specified — coherence, not predictive validity, since the same moments both estimate and evaluate the parameters.

**Level B — same-setting hold-out.** A deterministic hash assigns 210 of 861 SKUs to hold-out; re-estimated on the remaining 80% ($\beta = -0.2826$, $\sigma = 0.8560$), I predict held-out shares. Against a naive market-mean benchmark, the nested model loses: RMSPE $= 0.0871$ versus $0.0719$, correlation $0.587$ versus $0.686$. Plain logit does better still ($0.065$ / $0.752$): with $\hat{\sigma} \approx 0.86$, unequal within-nest shares over-share the wrong product in small, partially observed markets. The nested structure that earns the high $R^2$ overfits (Exhibit E6a).

**Levels C/D — promotional regime.** Withholding the 105 weeks averaging at least 13% discount (1,680 product-weeks) and estimating on regular weeks, hold-out RMSPE is $0.0053$ (MAPE $6.77\%$), but correlation is $-0.0039$ against a naive benchmark of $0.0051$: levels are captured, product-level movement is not (Exhibit E6b). The sharpest test is promotional lift — actual mean lift is $+0.000032$ while predicted lift is $-5\times10^{-9}$ (RMSE $= 0.0812$), four orders of magnitude too small at the mean. The cause is identifiable: the panel's own price coefficient is $-0.0004$ ($p = 0.87$) and the raw price–units correlation is $0.014$. The panel contains no detectable price response from which any choice model could learn (Exhibit E7).

**Table 3 — Validation ladder summary (Dong, 2026)**

| Tier | Exercise | Model | Benchmark | Verdict |
|---|---|---|---|---|
| A — Internal fit | In-sample; $N = 861$ | $R^2 = 0.806$; RMSPE $= 0.032$ | — | Coherence, not prediction |
| B — Hold-out | 210 SKUs | 0.0871; corr 0.587 | 0.0719; 0.686 | Loses to naive; plain logit wins |
| C/D — Regime | 1,680 product-weeks | 0.0053; corr $-0.0039$ | 0.0051 | Levels only; lift unpredicted |

The pattern — accurate levels, absent spikes — is what static discrete-choice theory predicts: nested logit is a *stock* model mapping current prices and fixed taste to shares, while promotions are *flow* phenomena driven by stockpiling, forward-buying, and anticipation, none of which enters $\delta_{jt}$. I therefore validate the model for structural uses — ranking elasticities, mapping substitution, setting long-run list prices — and invalidate it for forecasting event-level promotional lift.

## 6. Recommendation

**Asymmetric format pricing.** Hold baseline prices on powders ($-3.09$), which are inelastic relative to bars: broad powder discounting buys little volume while conceding margin on every unit, and implied Lerner margins ($1/\lvert E \rvert \approx 0.36$) leave room to price at or above current levels. Direct discount funding to bars ($-5.22$), where a 10% cut buys roughly 50% more volume — and because inter-nest cross-elasticities are about $0.0001$, a bar promotion will not erode powder margins. Rival bar price moves, felt at cross-elasticities up to $4.3$, should be matched within the bar line.

**Flavour portfolio.** Treat Salted Caramel ($-4.10$, cross-switching $0.317$) as a customer-acquisition loss leader: price efficiently pulls trial from competitors, but loyalty is shallow, so keep discounts shallow and short. Treat Mocha/Coffee ($-3.19$, lowest cross-substitution at $0.099$) as the margin-harvesting set that supports list-price increases and premium price-pack architecture. Birthday Cake powder ($-2.50$, the least elastic focal cell) completes the ladder as a premium-margin contributor rather than a promotional vehicle.

**Promotional calendar.** Plan calendars on the model's validated structural claims, not its invalidated forecasts: price cuts induce temporary switching *within* nests — concentrated in bars — and do not expand cross-category consumption. So (i) never discount competing same-nest brands in the same week, which would split a fixed demand pool at double margin cost; (ii) use format-exclusive deal periods so each promotion adds incremental nest volume; and (iii) measure lift with experimental price variation — A/B tests or hold-out stores — because a panel whose price coefficient is statistically zero cannot support model-based lift forecasting.

**Boundary conditions.** I treat all magnitudes as directional upper bounds: the price coefficient is identified without cost instruments, shares proxy reviews rather than units, and both hold-out exercises temper point predictions. Within those bounds — ranking SKUs, structuring price relationships, and allocating promotion depth by format — the estimates give a defensible quantitative foundation for the pricing decisions framed in Section 1.

## References

- Berry, S. (1994). Estimating Discrete-Choice Models of Product Differentiation. *RAND Journal of Economics*, 25(2), 242–262.
- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in Market Equilibrium. *Econometrica*, 63(4), 841–890.
- Dong, T. (2026). *Supply and Demand Functions: A Survey by Data Availability*. University at Albany. SSRN Working Paper 7399858. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7399858
- McFadden, D. (1974). Conditional Logit Analysis of Qualitative Choice Behavior. In P. Zarembka (Ed.), *Frontiers in Econometrics* (pp. 105–142). Academic Press.
- Train, K. (2009). *Discrete Choice Methods with Simulation* (2nd ed.). Cambridge University Press.
- Kaggle dataset: *Sports nutrition supplements* (daniyaljavaid). https://www.kaggle.com/datasets/daniyaljavaid/sports-nutrition-supplements (licence not stated by publisher).
- Kaggle dataset: *Supplement Sales Data* (zahidmughal2343). https://www.kaggle.com/datasets/zahidmughal2343/supplement-sales-data (Apache 2.0).

## Exhibits & Appendix

All figures and supplementary matrices referenced in Sections 1–6 are collected here (excluded from the page count).

### Figures

![Review mass distribution across flavour profiles and nests](output/figures/fig1_market_shares.png){width=60%}

*Exhibit E1: Review mass distribution across flavour profiles and product nests. Volume concentrates in Chocolate/Vanilla; Mocha, Birthday Cake, and Salted Caramel form the long tail.*

![Price per serving against log market share by nest](output/figures/fig2_price_share.png){width=60%}

*Exhibit E2: Price per serving against log market share by nest, showing the distinct price clusters of bars ($\$2.58$ mean) versus powders ($\$1.45$ mean).*

![Elasticity heatmap — Mocha/Coffee](output/figures/fig3_elasticity_heatmap_Mocha_Coffee.png){width=55%}

*Exhibit E3a: Own- and cross-price elasticity matrix — Mocha/Coffee (26 products). Deep-blue diagonal = own-elasticity; red blocks = intra-nest substitution.*

![Birthday Cake elasticity heatmap](output/figures/fig3_elasticity_heatmap_Birthday_Cake.png){width=55%}

*Exhibit E3b: Elasticity matrix — Birthday Cake (17 products).*

![Salted Caramel elasticity heatmap](output/figures/fig3_elasticity_heatmap_Salted_Caramel.png){width=55%}

*Exhibit E3c: Elasticity matrix — Salted Caramel (23 products). The bar block is uniformly deeper red, reflecting the higher switching rates of Exhibit ET2.*

![Counterfactual price sweep](output/figures/fig4_price_sweep.png){width=60%}

*Exhibit E4: Counterfactual price simulation — share diversion from Optimum Nutrition Mocha Cappuccino across a $\pm 40\%$ price sweep, all other prices fixed.*

![Actual versus fitted market shares](output/figures/fig5_fit.png){width=60%}

*Exhibit E5: Actual versus fitted market shares (in-sample, $N = 861$, normalised within flavour market; RMSPE $= 0.032$).*

![Level B hold-out: 210 held-out SKUs](output/figures/fig6a_holdout_skus.png){width=60%}

*Exhibit E6a: Level B — actual versus predicted shares for 210 held-out catalogue SKUs (RMSPE $= 0.0871$, correlation $= 0.587$; naive benchmark $0.0719$ / $0.686$).*

![Level C/D hold-out: promotional weeks](output/figures/fig6b_holdout_promotions.png){width=60%}

*Exhibit E6b: Level C/D — actual versus predicted shares across 1,680 held-out promotional product-weeks (RMSPE $= 0.0053$, MAPE $= 6.77\%$, correlation $= -0.0039$).*

![Promotional lift: actual versus predicted](output/figures/fig7_promo_lift.png){width=60%}

*Exhibit E7: Actual versus model-predicted promotional share lift (actual mean $+0.000032$, predicted mean $-0.000000005$, RMSE $= 0.0812$).*

### Supplementary tables

**Exhibit ET1 — The two datasets mapped to Level-2 fields**

| Feature | Catalogue (estimation) | Panel (validation) |
|---|---|---|
| Source | Kaggle: *Sports nutrition supplements* (bodybuilding.com listings; licence not stated) | Kaggle: *Supplement Sales Data* (Apache 2.0) |
| Units | 2,719 SKUs, one cross-section | 4,384 product-weeks (16 products × 274 weeks) |
| Shares | Review-mass proxy; $s_0 = 0.884$–$0.997$ | Weekly units → share (outside share fixed at 0.5) |
| Prices | Price and price per serving | Weekly retail price and % discount |
| Characteristics | Format, flavour, servings, ratings | Product name, category, platform, location |
| Ownership | Brand (108 source brands; 61 in sample) | Product identity |
| Role | Estimation, $N = 861$ | Validation, 1,680 promotional product-weeks |

**Exhibit ET2 — Focal flavour markets**

| | Mocha / Coffee | Birthday Cake | Salted Caramel |
|---|---|---|---|
| SKUs / brands | 26 / 15 | 17 / 9 | 23 / 14 |
| Powder / bar SKUs | 23 / 3 | 12 / 5 | 12 / 11 |
| Mean own-elasticity | $-3.194$ | $-3.236$ | $-4.095$ |
| — powder / bar | $-3.146$ / $-3.566$ | $-2.504$ / $-4.993$ | $-2.412$ / $-5.930$ |
| Intra-nest cross (mean) | $0.099$ | $0.276$ | $0.317$ |
| Inter-nest cross (mean) | $0.00013$ | $0.00013$ | $0.00009$ |
| Largest brand-level cross | $4.315$ (Convenient Nutrition $\leftarrow$ Quest) | $2.481$ (bars $\leftarrow$ Quest) | $3.538$ (Quest $\leftarrow$ OhYeah!) |

**Exhibit ET3 — Brand-level own elasticities and largest cross-exposures, Mocha/Coffee**

| Brand | Own (brand-level) | Largest cross-in (source; value) |
|---|---|---|
| Animal | $-2.459$ | Optimum Nutrition; 0.913 |
| Ascent | $-3.005$ | Optimum Nutrition; 0.913 |
| Body Nutrition | $-3.174$ | Optimum Nutrition; 0.913 |
| Bodybuilding.com Signature | $-1.298$ | Optimum Nutrition; 0.913 |
| Convenient Nutrition | $-4.282$ | Quest Nutrition; 4.315 |
| Dymatize | $-1.822$ | Optimum Nutrition; 0.913 |
| Isopure | $-4.059$ | Optimum Nutrition; 0.913 |
| Kaged Muscle | $-2.221$ | Optimum Nutrition; 0.913 |
| MuscleTech | $-2.780$ | Optimum Nutrition; 0.913 |
| Optimum Nutrition | $-1.085$ | MuscleTech; 0.317 |
| PEScience | $-2.399$ | Optimum Nutrition; 0.913 |
| Pro Supps | $-1.561$ | Optimum Nutrition; 0.913 |
| Quest Nutrition | $-1.063$ | Convenient Nutrition; 0.020 |
| Universal Nutrition | $-1.324$ | Optimum Nutrition; 0.913 |
| Vega | $-8.476$ | Optimum Nutrition; 0.913 |

### A. Flavour-profile market sizes

| Market | SKUs | Market | SKUs |
|---|---|---|---|
| Chocolate | 305 | Salted Caramel | 23 |
| Vanilla | 191 | Peanut Butter | 24 |
| Other | 116 | Mocha / Coffee | 26 |
| Cookies & Cream | 82 | Birthday Cake | 17 |
| Strawberry/Berry | 54 | Banana | 13 |
| | | Unflavored | 10 |

*Total $N = 861$ (715 powder, 146 bar); two single-nest markets (Banana, Unflavored) are powders only.*

### B. Full brand-level elasticity matrix — Mocha/Coffee

Rows: responding brand (aggregate share); columns: brand whose uniform price changes. Diagonal = brand-level own elasticity.

| | Animal | Ascent | Body Nutr. | Signature | Conv. Nutr. | Dymatize | Isopure | Kaged | MuscleTech | Opt. Nutr. | PEScience | Pro Supps | Quest | Universal | Vega |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Animal | $-2.459$ | 0.007 | 0.074 | 0.111 | 0.000 | 0.093 | 0.039 | 0.081 | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | 0.164 |
| Ascent | 0.057 | $-3.005$ | 0.074 | 0.111 | 0.000 | 0.093 | 0.039 | 0.081 | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | 0.164 |
| Body Nutrition | 0.057 | 0.007 | $-3.174$ | 0.111 | 0.000 | 0.093 | 0.039 | 0.081 | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | 0.164 |
| Bodybuilding.com Signature | 0.057 | 0.007 | 0.074 | $-1.298$ | 0.000 | 0.093 | 0.039 | 0.081 | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | 0.164 |
| Convenient Nutrition | 0.000 | 0.000 | 0.000 | 0.000 | $-4.282$ | 0.000 | 0.000 | 0.000 | 0.001 | 0.135 | 0.000 | 0.000 | 4.315 | 0.000 | 0.000 |
| Dymatize | 0.057 | 0.007 | 0.074 | 0.111 | 0.000 | $-1.822$ | 0.039 | 0.081 | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | 0.164 |
| Isopure | 0.057 | 0.007 | 0.074 | 0.111 | 0.000 | 0.093 | $-4.059$ | 0.081 | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | 0.164 |
| Kaged Muscle | 0.057 | 0.007 | 0.074 | 0.111 | 0.000 | 0.093 | 0.039 | $-2.221$ | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | 0.164 |
| MuscleTech | 0.057 | 0.007 | 0.074 | 0.111 | 0.000 | 0.093 | 0.039 | 0.081 | $-2.780$ | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | 0.164 |
| Optimum Nutrition | 0.057 | 0.007 | 0.074 | 0.110 | 0.000 | 0.092 | 0.039 | 0.081 | 0.317 | $-1.085$ | 0.010 | 0.009 | 0.009 | 0.031 | 0.164 |
| PEScience | 0.057 | 0.007 | 0.074 | 0.111 | 0.000 | 0.093 | 0.039 | 0.081 | 0.317 | 0.913 | $-2.399$ | 0.009 | 0.000 | 0.031 | 0.164 |
| Pro Supps | 0.057 | 0.007 | 0.074 | 0.111 | 0.000 | 0.093 | 0.039 | 0.081 | 0.317 | 0.913 | 0.010 | $-1.561$ | 0.000 | 0.031 | 0.164 |
| Quest Nutrition | 0.000 | 0.000 | 0.000 | 0.000 | 0.020 | 0.000 | 0.000 | 0.000 | 0.001 | 0.001 | 0.000 | 0.000 | $-1.063$ | 0.000 | 0.000 |
| Universal Nutrition | 0.057 | 0.007 | 0.074 | 0.111 | 0.000 | 0.093 | 0.039 | 0.081 | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | $-1.324$ | 0.164 |
| Vega | 0.057 | 0.007 | 0.074 | 0.111 | 0.000 | 0.093 | 0.039 | 0.081 | 0.317 | 0.913 | 0.010 | 0.009 | 0.000 | 0.031 | $-8.476$ |

### C. Additional output files

Full product-level elasticity matrices, brand matrices for Birthday Cake and Salted Caramel, per-SKU own elasticities, and hold-out prediction files are provided in `output/` (`elasticity_matrix_*.csv`, `brand_matrix_*.csv`, `own_elasticities.csv`, `panel_holdout_predictions.csv`, `estimation_summary.csv`, `reference_results.json`). Reproduction scripts: `scripts/nested_logit_analysis.m` (MATLAB/GNU Octave) and `scripts/reference_analysis.py` (Python), which independently reproduce all reported statistics to within $5\times10^{-4}$.
