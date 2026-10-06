"""
Reference implementation (Python) of the nested-logit (BLP Level-2 aggregate)
analysis used to generate every number and figure reported in the paper.

Part A : product catalogue -> nests (powder vs bar), flavour-profile markets,
                              NL estimation, own/cross-price elasticities,
                              brand-level substitution matrices, figures.
Part B : weekly panel      -> hold-out test on the promotional subset.
Part C : catalogue hold-out-> out-of-sample share fit for 20% held-out SKUs.

Run:  python3 scripts/reference_analysis.py
"""
import os
import re
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "output")
FIG = os.path.join(OUT, "figures")
for d in (OUT, FIG):
    os.makedirs(d, exist_ok=True)

RNG = np.random.default_rng(20240101)

POWDER_CATS = ["WHEY PROTEIN", "WHEY PROTEIN ISOLATE", "MICELLAR CASEIN PROTEIN",
               "PLANT PROTEIN", "LOW CARB PROTEIN"]
BAR_CATS = ["PROTEIN BARS"]

FLAVOR_RULES = [
    ("Birthday Cake", r"birthday"),
    ("Salted Caramel", r"salted caramel|caramel"),
    ("Mocha / Coffee", r"mocha|coffee|cappuccino|espresso|cafe|latte"),
    ("Cookies & Cream", r"cookie|oreo"),
    ("Chocolate", r"chocolate|cocoa|fudge|brownie|double rich|milk chocolate"),
    ("Vanilla", r"vanilla|creme|cream$"),
    ("Strawberry/Berry", r"strawberr|berry|raspberry"),
    ("Peanut Butter", r"peanut|pb&j"),
    ("Banana", r"banana"),
    ("Mint", r"\bmint\b"),
    ("Unflavored", r"unflavored|unflavoured|flavorless"),
]
FOCAL_MARKETS = ["Mocha / Coffee", "Birthday Cake", "Salted Caramel"]


def flavor_market(name):
    s = str(name).lower()
    for label, pat in FLAVOR_RULES:
        if re.search(pat, s):
            return label
    return "Other"


def money(x):
    try:
        return float(str(x).replace("$", "").replace(",", "").strip())
    except Exception:
        return np.nan


# ------------------------------------------------------------------ NL machinery
def nl_matrix(beta_price, sigma, prices, nests, s, s_g, s_jg):
    """Elasticity matrix E[j,m] = d ln s_j / d ln p_m for one market.

    Price enters utility linearly, so d delta_m / d p_m = beta_price and the
    nested-logit share derivatives d ln s_j / d delta_m are (Berry 1994):
      own          : (1 - sigma*s_jg)/(1-sigma) - s_jg*s_g
      cross, same  : -s_mg*(sigma/(1-sigma) + s_g)
      cross, other : -s_m
    Multiplying by p_m converts d ln s_j / d p_m into the standard
    elasticity d ln s_j / d ln p_m.
    """
    n = len(prices)
    one_m_sigma = 1.0 - sigma
    diag = (1.0 - sigma * s_jg) / one_m_sigma - s_jg * s_g
    E = np.zeros((n, n))
    for m in range(n):
        for j in range(n):
            if j == m:
                deriv = diag[j]
            elif nests[j] == nests[m]:
                deriv = -s_jg[m] * (sigma / one_m_sigma + s_g[m])
            else:
                deriv = -s[m]
            E[j, m] = beta_price * prices[m] * deriv
    return E


def nl_shares(delta, nests, sigma):
    """Closed-form nested-logit share mapping for one market.

    delta : vector of mean utilities (x*beta + FE), nests : nest ids,
    sigma : within-nest correlation parameter. Returns (s_inside, s0) where
    s_inside sums to 1-s0.
    """
    delta = np.asarray(delta, dtype=float)
    nests = np.asarray(nests)
    sigma = float(sigma)
    one_m_sigma = 1.0 - sigma
    v = delta / one_m_sigma
    uniq = np.unique(nests)
    lnD = {}
    for g in uniq:
        idx = np.where(nests == g)[0]
        mx = v[idx].max()
        lnD[g] = mx + np.log(np.exp(v[idx] - mx).sum())     # ln D_g
    # ln(1 + sum_g D_g^(1-sigma))
    terms = np.array([0.0] + [(one_m_sigma) * lnD[g] for g in uniq])
    mx = terms.max()
    lnDen = mx + np.log(np.exp(terms - mx).sum())
    s = np.zeros(len(delta))
    for g in uniq:
        idx = np.where(nests == g)[0]
        lnDg = lnD[g]
        for j in idx:
            s[j] = np.exp(v[j] - sigma * lnDg - lnDen)
    s0 = np.exp(-lnDen)
    return s, s0


# ------------------------------------------------------------------ Part A: data
def load_catalog():
    df = pd.read_csv(os.path.join(DATA, "sports_nutrition_supplements.csv"))
    df["price"] = df["Price"].map(money)
    df["pps"] = df["PricePerServing"].map(money)
    df["reviews"] = pd.to_numeric(df["Reviews"], errors="coerce").fillna(0)
    df["flavorrating"] = pd.to_numeric(df["FlavorRating"], errors="coerce")
    df["rating"] = pd.to_numeric(df["Rating"], errors="coerce")
    df["prod"] = df["Brand"].str.strip() + " | " + df["Name"].str.strip()
    df["nsku"] = df.groupby("prod")["prod"].transform("size")
    # a parent product's cumulative review mass is split equally over its SKUs
    df["mass"] = df["reviews"] / df["nsku"]
    df["s_all"] = df["mass"] / df["mass"].sum()
    df["market"] = df["Flavor"].fillna("NA").str.strip().map(flavor_market)
    # Nests: powder vs bar. The source catalogue lists several bar-form
    # products (Quest Bars, Protein Crisp Bars, ...) under powder categories,
    # so product-name form rules are applied inside the powder categories.
    is_bar_cat = df["Category"].isin(BAR_CATS)
    is_powder_cat = df["Category"].isin(POWDER_CATS)
    nm = df["Name"].fillna("").str.strip().str.lower()
    # substring form rules (portable to Octave, whose regexp lacks \b)
    bar_form = nm.str.contains("bar") | nm.str.contains("wafer") | (
        nm.str.contains("cookie") & ~nm.str.contains("shake"))
    df["nest"] = np.where(is_bar_cat | (is_powder_cat & bar_form), "bar",
                          np.where(is_powder_cat, "powder", ""))
    return df


def prepare_inside(df):
    inside = df[df["Category"].isin(POWDER_CATS + BAR_CATS) & (df["mass"] > 0)].copy()
    inside["s0"] = 1.0 - inside.groupby("market")["s_all"].transform("sum")
    inside["sg"] = inside.groupby(["market", "nest"])["s_all"].transform("sum")
    inside["s_jg"] = inside["s_all"] / inside["sg"]
    inside["dep"] = np.log(inside["s_all"] / inside["s0"])
    inside["lns_jg"] = np.log(inside["s_jg"])
    inside["lnpps"] = np.log(inside["pps"])
    return inside

def estimate_catalog(d, cols=("lnpps", "lns_jg"), cluster=True):
    mkt = pd.get_dummies(d["market"], prefix="mkt", drop_first=True, dtype=float)
    X = pd.concat([d[list(cols)], mkt], axis=1)
    X = sm.add_constant(X, has_constant="add")
    kwargs = {}
    if cluster:
        kwargs = {"cov_type": "cluster", "cov_kwds": {"groups": d["prod"]}}
    return sm.OLS(d["dep"], X).fit(**kwargs)


def market_elasticities(g, beta, sigma):
    """Product-level elasticity matrix for one flavour market."""
    g = g.copy().reset_index(drop=True)
    g["s_g"] = g.groupby("nest")["s_all"].transform("sum")
    E = nl_matrix(beta, sigma, g["pps"].values, g["nest"].values,
                  g["s_all"].values, g["s_g"].values, g["s_jg"].values)
    g["own"] = np.diag(E)
    return g, E


def brand_matrix(mkt, E):
    """Elasticity of brand a's aggregate share w.r.t. a uniform price
    change of +/-1%% in every SKU of brand c:

        Mb[a,c] = sum_{i in a} (s_i/S_a) * sum_{j in c} E[i,j]
    """
    brands = sorted(mkt["Brand"].unique())
    w = mkt["s_all"].values
    b = mkt["Brand"].values
    B = len(brands)
    Mb = np.zeros((B, B))
    for a, ba in enumerate(brands):
        ia = np.where(b == ba)[0]
        wa = w[ia] / w[ia].sum()
        for c, bc in enumerate(brands):
            ic = np.where(b == bc)[0]
            row = E[np.ix_(ia, ic)].sum(axis=1)   # effect on SKU i of brand-c price change
            Mb[a, c] = float(wa @ row)
    return brands, Mb


# ------------------------------------------------------------------ Part B: panel
POWDER_PANEL = {"Whey Protein", "BCAA", "Creatine", "Pre-Workout",
                "Electrolyte Powder", "Collagen Peptides"}


def load_panel(s0_assumed=0.5):
    p = pd.read_csv(os.path.join(DATA, "supplement_sales_weekly.csv"))
    p["Date"] = pd.to_datetime(p["Date"])
    p["market"] = p["Date"].dt.strftime("%Y-%m-%d")
    p["nest"] = np.where(p["Product Name"].isin(POWDER_PANEL), "powder", "solidose")
    tot = p.groupby("market")["Units Sold"].transform("sum")
    p["s"] = (1.0 - s0_assumed) * p["Units Sold"] / tot
    p["s_norm"] = p["s"] / p.groupby("market")["s"].transform("sum")
    p["s0"] = s0_assumed
    p["sg"] = p.groupby(["market", "nest"])["s"].transform("sum")
    p["s_jg"] = p["s"] / p["sg"]
    p["dep"] = np.log(p["s"] / p["s0"])
    p["lns_jg"] = np.log(p["s_jg"])
    p["lnp"] = np.log(p["Price"])
    return p


def main():
    report = {}

    # ================================================================== Part A
    cat = load_catalog()
    inside = prepare_inside(cat)
    model = estimate_catalog(inside, cols=("pps", "lns_jg"))
    beta, sigma = model.params["pps"], model.params["lns_jg"]
    d = inside.dropna(subset=["pps", "lns_jg", "dep"]).copy()
    # robustness: log price per serving
    m_log = estimate_catalog(d, cols=("lnpps", "lns_jg"))

    report["catalog"] = dict(
        n_obs=int(model.nobs), n_products=int(d["prod"].nunique()),
        n_markets=int(d["market"].nunique()), n_skus_source=int(len(cat)),
        beta_pps=float(beta), beta_se=float(model.bse["pps"]),
        beta_p=float(model.pvalues["pps"]),
        beta_lnpps_robust=float(m_log.params["lnpps"]),
        beta_lnpps_se=float(m_log.bse["lnpps"]),
        beta_lnpps_p=float(m_log.pvalues["lnpps"]),
        sigma_lnpps=float(m_log.params["lns_jg"]),
        sigma=float(sigma), sigma_se=float(model.bse["lns_jg"]),
        sigma_p=float(model.pvalues["lns_jg"]),
        r2=float(model.rsquared), adj_r2=float(model.rsquared_adj),
        s0_min=float(d["s0"].min()), s0_max=float(d["s0"].max()),
        powder_share=float(d.loc[d.nest == "powder", "s_all"].sum()),
        bar_share=float(d.loc[d.nest == "bar", "s_all"].sum()),
        mean_pps=float(d["pps"].mean()),
        mean_pps_powder=float(d.loc[d.nest == "powder", "pps"].mean()),
        mean_pps_bar=float(d.loc[d.nest == "bar", "pps"].mean()))

    # ---- own elasticities, all markets
    own_rows = []
    market_sheets = {}
    for name, g in d.groupby("market"):
        mkt, E = market_elasticities(g, beta, sigma)
        market_sheets[name] = (mkt, E)
        for _, r in mkt.iterrows():
            own_rows.append(dict(market=name, brand=r["Brand"], product=r["prod"],
                                 flavor=r["Flavor"], nest=r["nest"], pps=r["pps"],
                                 price=r["price"], share=r["s_all"], own=r["own"]))
    own = pd.DataFrame(own_rows)
    own.to_csv(os.path.join(OUT, "own_elasticities.csv"), index=False)
    report["own_elasticity"] = dict(
        mean=float(own["own"].mean()),
        mean_powder=float(own.loc[own.nest == "powder", "own"].mean()),
        mean_bar=float(own.loc[own.nest == "bar", "own"].mean()),
        min=float(own["own"].min()), max=float(own["own"].max()),
        p25=float(own["own"].quantile(.25)), p75=float(own["own"].quantile(.75)))

    # ---- focal flavour markets
    focal = {}
    for name in FOCAL_MARKETS:
        mkt, E = market_sheets[name]
        mkt = mkt.copy()
        mkt["label"] = (mkt["Brand"] + " - " + mkt["Flavor"].astype(str)).str.strip()
        labels = list(mkt["label"])
        stem = re.sub(r"[^a-zA-Z0-9]+", "_", name)
        pd.DataFrame(np.round(E, 4), index=labels, columns=labels).to_csv(
            os.path.join(OUT, "elasticity_matrix_%s.csv" % stem))
        brands, Mb = brand_matrix(mkt, E)
        pd.DataFrame(np.round(Mb, 4), index=brands, columns=brands).to_csv(
            os.path.join(OUT, "brand_matrix_%s.csv" % stem))
        off = ~np.eye(len(mkt), dtype=bool)
        same = np.array([[mkt["nest"][j] == mkt["nest"][m] for m in range(len(mkt))]
                         for j in range(len(mkt))]) & off
        cross = off & ~same
        focal[name] = dict(
            n_skus=int(len(mkt)), n_brands=int(mkt["Brand"].nunique()),
            n_powder=int((mkt["nest"] == "powder").sum()),
            n_bar=int((mkt["nest"] == "bar").sum()),
            own_mean=float(mkt["own"].mean()),
            own_powder=float(mkt.loc[mkt.nest == "powder", "own"].mean()),
            own_bar=float(mkt.loc[mkt.nest == "bar", "own"].mean()),
            cross_same_nest_mean=float(E[same].mean()),
            cross_cross_nest_mean=float(E[cross].mean()),
            cross_max=float(E[off].max()),
            own_min=float(mkt["own"].min()), own_max=float(mkt["own"].max()),
            brands=brands, brand_matrix=np.round(Mb, 3).tolist())

    # ---- substitution of focal flavours: cross-market brand aggregate table
    report["focal"] = focal

    # ============================================== Part C: catalogue hold-out
    d2 = d.copy()
    # deterministic, reproducible split identical to the MATLAB script
    d2["key"] = [sum((i + 1) * ord(c) for i, c in enumerate(p)) % 5
                 for p in d2["prod"]]
    d2["holdout"] = d2["key"] < round(5 * 0.20)
    tr, ho = d2[~d2["holdout"]].copy(), d2[d2["holdout"]].copy()
    m_tr = estimate_catalog(tr, cols=("pps", "lns_jg"), cluster=False)
    est = m_tr.params
    resid = d2.loc[~d2["holdout"], "dep"] - m_tr.fittedvalues
    # market level = constant + market dummy, plus mean training residual
    mkt_lvl = {}
    for m in tr["market"].unique():
        col = "mkt_" + m
        lvl = float(est["const"]) + (float(est[col]) if col in est.index else 0.0)
        sel = tr["market"].values == m
        if sel.any():
            lvl += float(resid.values[sel].mean())
        mkt_lvl[m] = lvl
    # market FE enter delta_hat so held-out products inherit their market level
    ho["delta_hat"] = (ho["market"].map(mkt_lvl).fillna(float(est["const"]))
                       + est["pps"] * ho["pps"])
    rows = []
    for name, g in ho.groupby("market"):
        # strict out-of-sample: share mapping uses sigma estimated on training SKUs only
        s_hat, s0_hat = nl_shares(g["delta_hat"].values, g["nest"].values,
                                  float(est["lns_jg"]))
        tmp = g[["prod", "market", "nest", "s_all"]].copy()
        tmp["s_pred"] = s_hat / s_hat.sum()
        tmp["s_actual"] = g["s_all"].values / g["s_all"].sum()
        rows.append(tmp)
    pr = pd.concat(rows)
    e = pr["s_pred"] - pr["s_actual"]
    # fair benchmark: average share by market computed from the *training* SKUs
    naive = tr.groupby("market")["s_all"].mean()
    pr["s_naive"] = pr["market"].map(naive)
    pr["s_naive"] = pr["s_naive"] / pr.groupby("market")["s_naive"].transform("sum")
    en_ho = pr["s_naive"] - pr["s_actual"]
    report["catalog_holdout"] = dict(
        n_holdout=int(len(pr)), n_markets=int(pr["market"].nunique()),
        rmspe=float(np.sqrt((e ** 2).mean())),
        mape=float((e.abs() / pr["s_actual"]).mean()),
        corr=float(pr["s_pred"].corr(pr["s_actual"])),
        naive_rmspe=float(np.sqrt((en_ho ** 2).mean())),
        naive_corr=float(pr["s_naive"].corr(pr["s_actual"])),
        beta_oos=float(est["pps"]), sigma_oos=float(est["lns_jg"]))

    # ================================================================== Part B
    panel = load_panel()
    promo_thr = 0.13
    # a week is a *promotional period* when its mean discount rate crosses the
    # threshold: the whole week (all 16 products) is then held out.
    week_disc = panel.groupby("market")["Discount"].mean()
    promo_weeks = set(week_disc[week_disc >= promo_thr].index)
    panel["promo"] = panel["market"].isin(promo_weeks)
    trp = panel[~panel["promo"]].copy()

    # (i) pooled NL with market (week) fixed effects on regular weeks
    mkt_dum = pd.get_dummies(trp["market"], prefix="w", drop_first=True, dtype=float)
    X = pd.concat([trp[["lnp", "lns_jg"]], mkt_dum], axis=1)
    X = sm.add_constant(X, has_constant="add")
    p_model = sm.OLS(trp["dep"], X).fit(cov_type="cluster",
                                        cov_kwds={"groups": trp["Product Name"]})
    # (ii) prediction model: product fixed effects (usable for unseen weeks)
    pdum = pd.get_dummies(trp["Product Name"], prefix="p", drop_first=True, dtype=float)
    X2 = pd.concat([trp[["lnp"]], pdum], axis=1)
    X2 = sm.add_constant(X2, has_constant="add")
    p_pred = sm.OLS(trp["dep"], X2).fit()
    sigma_pnl = p_model.params["lns_jg"]

    hold = panel[panel["promo"]].copy()
    hold["delta_hat"] = p_pred.params["const"] + p_pred.params["lnp"] * hold["lnp"]
    for c in p_pred.params.index:
        if c.startswith("p_"):
            hold.loc[hold["Product Name"] == c[2:], "delta_hat"] += p_pred.params[c]
    rows = []
    for wk, g in hold.groupby("market"):
        s_hat, _ = nl_shares(g["delta_hat"].values, g["nest"].values, sigma_pnl)
        tmp = g[["Product Name", "market", "nest", "s"]].copy()
        tmp["s_pred"] = s_hat / s_hat.sum()
        tmp["s_actual"] = g["s"].values / g["s"].sum()
        rows.append(tmp)
    pv = pd.concat(rows)
    naive_by_prod = trp.groupby("Product Name")["s"].mean()
    pv["s_naive"] = pv["Product Name"].map(naive_by_prod)
    pv["s_naive"] = pv["s_naive"] / pv.groupby("market")["s_naive"].transform("sum")
    em, en = pv["s_pred"] - pv["s_actual"], pv["s_naive"] - pv["s_actual"]
    pv.to_csv(os.path.join(OUT, "panel_holdout_predictions.csv"), index=False)

    # promo lift: actual vs predicted, promoted product-weeks
    base_actual = trp.groupby("Product Name")["s_norm"].mean()
    base_price = trp.groupby("Product Name")["lnp"].median()
    pv["s_base_actual"] = pv["Product Name"].map(base_actual)
    pv["lift_actual"] = pv["s_actual"] / pv["s_base_actual"] - 1.0
    # predicted shares if the same week were sold at regular (median) prices
    base = hold.copy()
    base["lnp"] = base["Product Name"].map(base_price)
    base["delta_hat"] = p_pred.params["const"] + p_pred.params["lnp"] * base["lnp"]
    for c in p_pred.params.index:
        if c.startswith("p_"):
            base.loc[base["Product Name"] == c[2:], "delta_hat"] += p_pred.params[c]
    rows_b = []
    for wk, g in base.groupby("market"):
        s_hat, _ = nl_shares(g["delta_hat"].values, g["nest"].values, sigma_pnl)
        tmp = g[["Product Name", "market"]].copy()
        tmp["s_base_pred"] = s_hat / s_hat.sum()
        rows_b.append(tmp)
    pb = pd.concat(rows_b)
    pv = pv.merge(pb, on=["Product Name", "market"], how="left")
    pv["lift_pred"] = pv["s_pred"] / pv["s_base_pred"] - 1.0

    report["panel"] = dict(
        n_obs=int(len(panel)), n_weeks=int(panel["market"].nunique()),
        n_products=int(panel["Product Name"].nunique()),
        promo_threshold=promo_thr, promo_rows=int(panel["promo"].sum()),
        promo_fraction=float(panel["promo"].mean()),
        beta_ln_price=float(p_model.params["lnp"]),
        beta_se=float(p_model.bse["lnp"]), beta_p=float(p_model.pvalues["lnp"]),
        sigma=float(sigma_pnl), sigma_se=float(p_model.bse["lns_jg"]),
        sigma_p=float(p_model.pvalues["lns_jg"]), r2_train=float(p_model.rsquared),
        corr_price_units=float(panel["Price"].corr(panel["Units Sold"])),
        holdout_n=int(len(pv)),
        holdout_rmspe=float(np.sqrt((em ** 2).mean())),
        holdout_mape=float((em.abs() / pv["s_actual"]).mean()),
        holdout_corr=float(pv["s_pred"].corr(pv["s_actual"])),
        naive_rmspe=float(np.sqrt((en ** 2).mean())),
        naive_corr=float(pv["s_naive"].corr(pv["s_actual"])),
        lift_actual_mean=float(pv["lift_actual"].mean()),
        lift_pred_mean=float(pv["lift_pred"].mean()),
        lift_actual_absmean=float(pv["lift_actual"].abs().mean()),
        lift_pred_absmean=float(pv["lift_pred"].abs().mean()),
        lift_corr=float(pv["lift_actual"].corr(pv["lift_pred"])),
        lift_rmse=float(np.sqrt(((pv["lift_actual"] - pv["lift_pred"]) ** 2).mean())))

    with open(os.path.join(OUT, "reference_results.json"), "w") as f:
        json.dump(report, f, indent=2, default=float)
    print(json.dumps(report, indent=2, default=float))


if __name__ == "__main__":
    main()
