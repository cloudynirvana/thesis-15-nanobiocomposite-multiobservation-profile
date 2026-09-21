#!/usr/bin/env python3
"""Legal multi-channel profile versus an illegal sensing-and-killing merge.

The generator is a compartmental concept, not a patient and not a device.
Three known preparation loadings multiply three map coefficients:

    y_R = rho * s      Raman / SERS reporter
    y_E = lam  * e     exosome compartment score
    y_D = kap  * c     photothermal or ROS assay readout

The legal object keeps those coefficients on separate channels and refuses
to place them in a therapeutic parameter vector. The illegal object forces
rho = lam = kap = theta and calls theta the platform parameter.

A fused score m = rho * kap * s * c is recorded as a second illegal map.
Only the product is visible, so the profile of rho is flat once kap is
refit. A slice that freezes kap is not that profile.

Research computation only. Not a medical device, dose, or clinical tool.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
FIGDIR = ROOT / "figures"
FIGDIR.mkdir(exist_ok=True)
PROFILE_DIR = REPO / "profiles"

SEED = 20260921
SIGMA = 0.08
N_REP = 6
CHI2_95 = 3.841
FLAT_SPREAD = 0.5
RANK_TOL = 1e-8
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")

# Crossref records checked on 21 September 2026. Syntax is not enough.
VERIFIED_DOIS = {
    "10.1016/j.jare.2015.02.007",
    "10.1016/j.canlet.2011.06.022",
    "10.1126/science.aau6977",
    "10.1080/20013078.2018.1535750",
    "10.1038/s41467-019-11718-4",
    "10.1021/acs.chemrev.5b00265",
    "10.1002/smll.201600393",
    "10.1038/natrevmats.2017.24",
    "10.1038/s41565-018-0246-4",
    "10.1093/bioinformatics/btp358",
}

DISCLAIMER = (
    "Research object. Not a medical device, not clinical decision support, "
    "not a diagnostic, not a dose, and not a cure. No document DOI."
)

S = np.array([0.40, 0.80, 1.20, 1.60], dtype=float)
E = np.array([1.00, 1.00, 0.20, 0.20], dtype=float)
C = np.array([1.50, 1.10, 0.70, 0.40], dtype=float)
RHO = 1.70
LAM = 0.85
KAP = 0.42
REQUIRED_CHANNELS = ("ch_sers", "ch_exosome", "ch_damage")
REQUIRED_NONPARAM = "merged_sensing_killing_theta"


def scrub(value):
    if isinstance(value, dict):
        return {str(k): scrub(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [scrub(v) for v in value]
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if not np.isfinite(number) or abs(number) < 1e-10:
            return 0.0
        return float(f"{number:.8g}")
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def rank_of(eigenvalues: np.ndarray) -> int:
    vals = np.asarray(eigenvalues, dtype=float)
    vals = np.sort(np.abs(vals))[::-1]
    if vals.size == 0 or vals[0] <= 0.0:
        return 0
    return int(np.sum(vals > RANK_TOL * vals[0]))


def fisher_legal() -> np.ndarray:
    scale = N_REP / SIGMA**2
    return scale * np.diag([np.dot(S, S), np.dot(E, E), np.dot(C, C)])


def fisher_single(design: np.ndarray, index: int) -> np.ndarray:
    scale = N_REP / SIGMA**2
    out = np.zeros((3, 3), dtype=float)
    out[index, index] = scale * float(np.dot(design, design))
    return out


def fisher_product() -> np.ndarray:
    w = S * C
    scale = (N_REP / SIGMA**2) * float(np.dot(w, w))
    return scale * np.array(
        [[KAP**2, RHO * KAP], [RHO * KAP, RHO**2]],
        dtype=float,
    )


def fisher_illegal_scalar() -> float:
    scale = N_REP / SIGMA**2
    return float(scale * (np.dot(S, S) + np.dot(E, E) + np.dot(C, C)))


def theta_hat_from_means(y_r: np.ndarray, y_e: np.ndarray, y_d: np.ndarray) -> float:
    numer = float(np.dot(S, y_r) + np.dot(E, y_e) + np.dot(C, y_d))
    denom = float(np.dot(S, S) + np.dot(E, E) + np.dot(C, C))
    return numer / denom


def draw_means(rng: np.random.Generator) -> dict[str, np.ndarray]:
    def mean(level: np.ndarray) -> np.ndarray:
        noise = rng.normal(0.0, SIGMA, size=(N_REP, level.size))
        return (level + noise).mean(axis=0)

    return {
        "y_R": mean(RHO * S),
        "y_E": mean(LAM * E),
        "y_D": mean(KAP * C),
    }


def grid_multipliers() -> np.ndarray:
    base = np.geomspace(0.35, 2.8, 16)
    return np.unique(np.sort(np.concatenate([base, [1.0]])))


def classify(delta: np.ndarray) -> str:
    spread = float(np.max(delta) - np.min(delta))
    if spread < FLAT_SPREAD:
        return "flat"
    if float(delta[0]) > CHI2_95 and float(delta[-1]) > CHI2_95:
        return "closed"
    return "open"


def profile_record(values: np.ndarray, delta: np.ndarray, reference: float) -> dict:
    return {
        "reference": reference,
        "grid": values.tolist(),
        "delta_chi2": delta.tolist(),
        "endpoint_delta_chi2": [float(delta[0]), float(delta[-1])],
        "spread": float(np.max(delta) - np.min(delta)),
        "call": classify(delta),
    }


def legal_coordinate_profile(design: np.ndarray, truth: float, mean: np.ndarray) -> dict:
    """Profile one map coefficient. The other channels do not enter."""
    scale = N_REP / SIGMA**2
    hat = float(np.dot(design, mean) / np.dot(design, design))
    rss_hat = scale * float(np.sum((mean - hat * design) ** 2))
    multipliers = grid_multipliers()
    values = multipliers * truth
    delta = []
    for value in values:
        rss = scale * float(np.sum((mean - value * design) ** 2))
        delta.append(rss - rss_hat)
    return profile_record(values, np.asarray(delta), truth)


def absent_parameter_profile(truth: float) -> dict:
    """A channel that does not contain the symbol. The profile is identically flat."""
    multipliers = grid_multipliers()
    values = multipliers * truth
    delta = np.zeros(values.size, dtype=float)
    return profile_record(values, delta, truth)


def illegal_profile(y_r: np.ndarray, y_e: np.ndarray, y_d: np.ndarray) -> dict:
    hat = theta_hat_from_means(y_r, y_e, y_d)
    scale = N_REP / SIGMA**2

    def rss(theta: float) -> float:
        return scale * float(
            np.sum((y_r - theta * S) ** 2)
            + np.sum((y_e - theta * E) ** 2)
            + np.sum((y_d - theta * C) ** 2)
        )

    rss_hat = rss(hat)
    values = grid_multipliers() * hat
    delta = np.asarray([rss(float(v)) - rss_hat for v in values])
    record = profile_record(values, delta, hat)
    record["theta_hat"] = hat
    record["abs_error"] = {
        "rho": abs(hat - RHO),
        "lambda": abs(hat - LAM),
        "kappa": abs(hat - KAP),
    }
    record["signed_error"] = {
        "rho": hat - RHO,
        "lambda": hat - LAM,
        "kappa": hat - KAP,
    }
    return record


def product_profile_and_slice() -> dict:
    """Noiseless fused score. Profiling rho with kappa free is flat. The slice is not."""
    w = S * C
    pi = RHO * KAP
    mean = pi * w
    scale = N_REP / SIGMA**2
    multipliers = grid_multipliers()
    rho_values = multipliers * RHO

    # Profile: kappa = pi / rho restores the product.
    delta_profile = np.zeros(rho_values.size, dtype=float)

    # Slice: kappa frozen at the generating value.
    delta_slice = []
    rss_hat = 0.0
    for rho in rho_values:
        predicted = float(rho) * KAP * w
        rss = scale * float(np.sum((mean - predicted) ** 2))
        delta_slice.append(rss - rss_hat)
    return {
        "product_pi": pi,
        "profile_rho": profile_record(rho_values, delta_profile, RHO),
        "slice_rho_kappa_frozen": profile_record(
            rho_values, np.asarray(delta_slice), RHO
        ),
    }


def spectrum(name: str, matrix: np.ndarray, n_param: int) -> dict:
    vals = np.linalg.eigvalsh(np.atleast_2d(matrix))
    vals = np.sort(np.real(vals))[::-1]
    return {
        "name": name,
        "n_param": n_param,
        "eigenvalues": vals.tolist(),
        "rank": rank_of(vals),
    }


def validate(profile: dict) -> list[str]:
    failed: list[str] = []

    def need(ok: bool, rule: str) -> None:
        if not ok and rule not in failed:
            failed.append(rule)

    need(profile.get("schema_version") == "1.1.0-multiobs", "schema_version")
    need(bool(profile.get("profile_id")), "profile_id")
    extends = profile.get("extends") or {}
    need(
        extends.get("object") == "DiseaseProfile" and extends.get("schema_version") == "1.0.0",
        "extends_disease_profile_1_0_0",
    )
    need(bool(profile.get("disease_id")) and "patient" not in profile, "disease_class_not_patient")
    need(profile.get("disclaimer") == DISCLAIMER, "disclaimer")

    answers = profile.get("answers") or []
    answer_ids = {row.get("id") for row in answers if isinstance(row, dict)}
    need(answer_ids.issuperset({"Q1", "Q2", "Q3", "Q4"}), "four_answers")

    channels = profile.get("observation_channels") or []
    ids = [row.get("channel_id") for row in channels if isinstance(row, dict)]
    need(set(REQUIRED_CHANNELS).issubset(set(ids)), "required_channels")
    symbols = []
    coefficients = []
    for row in channels:
        if not isinstance(row, dict):
            need(False, "channel_record")
            continue
        symbols.append(row.get("symbol"))
        coefficients.append(row.get("map_coefficient"))
        need(bool(row.get("observes")), f"observes:{row.get('channel_id')}")
        need(bool(row.get("refusal")), f"refusal_text:{row.get('channel_id')}")
        need(row.get("parameter_status") == "not_in_theta", "channel_parameter_status")
    need(len(symbols) == len(set(symbols)), "symbol_not_unique")
    need(len(coefficients) == len(set(coefficients)), "map_coefficient_not_unique")

    theta = profile.get("theta")
    need(theta == [], "theta_not_empty")

    non_ids = [
        row.get("id")
        for row in (profile.get("non_parameters") or [])
        if isinstance(row, dict)
    ]
    need(REQUIRED_NONPARAM in non_ids, f"missing_non_parameter:{REQUIRED_NONPARAM}")

    hypotheses = profile.get("admitted_hypotheses") or []
    need(len(hypotheses) >= 1, "admitted_hypothesis")
    for row in hypotheses:
        if not isinstance(row, dict):
            need(False, "hypothesis_record")
            continue
        need(bool(row.get("falsifier")), "hypothesis_falsifier")
        need(
            row.get("parameter_status") == "forbidden_to_enter_theta",
            "hypothesis_parameter_status",
        )

    for row in profile.get("candidate_mechanisms") or []:
        if not isinstance(row, dict):
            continue
        if row.get("id") == "M-PLATFORM-WORKS":
            need(row.get("status") == "refused", "platform_claim_not_refused")
        if row.get("status") == "admitted":
            need(bool(row.get("falsifier")), "candidate_admitted_without_falsifier")

    citations = profile.get("citations") or []
    need(len(citations) >= 1, "citations")
    doi_seen = False
    for row in citations:
        if not isinstance(row, dict):
            need(False, "citation_record")
            continue
        doi = row.get("doi")
        if doi is None:
            continue
        doi_seen = True
        need(bool(DOI_RE.match(str(doi))), "doi_syntax")
        need(str(doi) in VERIFIED_DOIS, "doi_not_in_crossref_list")
    need(doi_seen, "at_least_one_doi")
    return failed


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise SystemExit(f"{path} did not contain a mapping")
    return data


def write_figures(noisy: dict, illegal_noisy: dict, illegal_quiet: dict, profiles: dict) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
        }
    )
    prep = np.arange(1, S.size + 1)
    sem = SIGMA / np.sqrt(N_REP)
    hat = illegal_noisy["theta_hat"]
    panels = [
        ("Raman reporter", noisy["y_R"], RHO * S, hat * S, "#1f4e79"),
        ("Exosome compartment", noisy["y_E"], LAM * E, hat * E, "#2f6f4e"),
        ("Damage assay", noisy["y_D"], KAP * C, hat * C, "#b85c38"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.3), sharex=True)
    for ax, (title, observed, legal, illegal, color) in zip(axes, panels):
        ax.errorbar(
            prep,
            observed,
            yerr=sem,
            fmt="o",
            color=color,
            ms=5,
            lw=0.8,
            capsize=2,
            label="Noisy mean",
        )
        ax.plot(prep, legal, color="#222222", lw=1.3, label="Legal map")
        ax.plot(prep, illegal, color="#222222", lw=1.1, ls="--", label="Illegal θ")
        ax.set_title(title)
        ax.set_xlabel("Preparation")
        ax.set_xticks(prep)
    axes[0].set_ylabel("Channel mean")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / "channels_by_preparation.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))
    left = axes[0]
    for key, label, color in (
        ("rho", r"Legal $\rho$", "#1f4e79"),
        ("lambda", r"Legal $\lambda$", "#2f6f4e"),
        ("kappa", r"Legal $\kappa$", "#b85c38"),
    ):
        rec = profiles["legal_noiseless"][key]
        mult = np.asarray(rec["grid"]) / rec["reference"]
        left.plot(mult, rec["delta_chi2"], color=color, lw=1.4, label=label)
    flat = profiles["raman_only_kappa"]
    mult = np.asarray(flat["grid"]) / flat["reference"]
    left.plot(mult, flat["delta_chi2"], color="#666666", lw=1.6, label=r"$\kappa$ on Raman only")
    prod = profiles["product"]["profile_rho"]
    mult = np.asarray(prod["grid"]) / prod["reference"]
    left.plot(
        mult,
        prod["delta_chi2"],
        color="#111111",
        lw=1.2,
        ls=":",
        label=r"Product profile of $\rho$",
    )
    left.axhline(CHI2_95, color="#888888", lw=0.7, ls="--")
    left.set_ylim(-0.5, 28)
    left.set_xlabel("Multiplier of the generating value")
    left.set_ylabel(r"Profile $\Delta\chi^{2}$")
    left.legend(frameon=False, fontsize=7.5)
    left.set_title("Separate channels")

    right = axes[1]
    rec = illegal_quiet
    right.plot(rec["grid"], rec["delta_chi2"], color="#6b3fa0", lw=1.5, label=r"Illegal $\theta$")
    right.axhline(CHI2_95, color="#888888", lw=0.7, ls="--")
    right.axvline(rec["theta_hat"], color="#6b3fa0", lw=0.8, ls="--")
    for value, name, color in (
        (RHO, r"$\rho$", "#1f4e79"),
        (LAM, r"$\lambda$", "#2f6f4e"),
        (KAP, r"$\kappa$", "#b85c38"),
    ):
        right.axvline(value, color=color, lw=0.9, label=name)
    right.set_ylim(-0.5, 28)
    right.set_xlabel(r"Value of $\theta$")
    right.set_title("One symbol for three maps")
    right.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / "profiles_legal_illegal.png")
    plt.close(fig)

    blocks = profiles["spectra"]
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    xpos = []
    heights = []
    colors = []
    labels = []
    cursor = 0
    palette = {
        "legal_joint": "#1f4e79",
        "raman_only": "#7aa0c4",
        "exosome_only": "#8fbfa8",
        "damage_only": "#e0b09a",
        "product_rho_kappa": "#111111",
        "illegal_scalar": "#6b3fa0",
    }
    pretty = {
        "legal_joint": "Legal",
        "raman_only": "Raman",
        "exosome_only": "Exosome",
        "damage_only": "Damage",
        "product_rho_kappa": "Product",
        "illegal_scalar": "Illegal θ",
    }
    for block in blocks:
        vals = np.asarray(block["eigenvalues"], dtype=float)
        largest = np.max(vals) if vals.size else 1.0
        norm = vals / largest if largest > 0 else vals
        for height in norm:
            xpos.append(cursor)
            heights.append(max(float(height), 1e-16))
            colors.append(palette[block["name"]])
            cursor += 1
        labels.append((np.mean(xpos[-len(vals) :]), pretty[block["name"]]))
        cursor += 0.8
    ax.bar(xpos, heights, color=colors, width=0.8)
    ax.set_yscale("log")
    ax.set_ylim(1e-16, 3)
    ax.axhline(RANK_TOL, color="#888888", lw=0.7, ls="--")
    ax.set_xticks([pos for pos, _ in labels])
    ax.set_xticklabels([name for _, name in labels])
    ax.set_ylabel("Eigenvalue / largest in the block")
    ax.set_title("Fisher spectra")
    fig.tight_layout()
    fig.savefig(FIGDIR / "fisher_spectra.png")
    plt.close(fig)


def example_content(profile: dict) -> list[str]:
    failed = []
    coeffs = [row["map_coefficient"] for row in profile["observation_channels"]]
    if set(coeffs) != {"rho", "lambda", "kappa"}:
        failed.append("example_coefficients")
    non_ids = {row["id"] for row in profile["non_parameters"]}
    required = {
        "merged_sensing_killing_theta",
        "raman_peak_as_kill_coefficient",
        "exosome_tropism_as_delivery_rate",
        "spr_peak_as_anticancer_efficacy",
        "laser_fluence_or_infusion",
        "theranostic_label",
    }
    if not required.issubset(non_ids):
        failed.append("example_non_parameters")
    return failed


def main() -> None:
    legal = load_yaml(PROFILE_DIR / "agnp_exosome_raman.profile.yaml")
    illegal = load_yaml(PROFILE_DIR / "illegal_merged_theta.profile.yaml")
    legal_failures = validate(legal)
    illegal_failures = validate(illegal)
    content_failures = example_content(legal)

    json_path = PROFILE_DIR / "agnp_exosome_raman.profile.json"
    json_path.write_text(json.dumps(legal, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    legal_json = json.loads(json_path.read_text(encoding="utf-8"))
    json_failures = validate(legal_json)

    quiet = {
        "y_R": RHO * S,
        "y_E": LAM * E,
        "y_D": KAP * C,
    }
    rng = np.random.default_rng(SEED)
    noisy = draw_means(rng)

    spectra = [
        spectrum("legal_joint", fisher_legal(), 3),
        spectrum("raman_only", fisher_single(S, 0), 3),
        spectrum("exosome_only", fisher_single(E, 1), 3),
        spectrum("damage_only", fisher_single(C, 2), 3),
        spectrum("product_rho_kappa", fisher_product(), 2),
        {
            "name": "illegal_scalar",
            "n_param": 1,
            "eigenvalues": [fisher_illegal_scalar()],
            "rank": 1,
        },
    ]

    profiles = {
        "legal_noiseless": {
            "rho": legal_coordinate_profile(S, RHO, quiet["y_R"]),
            "lambda": legal_coordinate_profile(E, LAM, quiet["y_E"]),
            "kappa": legal_coordinate_profile(C, KAP, quiet["y_D"]),
        },
        "legal_noisy": {
            "rho": legal_coordinate_profile(S, RHO, noisy["y_R"]),
            "lambda": legal_coordinate_profile(E, LAM, noisy["y_E"]),
            "kappa": legal_coordinate_profile(C, KAP, noisy["y_D"]),
        },
        "raman_only_kappa": absent_parameter_profile(KAP),
        "damage_only_rho": absent_parameter_profile(RHO),
        "product": product_profile_and_slice(),
        "spectra": spectra,
    }
    illegal_quiet = illegal_profile(quiet["y_R"], quiet["y_E"], quiet["y_D"])
    illegal_noisy = illegal_profile(noisy["y_R"], noisy["y_E"], noisy["y_D"])

    denom = float(np.dot(S, S) + np.dot(E, E) + np.dot(C, C))
    weights = {
        "rho": float(np.dot(S, S) / denom),
        "lambda": float(np.dot(E, E) / denom),
        "kappa": float(np.dot(C, C) / denom),
    }

    cramer = {}
    fisher = fisher_legal()
    for name, truth, idx in (("rho", RHO, 0), ("lambda", LAM, 1), ("kappa", KAP, 2)):
        se = float(1.0 / np.sqrt(fisher[idx, idx]))
        cramer[name] = {"se": se, "relative_se": se / truth}

    results = {
        "seed": SEED,
        "sigma": SIGMA,
        "n_rep": N_REP,
        "chi2_95_1df": CHI2_95,
        "flat_spread_below": FLAT_SPREAD,
        "rank_tol_relative": RANK_TOL,
        "truth": {"rho": RHO, "lambda": LAM, "kappa": KAP},
        "design": {"s": S.tolist(), "e": E.tolist(), "c": C.tolist()},
        "design_weights_in_illegal_average": weights,
        "cramer_rao_legal": cramer,
        "validator": {
            "legal_yaml_rules_failed": legal_failures,
            "legal_json_rules_failed": json_failures,
            "legal_example_content_failed": content_failures,
            "illegal_yaml_rules_failed": illegal_failures,
            "legal_accepted": legal_failures == [] and json_failures == [] and content_failures == [],
            "illegal_refused": illegal_failures != [],
        },
        "profile_doi_count": sum(1 for row in legal["citations"] if row.get("doi")),
        "spectra": spectra,
        "profiles_noiseless_calls": {
            name: profiles["legal_noiseless"][name]["call"] for name in ("rho", "lambda", "kappa")
        },
        "profiles_noisy_calls": {
            name: profiles["legal_noisy"][name]["call"] for name in ("rho", "lambda", "kappa")
        },
        "raman_only_kappa_call": profiles["raman_only_kappa"]["call"],
        "damage_only_rho_call": profiles["damage_only_rho"]["call"],
        "product_profile_rho_call": profiles["product"]["profile_rho"]["call"],
        "product_slice_rho_call": profiles["product"]["slice_rho_kappa_frozen"]["call"],
        "product_profile_rho_spread": profiles["product"]["profile_rho"]["spread"],
        "product_slice_endpoint_delta_chi2": profiles["product"]["slice_rho_kappa_frozen"][
            "endpoint_delta_chi2"
        ],
        "illegal_noiseless": {
            "theta_hat": illegal_quiet["theta_hat"],
            "call": illegal_quiet["call"],
            "endpoint_delta_chi2": illegal_quiet["endpoint_delta_chi2"],
            "abs_error": illegal_quiet["abs_error"],
            "signed_error": illegal_quiet["signed_error"],
        },
        "illegal_noisy": {
            "theta_hat": illegal_noisy["theta_hat"],
            "call": illegal_noisy["call"],
            "endpoint_delta_chi2": illegal_noisy["endpoint_delta_chi2"],
            "abs_error": illegal_noisy["abs_error"],
        },
        "legal_noiseless_endpoint_delta_chi2": {
            name: profiles["legal_noiseless"][name]["endpoint_delta_chi2"]
            for name in ("rho", "lambda", "kappa")
        },
    }

    write_figures(noisy, illegal_noisy, illegal_quiet, profiles)
    out = ROOT / "results.json"
    out.write_text(json.dumps(scrub(results), indent=2) + "\n", encoding="utf-8")

    expected_illegal = [
        "channel_parameter_status",
        "symbol_not_unique",
        "map_coefficient_not_unique",
        "theta_not_empty",
        "missing_non_parameter:merged_sensing_killing_theta",
        "hypothesis_parameter_status",
        "platform_claim_not_refused",
        "candidate_admitted_without_falsifier",
    ]
    problems = []
    if legal_failures or json_failures or content_failures:
        problems.append(f"legal object refused: {legal_failures} {json_failures} {content_failures}")
    if illegal_failures != expected_illegal:
        problems.append(f"illegal rules {illegal_failures} != {expected_illegal}")
    for block, expect_rank in (
        ("legal_joint", 3),
        ("raman_only", 1),
        ("exosome_only", 1),
        ("damage_only", 1),
        ("product_rho_kappa", 1),
        ("illegal_scalar", 1),
    ):
        found = next(row for row in spectra if row["name"] == block)
        if found["rank"] != expect_rank:
            problems.append(f"{block} rank {found['rank']} != {expect_rank}")
    if profiles["legal_noiseless"]["rho"]["call"] != "closed":
        problems.append("legal rho profile")
    if profiles["legal_noiseless"]["lambda"]["call"] != "closed":
        problems.append("legal lambda profile")
    if profiles["legal_noiseless"]["kappa"]["call"] != "closed":
        problems.append("legal kappa profile")
    if profiles["raman_only_kappa"]["call"] != "flat":
        problems.append("raman-only kappa")
    if profiles["product"]["profile_rho"]["call"] != "flat":
        problems.append("product profile")
    if profiles["product"]["slice_rho_kappa_frozen"]["call"] != "closed":
        problems.append("product slice")
    if illegal_quiet["call"] != "closed":
        problems.append("illegal profile")
    if min(illegal_quiet["abs_error"].values()) <= 0.15:
        problems.append("illegal theta too close to a true coefficient")
    if problems:
        raise SystemExit("checks failed: " + "; ".join(problems))

    print(f"wrote {out}")
    print(f"wrote {json_path}")
    print("illegal rules:", ", ".join(illegal_failures))
    print(
        "theta_hat noiseless",
        round(illegal_quiet["theta_hat"], 4),
        "errors",
        {k: round(v, 4) for k, v in illegal_quiet["abs_error"].items()},
    )


if __name__ == "__main__":
    main()
