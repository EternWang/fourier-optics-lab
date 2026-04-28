#!/usr/bin/env python3
"""Reproduce the core numerical results in the Fourier Optics Lab report.

This script starts from raw CSV measurements, performs regression and
uncertainty propagation, and writes both machine-readable summaries
(CSV/JSON) and publication-ready plots.

The calculations mirror Appendix B of the report: screen-angle
regression, camera calibration, slit cutoff, and Abbe-limit estimate.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
OUT = ROOT / "analysis" / "output"
BLUE = "#2F6B9A"
ORANGE = "#D97935"
GREEN = "#5B8C5A"
GRAY = "#4A5568"
LIGHT = "#EEF2F6"


@dataclass
class ScreenFitResult:
    slope_m: float  # tan(theta)
    slope_sigma: float
    theta_rad: float
    theta_sigma: float
    d_um: float
    d_sigma_um: float
    r2: float


@dataclass
class CameraResult:
    s_um_per_px: float
    s_sigma_um_per_px: float
    d_um: float
    d_sigma_um: float


@dataclass
class SlitResult:
    d_um: float
    d_sigma_um: float
    w_true_in: float
    delta_w_in: float
    delta_center_mm: float


@dataclass
class AbbeResult:
    NA: float
    dx_min_um: float


def ensure_dirs() -> None:
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)


def set_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 240,
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.titlesize": 14,
            "axes.labelsize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#D9DEE7",
            "grid.linewidth": 0.8,
            "grid.alpha": 0.75,
            "legend.frameon": False,
        }
    )


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def draw_card(ax: plt.Axes, x: float, y: float, w: float, h: float, title: str, body: str, color: str) -> None:
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        facecolor="white",
        edgecolor="#CBD5E1",
        linewidth=1.0,
    )
    ax.add_patch(box)
    title_y = y + h - 0.075 if body else y + h / 2
    title_va = "top" if body else "center"
    ax.text(x + 0.035, title_y, title, ha="left", va=title_va, fontsize=10.2, weight="bold", color=color)
    if body:
        ax.text(x + 0.035, y + h - 0.215, body, ha="left", va="top", fontsize=8.45, color="#172033", linespacing=1.18)


def through_origin_regression(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Fit y = m x through the origin and return (m, sigma_m, r2)."""
    n = len(x)
    if n < 2:
        raise ValueError("Need at least 2 points for regression.")

    m = float(np.sum(x * y) / np.sum(x * x))
    residuals = y - m * x
    sigma_m = math.sqrt(float(np.sum(residuals**2) / ((n - 1) * np.sum(x**2))))

    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum(y**2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return m, sigma_m, r2


def screen_angle_method(lam_nm: float = 550.0) -> ScreenFitResult:
    df = pd.read_csv(DATA_RAW / "screen_angle.csv")
    L = df["L_cm"].to_numpy(dtype=float)
    y = df["y_cm"].to_numpy(dtype=float)

    m, sigma_m, r2 = through_origin_regression(L, y)

    theta = math.atan(m)
    sigma_theta = sigma_m / (1.0 + m**2)

    lam = lam_nm * 1e-9
    d = lam / math.sin(theta)
    sigma_d = abs(lam * math.cos(theta) / (math.sin(theta) ** 2)) * sigma_theta

    return ScreenFitResult(
        slope_m=m,
        slope_sigma=sigma_m,
        theta_rad=theta,
        theta_sigma=sigma_theta,
        d_um=d * 1e6,
        d_sigma_um=sigma_d * 1e6,
        r2=r2,
    )


def camera_method() -> CameraResult:
    df = pd.read_csv(DATA_RAW / "camera_calibration.csv").iloc[0]

    dx_um = float(df["cal_spacing_mm"]) * 1e3
    dp = float(df["cal_pixels_px"])
    sig_dp = float(df["sigma_cal_pixels_px"])

    N = float(df["count_pixels_px"])
    P = float(df["periods_count_P"])
    sig_P = float(df["sigma_periods_count_P"])

    s = dx_um / dp
    sig_s = s * sig_dp / dp

    d = s * N / P
    sig_d = d * math.sqrt((sig_s / s) ** 2 + (sig_P / P) ** 2)

    return CameraResult(
        s_um_per_px=s,
        s_sigma_um_per_px=sig_s,
        d_um=d,
        d_sigma_um=sig_d,
    )


def slit_cutoff_method(d_cam_um: float, lam_nm: float = 550.0) -> SlitResult:
    df = pd.read_csv(DATA_RAW / "slit_cutoff.csv").iloc[0]

    f = float(df["f_obj_mm"]) * 1e-3
    w_in = float(df["w_min_in"])
    sig_w_in = float(df["sigma_w_min_in"])

    inch = 0.0254
    w = w_in * inch
    sig_w = sig_w_in * inch

    lam = lam_nm * 1e-9
    d = 2.0 * f * lam / w
    sig_d = d * (sig_w / w)

    d_cam = d_cam_um * 1e-6
    w_true = 2.0 * f * lam / d_cam
    w_true_in = w_true / inch
    delta_w_in = w_true_in - w_in
    delta_center_mm = abs(delta_w_in) * inch * 1e3 / 2.0

    return SlitResult(
        d_um=d * 1e6,
        d_sigma_um=sig_d * 1e6,
        w_true_in=w_true_in,
        delta_w_in=delta_w_in,
        delta_center_mm=delta_center_mm,
    )


def abbe_limit(lam_nm: float = 550.0) -> AbbeResult:
    df = pd.read_csv(DATA_RAW / "abbe_params.csv").iloc[0]
    D = float(df["D_mm"]) * 1e-3
    f = float(df["f_mm"]) * 1e-3
    n = float(df["n"])

    lam = lam_nm * 1e-9
    NA = n * (D / (2.0 * f))
    dx_min = 0.61 * lam / NA
    return AbbeResult(NA=NA, dx_min_um=dx_min * 1e6)


def plot_screen_fit(screen: ScreenFitResult) -> None:
    set_plot_style()
    df = pd.read_csv(DATA_RAW / "screen_angle.csv")
    L = df["L_cm"].to_numpy(dtype=float)
    y = df["y_cm"].to_numpy(dtype=float)
    xerr = df["sigma_L_cm"].to_numpy(dtype=float)
    yerr = df["sigma_y_cm"].to_numpy(dtype=float)

    L_line = np.linspace(0, L.max() * 1.05, 200)
    y_line = screen.slope_m * L_line

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.errorbar(
        L,
        y,
        xerr=xerr,
        yerr=yerr,
        fmt="o",
        ms=6,
        color=BLUE,
        ecolor="#7FA7C7",
        capsize=3,
        label="Measured separation",
    )
    ax.plot(L_line, y_line, color=ORANGE, lw=2.4, label=f"Through-origin fit: y = {screen.slope_m:.4f} L")
    ax.set_xlabel("Screen distance L (cm)")
    ax.set_ylabel("First-order separation y (cm)")
    ax.set_title("Screen-angle regression")
    ax.text(
        0.03,
        0.95,
        f"Estimated grating period\n"
        f"d = {screen.d_um:.2f} +/- {screen.d_sigma_um:.2f} um\n"
        f"R^2 = {screen.r2:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": LIGHT, "edgecolor": "#CBD5E1"},
    )
    ax.legend(loc="lower right")
    save_figure(fig, OUT / "screen_fit.png")

    fig, ax = plt.subplots(figsize=(7.2, 3.3))
    residuals = y - screen.slope_m * L
    ax.axhline(0, color=GRAY, linewidth=1)
    ax.errorbar(L, residuals, yerr=yerr, fmt="o", ms=5, color=BLUE, ecolor="#7FA7C7", capsize=3)
    ax.set_xlabel("Screen distance L (cm)")
    ax.set_ylabel("Residual y - mL (cm)")
    ax.set_title("Residual diagnostics")
    save_figure(fig, OUT / "screen_residuals.png")


def plot_uncertainty_budget(screen: ScreenFitResult, cam: CameraResult, slit: SlitResult) -> None:
    set_plot_style()
    methods = ["screen", "camera", "slit (random)"]
    sigmas = [screen.d_sigma_um, cam.d_sigma_um, slit.d_sigma_um]

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    bars = ax.bar(methods, sigmas, color=[BLUE, GREEN, ORANGE], width=0.62)
    ax.set_ylabel("Random uncertainty in d (um)")
    ax.set_title("Random uncertainty by measurement route")
    for bar, value in zip(bars, sigmas):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.006, f"{value:.3f}", ha="center", va="bottom")
    save_figure(fig, OUT / "random_uncertainty_budget.png")


def plot_grating_method_comparison(
    screen: ScreenFitResult, cam: CameraResult, slit: SlitResult
) -> None:
    set_plot_style()
    methods = ["Screen angle", "Camera calibration", "Slit cutoff"]
    estimates = np.array([screen.d_um, cam.d_um, slit.d_um], dtype=float)
    sigmas = np.array([screen.d_sigma_um, cam.d_sigma_um, slit.d_sigma_um], dtype=float)
    ypos = np.arange(len(methods))

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    ax.axvspan(cam.d_um - cam.d_sigma_um, cam.d_um + cam.d_sigma_um, color=GREEN, alpha=0.12, label="Camera +/- 1 sigma")
    ax.errorbar(
        estimates,
        ypos,
        xerr=sigmas,
        fmt="o",
        color=BLUE,
        ecolor="#7FA7C7",
        capsize=5,
        markersize=8,
    )
    ax.axvline(
        cam.d_um,
        color=GREEN,
        linestyle="--",
        linewidth=1.4,
        label="camera estimate",
    )
    ax.set_yticks(ypos, methods)
    ax.set_xlabel("Estimated grating period d (um)")
    ax.set_title("Independent grating-period estimates")
    ax.set_xlim(9.35, 11.05)
    for x, yv, sigma in zip(estimates, ypos, sigmas):
        ax.text(x + sigma + 0.035, yv, f"{x:.2f} +/- {sigma:.2f}", va="center", color=GRAY)
    ax.annotate(
        "random sigma is small;\nsystematics dominate",
        xy=(slit.d_um + slit.d_sigma_um, ypos[2]),
        xytext=(10.62, ypos[2] + 0.36),
        arrowprops={"arrowstyle": "-", "color": GRAY},
        fontsize=9.5,
        color=GRAY,
    )
    ax.legend(loc="lower right")
    save_figure(fig, OUT / "grating_method_comparison.png")


def plot_research_snapshot(screen: ScreenFitResult, cam: CameraResult, slit: SlitResult, abbe: AbbeResult) -> None:
    set_plot_style()
    fig = plt.figure(figsize=(11.2, 5.5), facecolor="white")
    grid = fig.add_gridspec(2, 3, height_ratios=[0.92, 1.08], width_ratios=[1.0, 1.08, 0.92])

    ax_cards = fig.add_subplot(grid[0, :])
    ax_cards.axis("off")
    ax_cards.set_xlim(0, 1)
    ax_cards.set_ylim(0, 1)
    fig.suptitle("Fourier optics analysis workflow", x=0.04, y=0.985, ha="left", fontsize=17, weight="bold", color="#172033")
    fig.text(
        0.04,
        0.925,
        "Three independent measurement routes with structured uncertainty propagation and method comparison.",
        ha="left",
        fontsize=10.5,
        color=GRAY,
    )

    draw_card(
        ax_cards,
        0.02,
        0.08,
        0.28,
        0.68,
        "Screen-angle route",
        f"d = {screen.d_um:.2f} +/- {screen.d_sigma_um:.2f} um\n"
        f"through-origin R^2 = {screen.r2:.3f}\n"
        "diffraction spacing vs distance",
        BLUE,
    )
    draw_card(
        ax_cards,
        0.36,
        0.08,
        0.28,
        0.68,
        "Camera route",
        f"d = {cam.d_um:.2f} +/- {cam.d_sigma_um:.2f} um\n"
        f"scale = {cam.s_um_per_px:.3f} um/px\n"
        "calibrated object-plane pixels",
        GREEN,
    )
    draw_card(
        ax_cards,
        0.70,
        0.08,
        0.28,
        0.68,
        "Slit-cutoff route",
        f"d = {slit.d_um:.2f} +/- {slit.d_sigma_um:.2f} um\n"
        f"Abbe dx_min ~= {abbe.dx_min_um:.2f} um\n"
        "random error reported\nsystematics discussed",
        ORANGE,
    )

    ax_compare = fig.add_subplot(grid[1, :2])
    labels = ["Screen angle", "Camera calibration", "Slit cutoff"]
    estimates = np.array([screen.d_um, cam.d_um, slit.d_um], dtype=float)
    sigmas = np.array([screen.d_sigma_um, cam.d_sigma_um, slit.d_sigma_um], dtype=float)
    ypos = np.arange(len(labels))
    ax_compare.errorbar(estimates, ypos, xerr=sigmas, fmt="o", color=BLUE, ecolor="#7FA7C7", capsize=5, ms=8)
    ax_compare.axvline(cam.d_um, color=GREEN, ls="--", lw=1.5, label="camera estimate")
    ax_compare.set_yticks(ypos, labels)
    ax_compare.set_xlabel("Estimated grating period d (um)")
    ax_compare.set_title("Independent grating-period estimates", weight="bold", loc="left")
    ax_compare.set_xlim(9.35, 11.05)
    for x, y, sigma in zip(estimates, ypos, sigmas):
        ax_compare.text(x + sigma + 0.04, y, f"{x:.2f} +/- {sigma:.2f}", va="center", color=GRAY, fontsize=9.5)
    ax_compare.legend(loc="lower right")

    ax_flow = fig.add_subplot(grid[1, 2])
    ax_flow.axis("off")
    ax_flow.set_xlim(0, 1)
    ax_flow.set_ylim(0, 1)
    ax_flow.set_title("Artifacts", weight="bold", loc="left", pad=8)
    steps = ["raw CSV tables", "analytic propagation", "JSON + CSV outputs", "figures + PDF report"]
    y_positions = [0.82, 0.60, 0.38, 0.16]
    for idx, (step, y) in enumerate(zip(steps, y_positions)):
        draw_card(ax_flow, 0.06, y, 0.82, 0.13, f"{idx + 1}. {step}", "", [BLUE, GREEN, ORANGE, GRAY][idx])
        if idx < len(steps) - 1:
            ax_flow.annotate("", xy=(0.47, y - 0.015), xytext=(0.47, y - 0.07), arrowprops={"arrowstyle": "->", "color": "#94A3B8"})

    fig.tight_layout(rect=[0.03, 0.02, 0.99, 0.9])
    fig.savefig(OUT / "research_snapshot.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    ensure_dirs()

    screen = screen_angle_method(lam_nm=550.0)
    cam = camera_method()
    slit = slit_cutoff_method(d_cam_um=cam.d_um, lam_nm=550.0)
    abbe = abbe_limit(lam_nm=550.0)

    summary = {
        "screen": asdict(screen),
        "camera": asdict(cam),
        "slit": asdict(slit),
        "abbe": asdict(abbe),
    }

    (DATA_PROCESSED / "results.json").write_text(json.dumps(summary, indent=2))

    rows = [
        {"method": "screen", "d_um": screen.d_um, "sigma_rand_um": screen.d_sigma_um},
        {"method": "camera", "d_um": cam.d_um, "sigma_rand_um": cam.d_sigma_um},
        {"method": "slit", "d_um": slit.d_um, "sigma_rand_um": slit.d_sigma_um},
    ]
    pd.DataFrame(rows).to_csv(DATA_PROCESSED / "grating_results.csv", index=False)

    plot_screen_fit(screen)
    plot_uncertainty_budget(screen, cam, slit)
    plot_grating_method_comparison(screen, cam, slit)
    plot_research_snapshot(screen, cam, slit, abbe)

    print("=== Fourier Optics Lab: Reproduced results ===")
    print(f"Screen-angle: d = {screen.d_um:.2f} +/- {screen.d_sigma_um:.2f} um (random), R^2 = {screen.r2:.3f}")
    print(f"Camera:       d = {cam.d_um:.2f} +/- {cam.d_sigma_um:.2f} um (random)")
    print(f"Slit cutoff:  d = {slit.d_um:.2f} +/- {slit.d_sigma_um:.2f} um (random)")
    print(f"Abbe limit:   NA ~= {abbe.NA:.3f}, dx_min ~= {abbe.dx_min_um:.2f} um")
    print()
    print("Outputs written to:")
    print(f"  {DATA_PROCESSED}")
    print(f"  {OUT}")


if __name__ == "__main__":
    main()
