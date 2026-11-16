#!/usr/bin/env python3
"""Typesets the display equations of the manual and the lattice report as pictures.

    python3 make_equations.py            # writes eq_*.png next to this file

The Claude Doc shows a latex block as a formula, but its PDF, Word and Markdown exports and plain
Markdown viewers show the TeX source. Pictures read the same everywhere. Needs pdflatex (standalone,
amsmath), pdftoppm and Pillow. All pictures have the same width, so the type has one size when a
viewer scales them to the width of the page.
"""
import os, subprocess, tempfile

from PIL import Image, ImageChops

EQUATIONS = {
    # manual, part 2: the label equation
    "eq_label_transport":
        r"\frac{\partial (\alpha_i \rho_i X)}{\partial t} + \nabla\cdot(\alpha_i \rho_i \mathbf{u}_i X)"
        r" - X\left[\frac{\partial (\alpha_i \rho_i)}{\partial t} + \nabla\cdot(\alpha_i \rho_i \mathbf{u}_i)\right] = 0",
    "eq_label_discrete":
        r"(\alpha\rho)^{o}\,\frac{X_P^{n}-X_P^{o}}{\Delta t} + \sum_{\text{inflow faces}} |F_f|\,(X_P^{n}-X_N^{n}) = 0",
    # manual: displacement form across periodic sides
    "eq_displacement":
        r"\frac{\partial (m\,\xi)}{\partial t} + \nabla\cdot(F\,\xi) - \xi\,\mathrm{contErr}"
        r" = -\frac{1}{V}\sum_f F_f\,(x_f - x_c)",
    # manual, part 3: the run
    "eq_shift":
        r"c^{\mathrm{new}}(c_1) = \frac{\sum_k w_k \, c^{\mathrm{old}}(c_{0,k})}{\sum_k w_k}",
    "eq_diffusion_swap":
        r"\Delta m = D \cdot \min(V_0\alpha_0,\, V_1\alpha_1) \cdot \tfrac{1}{2}\,(c_1 - c_0)",
    "eq_drift":
        r"V_{\mathrm{drift}} = \mathbf{u}_{\mathrm{drift}}\cdot\mathbf{A}\,\Delta t, \qquad"
        r" \Delta m = c(\text{donor}) \; V_{\mathrm{drift}} \; \alpha(\text{donor}) \; \rho",
    "eq_mass_integral":
        r"M = \sum_{\text{cells}} c \; V \; \alpha \; \rho",
    "eq_recurrence_norm":
        r"R_{ij} = \sum_k \left| n_k(i) - n_k(j) \right|",
    # report "Reactions in a lattice reactor": implicit diffusion step with the wall reaction
    "eq_implicit_step":
        r"V\,\frac{c^{\mathrm{new}}-c^{\mathrm{old}}}{\Delta t} \;=\; \sum_{\text{faces}} \frac{D\,A_f}{d_f}\,"
        r"\bigl(c^{\mathrm{new}}_{\mathrm{nb}}-c^{\mathrm{new}}\bigr) \;-\; V\,k_{\mathrm{wall}}\,c^{\mathrm{new}} \;+\; S",
}

DPI, WIDTH, MARGIN = 300, 1900, 36            # 1900 px = 16 cm of text at 300 dpi
DOC = r"""\documentclass[11pt,border=1pt]{standalone}
\usepackage{amsmath,amssymb}
\begin{document}
$\displaystyle %s$
\end{document}
"""


def build(name, tex, out_dir):
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "eq.tex"), "w") as f:
            f.write(DOC % tex)
        subprocess.run(["pdflatex", "-interaction=batchmode", "-halt-on-error", "eq.tex"], cwd=tmp, check=True,
                       stdout=subprocess.DEVNULL)
        subprocess.run(["pdftoppm", "-png", "-r", str(DPI), "-singlefile", "eq.pdf", "eq"], cwd=tmp, check=True)
        img = Image.open(os.path.join(tmp, "eq.png")).convert("RGB")
    box = ImageChops.difference(img, Image.new("RGB", img.size, "white")).getbbox()
    img = img.crop(box)
    if img.width > WIDTH - 2 * MARGIN:
        raise SystemExit("%s is wider than the page (%d px)" % (name, img.width))
    page = Image.new("RGB", (WIDTH, img.height + 2 * MARGIN), "white")
    page.paste(img, ((WIDTH - img.width) // 2, MARGIN))
    path = os.path.join(out_dir, name + ".png")
    page.save(path, optimize=True)
    return path, img.size


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    for name, tex in EQUATIONS.items():
        path, size = build(name, tex, here)
        print("%-24s %4d x %3d px" % (os.path.basename(path), size[0], size[1]))
