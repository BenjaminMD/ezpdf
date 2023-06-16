from ezplot import create_basic_plot, scatter_w_outline
from dataclasses import dataclass, field
from typing import Dict, Union
import numpy.typing as npt
import numpy as np
import re


@dataclass
class G:
    recipe: object
    r: npt.NDArray = field(init=False)
    obs: npt.NDArray = field(init=False)
    calc: npt.NDArray = field(init=False)
    diff: npt.NDArray = field(init=False)
    baseline: float = field(init=False)
    composition: Dict[str, npt.NDArray] = field(init=False)
    zero: npt.NDArray = field(init=False)
    xlabel: str = r'$r\,/\,\mathrm{\AA}$'
    ylabel: str = r'$G(r)\,/\,\mathrm{\AA}^{-2}$'
    contrib: str = 'Conbtributing Phases:'
    legend: Dict[str, Union[str, int, float]] = field(default_factory=lambda: {
            'ncol': 3,
            'loc': 'upper right',
            'borderpad': 0.1,
            'labelspacing': 0.1,
            'columnspacing': 0.4,
            'handlelength': 1.0,
            'handletextpad': 0.1
        })

    def __post_init__(self):
        r, obs, calc, diff, baseline, composition = get_gr(self.recipe)
        norm = lambda y: normalize(-1, 1, obs, y)
        self.r = r
        self.obs = norm(obs)
        self.calc = norm(calc)
        self.diff = norm(diff) + norm(baseline)
        self.baseline = norm(baseline)
        self.composition = composition
        self.zero = norm(0)


def normalize(lower_bound, upper_bound, y_ref, y):
    """
    Normalizes y to y_ref between lower_bound and upper_bound
    """
    y_min = y_ref.min()
    y_max = y_ref.max()
    print(y_min, y_max, "y_min, y_max")
    y_norm = (upper_bound - lower_bound) * (y - y_min)
    y_norm = y_norm / (y_max - y_min) + lower_bound
    print(y_norm.min(), y_norm.max(), "y_norm.min(), y_norm.max()")
    return y_norm


def pdf_base(g: G, res):
    print(normalize(-1, 1, np.array([-10, 10]), np.array([-10, 10])))

    fig, ax = create_basic_plot(g.xlabel, g.ylabel)

    scatter_w_outline(ax, g.r, g.obs, label='obs', color="#C9B1D0")

    ax.plot(g.r, g.calc, label='calc', c="#4B0082")
    ax.plot(g.r, g.diff, label='diff', c="#BDBDBD")

    ax.axhline(y=g.baseline + g.zero,  linewidth=0.5, color='k', zorder=-1)
    ax.axhline(y=g.zero, linewidth=0.5, color='k', zorder=-1)

    ax.set_xlim(g.r.min(), g.r.max())
    ax.legend(title=f'Rw = {res.rw:.2f}', **g.legend)
    return fig, ax


def plot_PDF(fit, recipe, res):
    g = G(recipe)
    fig, ax = pdf_base(g, res)
    if len(g.composition) == 1:
        return fig, ax
    gr_p = {}
    ord: Dict[str, float] = {}
    name_gcal_kv = zip(fit.formulas.values(), g.composition.values())
    for i, (name, gcalc_p) in enumerate(name_gcal_kv):
        span = gcalc_p.min() - gcalc_p.max()
        ord[name] = span
        gr_p[name] = gcalc_p

    for i, name in enumerate(sorted(ord, key=ord.get)):  # type: ignore
        ax.plot(g.r, gr_p[name] + g.baseline * (i + 3))

        x_pos = g.r.max()*0.9
        y_pos = gr_p[name][g.r > g.r.max() * 0.8].max() + g.baseline * (i + 3)
        ax.text(x_pos, y_pos, f'{name}', va='bottom', ha='center')

    i = len(gr_p)
    name = list(gr_p.keys())[-1]
    y_min = gr_p[name].min()*1.5 + g.baseline * (i + 3) + g.baseline * 0.2
    y_max = g.calc[g.r > 1/3 * g.r.max()].max() * 2.0
    ax.set_ylim(y_min, y_max)
    ax.set_ylim(-1.5, 1.1)

    x_pos = g.r.max()*0.5
    y_pos = g.baseline * 1.8
    bbox = dict(facecolor='white',
                edgecolor='black',
                boxstyle='round,pad=0.2',
                alpha=0.5)
    ax.text(x_pos, y_pos, g.contrib, va='top', ha='center', bbox=bbox)

    ax.set_yticklabels(['']*10)
    ax.axhspan(g.baseline * 1.5, -1000, color='k', zorder=-5, alpha=0.15)
    return fig, ax


def get_gr(recipe):
    """
    g.t the gr of a recipe and for each phase contribution
    returns:
    - r: list of floats
    - gobs: list of floats
    - gcalc: list of floats
    - gdiff: list of floats
    - baseline: float
    - gr_composition: dict of list of floats
    """
    def remove_consecutive_duplicates(string, char):
        indices = [m.start() for m in re.finditer(char * 2, string)]
        if indices:
            for i in indices:
                string = string[:i] + string[i+1:]
            return remove_consecutive_duplicates(string, char)
        else:
            return string

    equation = recipe.PDF.getEquation()
    for char in [r'\)', r'\(']:
        equation = (remove_consecutive_duplicates(equation, char))

    prof = recipe._contributions['PDF'].profile
    r = prof.x
    gobs = prof.y
    gcalc = recipe._contributions['PDF'].evaluate()
    baseline = 1.10 * gobs.min()
    gdiff = gobs - gcalc

    gr_composition = {}
    for eq in equation.split(' + '):
        gr = recipe.PDF.evaluateEquation(eq[1:])
        gr_composition[eq[1:]] = gr

    return r, gobs, gcalc, gdiff, baseline, gr_composition
