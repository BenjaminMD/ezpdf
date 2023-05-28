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
    ylabel: str = r'$g(r)\,/\,\mathrm{\AA}^{-2}$'
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
        self.r = r
        self.obs = obs
        self.calc = calc
        self.diff = diff
        self.baseline = baseline
        self.composition = composition
        self.zero = np.zeros(len(r))


def pdf_base(g: G, res):

    fig, ax = create_basic_plot(g.xlabel, g.ylabel)

    scatter_w_outline(ax, g.r, g.obs, label='obs')

    ax.plot(g.r, g.calc, label='calc')
    ax.plot(g.r, g.diff + g.baseline, label='diff')

    ax.axhline(y=g.baseline,  linewidth=0.5, color='k', zorder=-1)
    ax.axhline(y=0, linewidth=0.5, color='k', zorder=-1)

    ax.set_xlim(g.r.min(), g.r.max())
    ax.legend(title=f'Rw = {res.rw:.2f}', **g.legend)
    return fig, ax


def plot_PDF(fit, recipe, res):
    g = G(recipe)
    fig, ax = pdf_base(g, res)
    if len(g.composition) == 1:
        return fig, ax
    gr_p = {}
    ord = {}
    for i, (name, gcalc_p) in enumerate(zip(fit.formulas.values(), g.composition.values())):
        span = gcalc_p.min() - gcalc_p.max()
        ord[name] = span
        gr_p[name] = gcalc_p 
    
    for i, name in enumerate(sorted(ord, key=ord.get)):
        ax.plot(g.r, gr_p[name] + g.baseline * (i + 3))
        ax.text(g.r.max()*0.9, gr_p[name][g.r  > g.r.max()* 0.8].max() + g.baseline * (i + 3), f'{name}', va='bottom', ha='center')
    ax.set_ylim(gr_p[name].min()*1.5 + g.baseline * (i + 3) + g.baseline * 0.2, g.calc[g.r > 1/3 * g.r.max()].max() * 2.0)
    ax.axhspan(g.baseline * 1.5, -1000, color='black', zorder=-5, alpha=0.15)
    ax.text(g.r.max()*0.5, g.baseline * 1.8, g.contrib, va='top', ha='center', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2', alpha=0.5))
    ax.set_yticklabels(['']*10)
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
    baseline = 1.35 * gobs.min()
    gdiff = gobs - gcalc

    gr_composition = {}
    for eq in equation.split(' + '):
        gr = recipe.PDF.evaluateEquation(eq[1:])
        gr_composition[eq[1:]] = gr

    return r, gobs, gcalc, gdiff, baseline, gr_composition
