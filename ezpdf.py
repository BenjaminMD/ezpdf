from ezplot import create_basic_plot, scatter_w_outline
import numpy as np
import re 


def plot_single_PDF(recipe, res):
    r, gobs, gcalc, gdiff, baseline, gr_composition = get_gr(recipe)

    xlabel = r'$r\,/\,\mathrm{\AA}$'
    ylabel = r'$G$($r$)$\,/\,\mathrm{\AA}^{-2}$'
    fig, ax = create_basic_plot(xlabel, ylabel)

    scatter_w_outline(ax, r, gobs, label='obs')

    ax.plot(r, gcalc, label='calc')
    ax.plot(r, gdiff + baseline, label='diff')

    ax.plot(r, np.zeros(len(r)) + baseline, '--', linewidth=1.0, color='black', zorder=-1)
    ax.plot(r, np.zeros(len(r)), '--',  color='black', linewidth=1.0, zorder=1)
 
    ax.set_xlim(r.min(), r.max())

    ax.legend(title=f'Rw = {res.rw:.2f}', ncol=3, loc='upper right')



def plot_PDF(fit, recipe, res):
    fig, ax = create_basic_plot('$r\,/\,\mathrm{\AA}$', '$G$($r$)$\,/\,\mathrm{\AA}^{-2}$')
    r, gobs, gcalc, gdiff, baseline, gr_composition = get_gr(recipe)

#    [ax.scatter(r, gobs, **kwargs) for kwargs in [dict(s=36, edgecolors='0.0', lw=1.5), dict(s=36, edgecolors='1.0', lw=0), dict(s=35, edgecolors='#f6a800', lw=0, alpha=0.1725)]]
#    ax.scatter(r, gobs, 36, "0.0", lw=1.5)
#    ax.scatter(r, gobs, 36, "1.0", lw=0)
#    ax.scatter(r, gobs, 35, "#f6a800", lw=0, alpha=0.1725)
#    ax.scatter([], [], 80, "#f6a800", lw=0, label='obs')
    
    ax = scatter_w_outline(ax, r, gobs, label='obs')

    ax.plot(r, gcalc, '-', label='calc', linewidth=1.0)
    ax.plot(r, gdiff + baseline, '-', label='diff', color='green')
    ax.plot(r, np.zeros(len(r)) + baseline, '--', linewidth=1.0, color='black', zorder=-1)
    ax.plot(r, np.zeros(len(r)), '--',  color='black', linewidth=1.0, zorder=1)
    ax.set_xlim(r.min(), r.max())
        
    ax.legend(title=f'Rw = {res.rw:.2f}', ncol=3, loc='upper right')
    gr_p = {}
    ord = {}
    for i, (name, gcalc_p) in enumerate(zip(fit.formulas.values(), gr_composition.values())):
        span = gcalc_p.min() - gcalc_p.max()
        ord[name] = span
        gr_p[name] = gcalc_p 
    
    for i, name in enumerate(sorted(ord, key=ord.get)):
        ax.plot(r, gr_p[name] + baseline * (i + 3))
        ax.text(r.max()*0.9, gr_p[name][r  > r.max()* 0.8].max() + baseline * (i + 3), f'{name}', va='bottom', ha='center')
    ax.set_ylim(gr_p[name].min()*1.5 + baseline * (i + 3) + baseline * 0.2, gcalc.max() * 1.5)
    ax.axhspan(baseline * 1.5, -1000, color='black', zorder=-5, alpha=0.15)
    ax.text(r.max()*0.5, baseline * 1.8, "Contributing Phases", va='top', ha='center', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2', alpha=0.5))
    ax.set_yticklabels(['']*10)
    


def get_gr(recipe):
    """
    Get the gr of a recipe and for each phase contribution
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
    for char in ['\)', '\(']:
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
