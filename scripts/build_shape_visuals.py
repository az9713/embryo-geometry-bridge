"""Analytic teaching geometries, not fitted tissue fields or model equilibria."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT/'docs/assets'
H, L, A = 1., 2., .35


def style3d(ax, title):
    ax.set_title(title, fontsize=12, color='#edf3fc', pad=14)
    ax.set_facecolor('#102035')
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((.08, .15, .24, 1))
        axis.pane.set_edgecolor('#526a84')
        axis.label.set_color('#d6e3f4')
    ax.tick_params(colors='#c8d7ec', labelsize=8, pad=0)
    ax.set_xlabel('x / H', labelpad=0)
    ax.set_ylabel('y / H', labelpad=0)
    ax.set_zlabel('z / H', labelpad=0)
    ax.view_init(elev=29, azim=-58)
    ax.set_box_aspect((1, 1, .54))


def save(fig, stem):
    for suffix in ['png', 'svg']:
        fig.savefig(ASSETS/(stem+'.'+suffix), dpi=170,
                    metadata={'Creator': 'Research Project'} if suffix == 'svg'
                    else {'Software': 'Research Project'})
    plt.close(fig)


def main():
    plt.rcParams.update({'figure.facecolor': '#102035', 'axes.facecolor': '#102035',
                         'text.color': '#edf3fc', 'axes.labelcolor': '#d6e3f4',
                         'axes.edgecolor': '#58708b', 'xtick.color': '#c8d7ec',
                         'ytick.color': '#c8d7ec', 'font.size': 11, 'svg.fonttype': 'none'})
    rr, tt = np.meshgrid(np.linspace(0, L, 61), np.linspace(0, 2*np.pi, 145))
    xx, yy = rr*np.cos(tt), rr*np.sin(tt)
    phi = (1-(rr/L)**2)**2
    fields = [(0., 0.), (A, 0.), (A, 4.)]
    titles = ['Reference: flat finite-height patch',
              'Circular depression: axisymmetry allowed',
              'Twofold depression: axisymmetry excludes it']
    fig = plt.figure(figsize=(15, 5.1), constrained_layout=True)
    saved = {}
    for i, ((depth, beta), title) in enumerate(zip(fields, titles)):
        ax = fig.add_subplot(1, 3, i+1, projection='3d')
        z = H-depth*phi*(1+beta*(rr/L)**2*np.cos(2*tt))
        ax.plot_surface(xx, yy, z, cmap='viridis', norm=Normalize(.55, 1.05),
                        rcount=53, ccount=37, alpha=.97, linewidth=0)
        ax.plot_wireframe(xx, yy, z, rstride=12, cstride=8,
                          color='#d1efff', linewidth=.5, alpha=.6)
        ax.plot_wireframe(xx, yy, np.zeros_like(z), rstride=24, cstride=15,
                          color='#537b9b', linewidth=.6, alpha=.65)
        for theta in np.linspace(0, 2*np.pi, 13)[:-1]:
            ax.plot([L*np.cos(theta)]*2, [L*np.sin(theta)]*2, [0, H],
                    color='#7e9fb8', lw=.6, alpha=.65)
        style3d(ax, title)
        ax.set_zlim(0, 1.08)
        if i:
            ax.text2D(.05, .03, 'Same annular mean contrast: Q = -0.30625 H',
                      transform=ax.transAxes, color='#f1d096', fontsize=9)
        saved['z_'+str(i)] = z
    fig.suptitle('Explicit teaching maps: a fold is an embedding; one scalar does not determine its shape',
                 fontsize=14)
    save(fig, 'shape-maps-concrete')
    np.savez_compressed(ASSETS/'shape-maps-concrete-data.npz', X=xx, Y=yy, **saved)

    r = np.linspace(0, L, 401)
    f = H-A*(1-r*r/(L*L))**2
    fp = 4*A*r/(L*L)*(1-r*r/(L*L))
    fpp = 4*A/(L*L)*(1-3*r*r/(L*L))
    js = np.sqrt(1+fp*fp)
    kr = fpp/(1+fp*fp)**1.5
    kt = np.empty_like(r)
    kt[1:] = fp[1:]/(r[1:]*np.sqrt(1+fp[1:]**2))
    kt[0] = 4*A/(L*L)
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.5), constrained_layout=True)
    axes[0].plot(r/L, f/H, color='#83dfdc', lw=2.5)
    axes[0].axhline(1, color='#c3cfde', ls=':', lw=1)
    axes[0].annotate('Depth A = 0.35 H', xy=(0, .65), xytext=(.15, .85),
                     arrowprops={'arrowstyle': '->', 'color': '#ecc88c'}, color='#ecc88c')
    axes[0].set_ylabel('Apical height / H')
    axes[0].set_title('Position: a visible depression')
    axes[1].plot(r/L, js, color='#83dfdc', lw=2.5)
    axes[1].axhline(1, color='#c3cfde', ls=':', lw=1)
    axes[1].set_ylabel('Surface area ratio J_s')
    axes[1].set_title('Slope: increased surface area')
    axes[2].plot(r/L, kr*H, color='#efb870', lw=2.5, label='Meridional curvature')
    axes[2].plot(r/L, kt*H, color='#83dfdc', lw=2.5, label='Circumferential curvature')
    axes[2].axhline(0, color='#c3cfde', lw=1)
    axes[2].axvline(1/np.sqrt(3), color='#c3cfde', ls=':', lw=1)
    axes[2].set_ylabel('Principal curvature × H')
    axes[2].set_title('Second derivatives: bowl to saddle collar')
    legend = axes[2].legend(facecolor='#18314a', edgecolor='#58708b', fontsize=9)
    for label in legend.get_texts(): label.set_color('#edf3fc')
    for ax in axes:
        ax.set_xlabel('Reference radius R / L')
        ax.grid(alpha=.13)
    fig.suptitle('One explicit circular fold; illustrative geometry, not a force-balanced or biological result', fontsize=13)
    save(fig, 'shape-features-concrete')
    with (ASSETS/'shape-features-concrete.csv').open('w', newline='', encoding='utf8') as out:
        writer = csv.writer(out)
        writer.writerow(['R_over_L', 'height_over_H', 'area_ratio',
                         'meridional_curvature_times_H', 'hoop_curvature_times_H'])
        writer.writerows(zip(r/L, f/H, js, kr*H, kt*H))

    s, t = np.meshgrid(np.linspace(-.95, .95, 51), np.linspace(-.8, .8, 35))
    radius = .8
    fig = plt.figure(figsize=(11.4, 4.6), constrained_layout=True)
    for i in range(2):
        ax = fig.add_subplot(1, 2, i+1, projection='3d')
        x = s if i == 0 else radius*np.sin(s/radius)
        z = np.zeros_like(s) if i == 0 else radius*(1-np.cos(s/radius))
        ax.plot_surface(x, t, z, color='#59b7b5', alpha=.9, linewidth=0)
        ax.plot_wireframe(x, t, z, rstride=4, cstride=5, color='#def5ff', linewidth=.65)
        style3d(ax, 'Flat patch: k1 = k2 = 0' if i == 0 else 'Rolled patch: k1 = 1/Rc, k2 = 0')
        ax.set_zlim(-.05, .7)
        ax.text2D(.04, .03, 'Same intrinsic metric g = identity; same area ratio J_s = 1',
                  transform=ax.transAxes, fontsize=9, color='#f1d096')
    fig.suptitle('Stretching and bending are different: the same lengths can form different 3D shapes', fontsize=13)
    save(fig, 'shape-metric-versus-curvature')

    c = F(1, 4)**2
    inner = 1-c+c*c/3
    outer = (1-F(3, 4)**2)**2/3
    assert inner-outer == F(7, 8)
    assert F(1)-2*F(7, 20) == F(3, 10)
    assert F(1)-2*F(7, 20)*F(125, 108) == F(41, 216)
    receipt = {'status': 'PASS_ANALYTIC_TEACHING_FIGURES',
               'scope': 'Explicit analytic maps and elementary exact checks. No force solve, experimental fit or new research result.',
               'parameters': {'H': 1, 'L': 2, 'A': '7/20', 'asymmetric_beta': 4},
               'exact_checks': {'mean_phi_center': str(inner), 'mean_phi_ring': str(outer),
                                'height_contrast_over_A': '-7/8',
                                'axisymmetric_body_J_lower_bound': '3/10',
                                'asymmetric_body_J_lower_bound': '41/216'},
               'sampling_checks': {'area_ratio_min': float(js.min()),
                                   'max_slope': float(fp.max()),
                                   'center_curvatures_times_H': [float(kr[0]*H), float(kt[0]*H)]},
               'files': {p.name: sha256(p.read_bytes()).hexdigest()
                         for p in sorted(ASSETS.glob('shape-*-concrete*'))},
               'generator_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    receipt['files'].update({p.name: sha256(p.read_bytes()).hexdigest()
                             for p in sorted(ASSETS.glob('shape-metric-versus-curvature.*'))})
    (ROOT/'docs/SHAPE-VISUALS-RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf8')
    print(json.dumps({'status': receipt['status'], 'figures': 3,
                      'Q_over_A': '-7/8', 'positive_J_bounds': ['3/10', '41/216']}))


if __name__ == '__main__':
    main()
