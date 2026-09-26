from pathlib import Path
import csv
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader((HERE / 'data/khatri_2014_fig6_digitized.csv').open()))
x = [float(r['continuous_equivalent_d_per_day']) for r in rows]
y = [float(r['colony_diameter_mm']) for r in rows]
e = [float(r['colony_diameter_sd_mm']) for r in rows]
pred = 0.20143503943516178
fig, ax = plt.subplots(figsize=(5.2, 3.8))
ax.errorbar(x, y, yerr=e, fmt='o', color='#166534', capsize=3, label='Khatri Fig. 6 digitization')
ax.axhline(pred, color='#9a3412', ls='--', label='Frozen Zhang PBE\n(uniform withdrawal)')
ax.set_xlabel('Continuous-equivalent dilution rate, $D$ (d$^{-1}$)')
ax.set_ylabel('Mean colony diameter (mm)')
ax.set_xlim(0.04, 0.145); ax.set_ylim(0, 0.48)
ax.legend(frameon=False, fontsize=8); fig.tight_layout()
fig.savefig(HERE / 'results/khatri_frozen_validation.png', dpi=220)
fig.savefig(HERE / 'results/khatri_frozen_validation.pdf')
