import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from matplotlib.patches import Patch

FONT_SIZE = 16
plt.rcParams.update({'font.size': FONT_SIZE})

def parse_to_long_form(file_path):
    data_list = []
    indices = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if re.search(r'^\s*\d+\s*\|', line):
                all_nums = re.findall(r'-?\d+', line)
                if not all_nums: continue
                
                trial = int(all_nums[0])
                vals = [int(n) for n in all_nums[1:]]
                
                for i, idx_val in enumerate(indices):
                    if i*2 + 1 < len(vals):
                        data_list.append({'Index': idx_val, 'Type': 'ct0', 'Value': vals[i*2]})
                        data_list.append({'Index': idx_val, 'Type': 'ct1', 'Value': vals[i*2+1]})
    return pd.DataFrame(data_list)

# The log to plot. The parser reads a run's own output, so a fresh
# indcpad_ctct log works as well as the stored one, which is the default.
import sys
log_path = sys.argv[1] if len(sys.argv) > 1 else 'paper_results/data-ctct-40x10.txt'
df_plot = parse_to_long_form(log_path)
if df_plot.empty:
    raise SystemExit("no trial rows found in %s" % log_path)
print("plotting %d trials from %s" % (df_plot['Value'].size // 20, log_path))

plt.figure(figsize=(10, 5.5))
unique_indices = sorted(df_plot['Index'].unique())
colors = sns.color_palette("husl", len(unique_indices))

# no median line: the two distributions are both centred at 0 and are told apart by their spread
ax = sns.boxplot(x='Index', y='Value', hue='Type', data=df_plot,
                 medianprops={'visible': False}, showfliers=False, legend=False)

for i, patch in enumerate(ax.patches):
    idx_group = i % len(unique_indices)
    color = colors[idx_group]
    
    if i < len(unique_indices):
        patch.set_facecolor(color)
        patch.set_alpha(0.9)
    else:
        patch.set_facecolor(color)
        patch.set_hatch('////')
        patch.set_alpha(0.6)

ax.axhline(0, color='black', linestyle='-', linewidth=2, zorder=10)
vmax = df_plot['Value'].abs().max()
ax.set_ylim(-1.12 * vmax, 1.45 * vmax)   # symmetric, with headroom for the legend
ax.set_xlim(-0.62, 9.72)  # a little margin on both sides

legend_elements = [
    Patch(facecolor='gray', alpha=0.9, label="$\\mathcal{D}_0$ (left)"),
    Patch(facecolor='gray', hatch='////', alpha=0.6, label="$\\mathcal{D}_1$ (right)")
]
plt.legend(handles=legend_elements, loc='upper center', ncol=2, frameon=False,
           fontsize=FONT_SIZE)

plt.xlabel("Coefficient index", fontsize=FONT_SIZE + 2)
plt.ylabel("Decryption result", fontsize=FONT_SIZE + 2)
plt.xticks(fontsize=FONT_SIZE)
plt.yticks(fontsize=FONT_SIZE)

plt.title("") 

plt.grid(axis='y', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.show()