import matplotlib.pyplot as plt
import myokit
import os

import modelling

data_dir = os.path.join(modelling.RESULT_DIR, 'background')
fig_dir = os.path.join(modelling.FIG_DIR, 'background')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Load AP signal
AP_log = myokit.DataLog.load_csv(os.path.join(data_dir, 'APclamp.csv'))

fig = plt.figure(figsize=(3, 3))
ax = fig.add_subplot(111)
ax.plot(AP_log.time(), AP_log['membrane.V'], 'k', linewidth=2)

ax.set_xlim(0, 350)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
plt.savefig(os.path.join(fig_dir, 'AP.svg'), format='svg')
