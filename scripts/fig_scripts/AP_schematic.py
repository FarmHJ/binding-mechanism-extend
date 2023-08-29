import matplotlib.pyplot as plt
import myokit


data_dir = '../../simulation_data/background/'
fig_dir = '../../figures/background/'

# Load AP signal
AP_log = myokit.DataLog.load_csv(data_dir + 'APclamp.csv')

fig = plt.figure(figsize=(3, 3))
ax = fig.add_subplot(111)
ax.plot(AP_log.time(), AP_log['membrane.V'], 'k', linewidth=2)

ax.set_xlim(0, 350)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
plt.savefig(fig_dir + 'AP.svg', format='svg')
