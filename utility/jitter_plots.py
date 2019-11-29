import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def jitter_plot(data, columns, y_ticks=None, x_label=None, y_label=None, figsize=None):
    """'data' is expeceted to be a dictionary {'category1':[val1, val2,...], 'category2':[...],...}. It can contain additional data, but only categories named in 'columns' will be displayed (in the given order)."""
    dfs = [pd.DataFrame({cat:data[cat]}) for cat in columns]
    df = pd.concat(dfs, axis=1) # Allows for unequal sizes 
    if figsize:
        plt.figure(figsize=figsize)
    ax = sns.barplot(data=df, ci="sd", errwidth=1.5, capsize=0.15, order=columns)
    ax = sns.swarmplot(data=df, size=5, edgecolor='black', linewidth=1.0, order=columns) # color='black'
    if y_ticks:
        ax.set_yticks(y_ticks)
    if x_label:
        ax.set_xlabel(x_label)
    if y_label:
        ax.set_ylabel(y_label)#, fontweight='bold') # invalidates italics
    sns.despine()
    plt.show()

# relative wolbachia counts
data_wol = {
'Control':[1.021276691, 1.155293362, 1.092554453, 0.730875494, 0.806550375, 0.442964197, 1.149042465, 1.601442964],
'Fosmidomycin':[0.49857347, 0.430636853, 0.608088411, 0.331144234, 0.492089925, 0.490840966, 1.401609492, 0.2037194, 0.296991048],
'MDL-29951':[0.09904602, 0.358616242, 0.145153822, 0.112707962, 0.121009992, 0.35886761, 0.294754371, 0.393283097],
'Tenofovir':[0.249472169, 0.79240825, 0.103185698, 0.387572833, 0.17110741, 0.365528726, 0.235914238, 0.102691923]
}
#  relative mf production
data_mf = {
'Control':[1.048824593, 0.940325497, 1.274864376, 1.03074141, 0.705244123, 1.264631043, 1.353689567, 0.801526718, 0.819338422, 0.605597964],
'Fosmidomycin':[0.587703436, 0.741410488, 0.623869801, 0.623869801, 0.777576854, 0.51653944, 0.890585242, 0.623409669, 0.712468193, 0.480916031],
'MDL-29951':[0.922242315, 0.777576854, 0.940325497, 0.795660036, 0.931283906, 0.730279898, 0.748091603, 1.193384224, 0.890585242, 0.445292621],
'Tenofovir':[0.605786618, 0.904159132, 0.714285714, 0.533453888, 0.867992767, 0.765903308, 0.694656489, 0.338422392, 0.409669211, 0.694656489]
}

figsize = (5.0, 4.5)
columns = ['Control', 'Fosmidomycin', 'MDL-29951', 'Tenofovir']
y_ticks = [0.0, 0.5, 1.0, 1.5]
x_label = 'Drug treatment'
#y_label = 'Relative $\it{Wolbachia}$ counts per worm'
y_label = 'Relative microfilaria released per worm'
jitter_plot(data_mf, columns, y_ticks=y_ticks, x_label=x_label, y_label=y_label, figsize=figsize)
