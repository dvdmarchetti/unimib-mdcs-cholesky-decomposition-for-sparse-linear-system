import math
import os


from adjustText import adjust_text
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def label_point(x, y, val, ax):
    points = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    texts = []
    for i, point in points.iterrows():
        texts.append(ax.text(point['x'], point['y'], str(point['val'])))
    adjust_text(texts)


# Plot style
sns.set_style('darkgrid')

# Results to parse
results = [
    'cpp/windows-output.csv',
    'cpp/unix-output.csv',
    'matlab/windows-output.csv',
]

# Build single csv file
frames = []
for filename in results:
    frame = pd.read_csv(filename)
    # Add source column
    path = filename.split('/')
    frame['source'] = path[0]
    frame['os'] = path[1].split('-')[0]
    frames.append(frame)

# Concatenate all results
data = pd.concat(frames)
print(data)

metrics = [
    {
        'df': pd.DataFrame({ 'size': data['size'], 'metric': 'memory', 'value': data['memory_delta'], 'source': data['source'], 'os': data['os'] }),
        'label': 'Process Memory (in Bytes)'
    },
    {
        'df': pd.DataFrame({ 'size': data['size'], 'metric': 'error', 'value': data['relative_error'], 'source': data['source'], 'os': data['os'] }),
        'label': 'Relative Error'
    },
    {
        'df': pd.DataFrame({ 'size': data['size'], 'metric': 'time', 'value': data['solve_time'], 'source': data['source'], 'os': data['os'] }),
        'label': 'Cholesky Resolution Time (in seconds)',
    }
]

# Print faceted grid plots (3x2)
dfs = pd.concat([metrics[0]['df'], metrics[1]['df'], metrics[2]['df']])
g = sns.FacetGrid(dfs, col='os', row='metric', hue='source', palette="Set1")
g = g.map(sns.lineplot, 'size', 'value', marker='o')
g.add_legend().set_xlabels('Size')

for i, ax in enumerate(g.axes.flat):
    metric = metrics[int(i/2)]
    ax.set_ylabel(metric['label'])
    ax.set_yscale('log')

# Plot with annotations
fig, axes = plt.subplots(3,1)
fig.set_figheight(15)
fig.set_figwidth(15)
axes = axes.reshape(-1)
[ax.set(yscale='log') for ax in axes]
for i, metric in enumerate(metrics):
    ax = sns.lineplot(x='size', y='value', hue='os', palette='Set1', data=pd.melt(metric['df'], ['size', 'source', 'os', 'metric']), marker="o", ax=axes[i])

    ax.set_xlabel('Size')
    ax.set_ylabel(metric['label'])

    # label_point(metric['df']['size'], metric['df']['value'], metric['df']['value'], ax)

plt.show()
