import math
import os

import matplotlib
matplotlib.use('WebAgg')


from adjustText import adjust_text
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# def label_point(x, y, val, ax):
#     points = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
#     texts = []
#     for i, point in points.iterrows():
#         texts.append(ax.text(point['x'], point['y'], str(point['val'])))
#     # adjust_text(texts)


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
results = pd.concat(frames)

metrics = [
    {
        'df': pd.DataFrame({ 'filename': results['filename'], 'size': results['size'], 'metric': 'memory', 'value': results['memory_delta'], 'source': results['source'], 'os': results['os'] }),
        'label': 'Process Memory (in Bytes)'
    },
    {
        'df': pd.DataFrame({ 'filename': results['filename'], 'size': results['size'], 'metric': 'error', 'value': results['relative_error'], 'source': results['source'], 'os': results['os'] }),
        'label': 'Relative Error'
    },
    {
        'df': pd.DataFrame({ 'filename': results['filename'], 'size': results['size'], 'metric': 'time', 'value': results['solve_time'], 'source': results['source'], 'os': results['os'] }),
        'label': 'Cholesky Resolution Time (in seconds)',
    }
]

# Plot style
sns.set_style('darkgrid')

# Print faceted grid plots (3x2)
# Each column is a different operating system, each row is a different metric
data = pd.concat([metric['df'] for metric in metrics])
graph = sns.FacetGrid(data, col='source', row='metric', hue='os', palette="Set1", sharey=False, sharex=False)
graph = graph.map(sns.lineplot, 'size', 'value', marker='o', ci=None)
graph.add_legend().set_xlabels('Size')

for i, ax in enumerate(graph.axes.flat):
    metric = metrics[int(i/2)]
    ax.set_ylabel(metric['label'])
    ax.set_yscale('log')

graph = sns.FacetGrid(data, col='os', row='metric', hue='source', palette="Set1", sharey=False, sharex=False)
graph = graph.map(sns.lineplot, 'size', 'value', marker='o', ci=None)
graph.add_legend().set_xlabels('Size')

for i, ax in enumerate(graph.axes.flat):
    metric = metrics[int(i/2)]
    ax.set_ylabel(metric['label'])
    ax.set_yscale('log')

plt.show()
