from bokeh.layouts import gridplot
from bokeh.models import NumeralTickFormatter, HoverTool
from bokeh.plotting import figure, output_file, show
import pandas as pd


Palette = ['#2970b0', '#b92731']

def plot(df, x=None, y=None, title=None, group=None, x_axis_label=None, y_axis_label=None, hue=None, fmt=''):
    hover_tool = HoverTool(
        tooltips=[
            ('size', '@size'),
            ('matrix', '@filename'),
            ('value', '@{}{}'.format(y, fmt))
        ],
        formatters={
            'value': 'printf'
        }
    )

    plots = []
    for (key, subdf) in df.groupby([group]):
        i = 0
        p = figure(
            title='{} / {}'.format(key.capitalize(), title),
            y_axis_type='log',
        )

        # Plot each unique value in the hue param column with a different color
        for (_, df_group) in subdf.groupby([hue]):
            p.line(x=x, y=y, line_width=2, color=Palette[i], legend_label=df_group[hue][0], source=df_group)
            p.circle(x=x, y=y, size=10, color=Palette[i], legend_label=df_group[hue][0], source=df_group)
            i += 1

        # Customize figure aesthetics
        p.title.text_font_size = '18pt'
        p.background_fill_color = '#eaeaf2'

        p.xaxis.axis_label = 'Size'
        p.xaxis.axis_label_text_font_size = '18pt'
        p.xaxis.major_label_text_font_size = '12pt'
        p.xaxis[0].formatter = NumeralTickFormatter(format=',')
        p.xgrid.grid_line_color = '#ffffff'
        p.xgrid.grid_line_width = 2
        p.xgrid.minor_grid_line_color = '#ffffff'
        p.xgrid.minor_grid_line_alpha = 0.6

        p.yaxis.axis_label = title
        p.yaxis.axis_label_text_font_size = '18pt'
        p.yaxis.major_label_text_font_size = '12pt'
        p.ygrid.grid_line_color = '#ffffff'
        p.ygrid.grid_line_width = 2
        p.ygrid.minor_grid_line_color = '#ffffff'
        p.ygrid.minor_grid_line_alpha = 0.6

        p.legend.location = 'bottom_right'
        p.legend.click_policy = 'hide'
        p.add_tools(hover_tool)

        plots.append(p)

    return plots


# Results to parse
results = [
    'cpp/windows-output.csv',
    'cpp/unix-output.csv',
    'matlab/windows-output.csv',
    'matlab/unix-output.csv',
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
df = pd.DataFrame(pd.concat(frames)).sort_values('size')


## Build results per os plots
output_file('results_per_os.html')

# Build plots
memory = plot(df, x='size', y='memory_delta', title='Memory Usage (bytes)', group='os', hue='source', fmt='{0.000 b}')
time = plot(df, x='size', y='solve_time', title='Cholesky Resolution Time (seconds)', group='os', hue='source', fmt='{:}')
error = plot(df, x='size', y='relative_error', title='Relative Error', group='os', hue='source')

# Arrage plots in a grid
grid = gridplot([memory, time, error], sizing_mode='stretch_width')
show(grid)


## Build results per implementation plots
output_file('results_per_implementation.html')

# Build plots
memory = plot(df, x='size', y='memory_delta', title='Memory Usage (bytes)', group='source', hue='os', fmt='{0.000 b}')
time = plot(df, x='size', y='solve_time', title='Cholesky Resolution Time (seconds)', group='source', hue='os', fmt='{:}')
error = plot(df, x='size', y='relative_error', title='Relative Error', group='source', hue='os')

# Arrage plots in a grid
grid = gridplot([memory, time, error], sizing_mode='stretch_width')
show(grid)
