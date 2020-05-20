from bokeh.layouts import gridplot
from bokeh.models import NumeralTickFormatter, HoverTool
from bokeh.plotting import figure, output_file, show
import pandas as pd


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
results = pd.concat(frames)
df = pd.DataFrame(results).sort_values('size')


Palette = ['#2970b0', '#b92731']

def plot(df, x=None, y=None, title=None, group=None, x_axis_label=None, y_axis_label=None, hue=None, unit=''):
    hover_tool = HoverTool(
        tooltips=[
            ('size', '@size'),
            ('matrix', '@filename'),
            ('value', '@'+y+'{,} '+unit)
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

        for (_, subdf_slice) in subdf.groupby([hue]):
            p.line(x=x, y=y, line_width=2, color=Palette[i], legend_label=subdf_slice[hue][0], source=subdf_slice)
            p.circle(x=x, y=y, size=10, color=Palette[i], legend_label=subdf_slice[hue][0], source=subdf_slice)
            i += 1

        p.title.text_font_size = '18pt'
        p.background_fill_color = '#eaeaf2'

        p.xaxis.axis_label = 'Size'
        p.xaxis.axis_label_text_font_size = '13pt'
        p.xaxis[0].formatter = NumeralTickFormatter(format=',')
        p.xgrid.grid_line_color = '#ffffff'
        p.xgrid.grid_line_width = 2
        p.xgrid.minor_grid_line_color = '#ffffff'
        p.xgrid.minor_grid_line_alpha = 0.5

        p.yaxis.axis_label = title
        p.yaxis.axis_label_text_font_size = '13pt'
        p.ygrid.grid_line_color = '#ffffff'
        p.ygrid.grid_line_width = 2
        p.ygrid.minor_grid_line_color = '#ffffff'
        p.ygrid.minor_grid_line_alpha = 0.5

        p.legend.location = 'bottom_right'
        p.legend.click_policy = 'hide'
        p.add_tools(hover_tool)

        plots.append(p)

    return plots

output_file('results.html')

memory = plot(df, x='size', y='memory_delta', title='Memory Usage (bytes)', group='os', hue='source', unit='Byte')
time = plot(df, x='size', y='solve_time', title='Cholesky Resolution Time (seconds)', group='os', hue='source', unit='s')
error = plot(df, x='size', y='relative_error', title='Relative Error', group='os', hue='source')

grid = gridplot(
    [memory, time, error],
    sizing_mode='stretch_width'
)

show(grid)
exit()
