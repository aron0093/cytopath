import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import matplotlib.cm as mc

# Make hexadecimal color per value
def hex_color(cmap, val):
    cmap = mc.get_cmap(cmap)
    color = [min(int(round(i*256)), 255) for i in cmap(val)[0:3]]
    return "#%02x%02x%02x" % tuple(color)
    #return '#{}'.format(''.join('{:02X}'.format(i) for i in color))

# Convolve values
def convolve(val, n):
    return np.convolve(val, np.ones(n)/n, mode="same")

# Standardize pandas series
def standardize_(val):
    return (val - val.mean())/(val.std()+1e-15)

# Base function
def make_radial_heatmap(ax, values, num_blocks, num_rings, ring_size, empty_size):

    radii = np.linspace(ring_size, 1, num_rings) + empty_size

    for count, radius in enumerate(radii):
        if values[count] == 'skip':
            colors = [(1,1,1)]*num_blocks
        else:
            cmap = values[count][1]
            colors = [hex_color(cmap, j) for j in values[count][0]]

        ax.pie([10/num_blocks]*num_blocks, radius = radius, colors = colors,
               counterclock=False, explode=[0.05]*num_blocks, startangle=0) 

    # Empty space for donut           
    ax.pie([100], radius = empty_size, colors = ["white"])  
    ax.axis("equal")

    return ax
   
def radial_heatmap(adata, var_names, end_points=None, trajectories=None, sortby='alignment', layer='Ms', 
                   col_color_key=None, col_colors=['Blues', 'Reds', 'Greens', 'Oranges','Purples'], 
                   color_map='viridis', n_convolve=10, standard_scale=True, mode='cut',
                   bins=100, figsize=(10, 10), ax=None, table_ax=None):
    
    # Check if data exists
    plotting_data = pd.DataFrame(np.zeros((5,5)), columns=['End point', 'Trajectory', 'Cell', 'Step', 'Allignment Score'])
    try:
        if sortby=='alignment':
            plotting_data = pd.DataFrame(adata.uns['trajectories']['cells_along_trajectories_each_step'])
        elif sortby=='pseudotime':
            plotting_data = pd.DataFrame(adata.uns['trajectories']['cells_along_trajectories'])       
    except:
        raise ValueError('Either cytopath was not run or the data has been corrupted.')

    # Add markers
    var_names_ = [var_name for var_name in var_names if var_name in adata.var_names]
    exp_data = pd.DataFrame(adata[:, var_names_].layers[layer], columns=var_names_)
    for var_name in var_names_:
        if standard_scale:
            exp_data[var_name] = standardize_(exp_data[var_name])
        plotting_data[var_name] = exp_data.loc[plotting_data.Cell.values, var_name].values
        
    # If end_point(s) specified
    if end_points is not None:
        if type(end_point) != list:
            end_point = [end_point]
        plotting_data = plotting_data.loc[plotting_data['End point'].isin(end_point)]

    # If trajectory specified
    if trajectories is not None:
        if type(trajectories) != list:
            trajectories = [trajectories]
        plotting_data = plotting_data.loc[plotting_data['Trajectory'].isin(trajectories)]

    # Run quantiles if specified
    if bins is not None:
        if mode=='quantile':
            plotting_data['binned_step'] = pd.qcut(plotting_data.Step, bins, labels=np.arange(bins))
        elif mode=='cut':
            plotting_data['binned_step'] = pd.cut(plotting_data.Step, bins, labels=np.arange(bins))
        dial_key = 'binned_step'
    else:
        dial_key = 'Step'
        bins=plotting_data.shape[0]

    # If categorical provided
    if col_color_key is not None:

        # Save category names for later
        col_color_categories = adata.obs[col_color_key].astype(str).unique()

        # Add col color data to plotting dataframe
        plotting_data[col_color_key] = adata.obs.iloc[plotting_data.Cell.values, adata.obs.columns.get_loc(col_color_key)].values
        plotting_data = plotting_data.merge(pd.get_dummies(plotting_data[col_color_key]), left_index=True, right_index=True)

        plotting_data = plotting_data.groupby(dial_key).mean().reset_index()
        plotting_data = plotting_data.sort_values(dial_key)

        # Standardize categorical and assemble values
        col_color_values = []
        for i, col_color_category in enumerate(col_color_categories):

            if col_color_category in plotting_data.columns:
                plotting_data[col_color_category] = standardize_(plotting_data[col_color_category])

                if n_convolve is not None:
                    plotting_data[col_color_category] = convolve(plotting_data[col_color_category].values, n_convolve)

                try:
                    if type(col_colors) == list:
                        col_color_values.append((plotting_data[col_color_category].values, col_colors[i]))
                    elif type(col_colors) == dict:
                        col_color_values.append((plotting_data[col_color_category].values, col_colors[col_color_category]))
                except:
                    raise ValueError('Not enough cmaps in col_colors')

    # Assemble expression values for the radial heatmap       
    var_values = [(plotting_data[var_name].values, color_map) for var_name in var_names_]

    # Assemble all values
    if col_color_key:
        values = col_color_values + ['skip'] + var_values
    else:
        values = var_values
    
    # Make table text
    table_text = []
    for i, value in enumerate(values):
        if value=='skip':
            table_text.append(['','',''])
        else:
            table_text.append([i, value[0], value[1]])
    table_text = pd.DataFrame(table_text, columns=['Ring', 'Marker', 'Color Map'])

    # Make radial heatmap
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax = make_radial_heatmap(ax, values, bins, len(values), 8, 10)
    
    if table_ax:
        pd.plotting.table(ax=table_ax, data=table_text, colLabels=table_text.columns, colWidths=[1/4]*3, fontsize=30)

    return ax