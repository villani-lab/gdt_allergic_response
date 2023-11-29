import pegasus as pg
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.colors as clr


# colormap for hexplots
blues_cmap = clr.LinearSegmentedColormap.from_list('gene_cmap', ["#e0e0e1", '#4576b8', '#02024a'], N=200)


def plot_umap(df, title, colors, cluster_col, show=False,
              fig_size=(4, 4), wspace=0.4, marker_multiplier=20,
              cluster_font_size=10, axis_font_size=14, ncol_legend=1,
              dpi=300, simple=False, save_fig=False, save_name=None):

    def get_umap_numbers(df, cluster_name_col):
        # get cluster names ranked by abundance
        n_cells = df.obs[cluster_name_col].astype('category').value_counts()

        # order cluster names by abundance, doublets last
        names = n_cells.index
        non_doublets = names[['doublets' not in x.lower() for x in names]]
        doublets = names[['doublets' in x.lower() for x in names]]
        names = list(non_doublets)
        names.extend(doublets)

        # assign umap numbers
        number_dict = dict(zip(names, [*range(1, len(names) + 1)]))
        df.obs['umap_numbers'] = [number_dict[x] for x in df.obs[cluster_name_col]]
        df.obs['umap_numbers'] = df.obs['umap_numbers'].astype('category')

        # get number, name, and n_cells for legend
        legend_dict = {}
        for n, name in enumerate(names, start=1):
            legend_dict[n] = {'name': name, 'n_cells': n_cells[name]}
        return df, legend_dict

    # number clusters by number of cells (with doublets last)
    df, legend_dict = get_umap_numbers(df, cluster_col)

    # get legend labels
    labels = []
    for n in legend_dict.keys():
        labels.append(f"{n}. {legend_dict[n]['name']} (n={legend_dict[n]['n_cells']:,})")

    # create base_figure
    loc = 'right margin' if simple else 'on data'
    palette = ','.join(colors)
    fig = pg.scatter(df,
                     basis="umap",
                     attrs="umap_numbers",
                     legend_loc=loc,
                     palettes=palette,
                     panel_size=fig_size,
                     wspace=wspace,
                     return_fig=True)
    fig.patch.set_visible(True)
    fig.text(x=.07, y=.09, s=f'{df.shape[0]:,} cells', fontsize=cluster_font_size)
    ax = fig.axes[0]
    # edit cluster number font
    ax_is_text = [isinstance(i, matplotlib.text.Text) for i in ax.get_children()]
    ax_text_idx = np.where(ax_is_text)[0]
    for idx in ax_text_idx:
        if ax.get_children()[idx]._text in df.obs['umap_numbers'].cat.categories.astype('string'):
            text = ax.get_children()[idx]
            text.set_path_effects([path_effects.Stroke(linewidth=4, foreground='white'),
                                   path_effects.Normal()])
            text.set_fontsize(cluster_font_size)
    # add legend
    ax.legend(labels,
              bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0,
              markerscale=marker_multiplier, frameon=False,
              prop={'size': cluster_font_size},
              ncol=ncol_legend)
    # edit font sizes
    ax.xaxis.label.set_fontsize(axis_font_size)
    ax.yaxis.label.set_fontsize(axis_font_size)
    ax.set_title(title, fontsize=axis_font_size)
    plt.tight_layout()
    if save_fig:
        plt.rcParams['svg.fonttype'] = 'none'
        plt.savefig(save_name, dpi=dpi)
    if show:
        plt.show()
    plt.close()


def make_dotplot(df, cluster_order, gene_order, title,
                 obs_name_col='cluster_names', figsize=(8, 4), save=False, save_name=None):
    # remove doublets if any
    df = df[['doublets' not in x.lower() for x in df.obs[obs_name_col]],].copy()
    df.obs[obs_name_col] = df.obs[obs_name_col].cat.remove_unused_categories()

    # set figure font
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    plt.rcParams.update({'font.size': 19})

    # create figure
    fig, ax = plt.subplots(figsize=(12, 6))
    dotplot = sc.pl.dotplot(df, gene_order, obs_name_col,
                            var_group_rotation=30, standard_scale='var', ax=ax,
                            categories_order=cluster_order,
                            figsize=figsize,
                            return_fig=True)
    dotplot.legend(width=2, colorbar_title='Scaled\nExpression')
    dotplot.make_figure()
    axes_dict = dotplot.get_axes()
    axes_dict["mainplot_ax"].set_axisbelow(True)
    axes_dict["mainplot_ax"].xaxis.tick_top()
    axes_dict["mainplot_ax"].set_title(f'{title}: Cluster Marker Genes')
    fig.tight_layout()
    if save:
        plt.savefig(save_name, dpi=300)
    plt.show()
    plt.close()


def signature_score_per_cell(data, gene_set, score_title):
    # Get rid of genes that aren't in data
    orig_len = len(gene_set)
    gene_set = [gene for gene in gene_set if gene in data.var_names]
    print(str(len(gene_set)) + "/" + str(orig_len) + " of the gene set genes are measured")

    # Limit the data to just those genes
    dat = data[:, gene_set].X
    dat = dat.toarray()
    mean = dat.mean(axis=0)
    var = dat.var(axis=0)
    std = np.sqrt(var)

    with np.errstate(divide="ignore", invalid="ignore"):
        dat = (dat - mean) / std
    dat[dat < -5] = -5
    dat[dat > 5] = 5

    scores = dat.mean(axis=1)
    data.obs[score_title] = scores


def hex_plot(df, title, col, gridsize=200, cmap='YlOrRd',
             save=False, save_name=None, dpi=300):
    # get umap coordinates
    umap_coords = pd.DataFrame(df.obsm['X_umap'], columns=['x', 'y'])
    x = umap_coords['x']
    y = umap_coords['y']

    # get data to plot
    data = df.obs[col]

    # create plot
    fig, ax = plt.subplots(figsize=(6, 4.5))
    hb = ax.hexbin(x=x,
                   y=y,
                   C=data,
                   cmap=cmap,
                   gridsize=gridsize,
                   edgecolors='none')
    cb = fig.colorbar(hb, ax=ax, shrink=.75, aspect=10)
    cb.ax.set_title(col, loc='left', fontsize=14)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)
    ax.set_title(title, fontsize=18)
    ax.tick_params(left=False, labelleft=False,
                   bottom=False, labelbottom=False)

    if save:
        plt.savefig(f'{save_name}_hexbin{gridsize}.pdf', dpi=dpi)
    plt.show()
    plt.close()