import os,sys
import pandas as pd
import argparse
import collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
import numpy as np

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('paf', help='Alignment file produced by all-vs-all minimap2', type=str)

    # Optional arguments
    parser.add_argument('-p', '--prefix', help='Output file prefix [output]', type=str, default='output')
    parser.add_argument('-s', '--min_score_frac', help='Minimum fraction of maximum score for a cluster [0.8]', type=float, default=0.8)
    parser.add_argument('-n', '--min_reads', help='Minimum reads to create cluster [5]', type=int, default=5)

    # Parse arguments
    args = parser.parse_args()

    return args

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the 'center' of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        # kwargs['color_threshold'] = max_d
        kwargs['color_threshold'] = None
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = sch.dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        # plt.title('Hierarchical Clustering Dendrogram (truncated)')
        # plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        plt.xticks([])
        # plt.yticks([])
        
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                # plt.annotate('%.3g' % y, (x, y), xytext=(0, -5),
                #              textcoords='offset points',
                #              va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    ax = plt.gca()
    # ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    return ddata

def run_clustering(df, args, default_sim):
    Z   = sch.linkage(df, 'ward', optimal_ordering=True)

    c, coph_dists = sch.cophenet(Z, pdist(df))

    args.max_dist = 3*np.median(Z[:,2])
    
    fig = plt.figure(figsize=(12, 12))

    buff          = 0.005
    x_start       = 0.03
    y_start       = 0.23
    cb_width      = 0.05
    cb_lab_buff   = 0.03
    ld_xstart     = x_start + cb_width + cb_lab_buff + buff
    ld_width      = 0.05
    mat_xstart    = ld_xstart + ld_width + buff
    mat_width     = 0.6
    mat_height    = 0.6
    ccolor_ystart = y_start + mat_height + buff
    ccolor_h      = 0.03
    td_ystart     = y_start + mat_height + buff + ccolor_h + buff
    td_height     = 0.12

    # left side dendrogram
    ax1 = fig.add_axes([ld_xstart, y_start, ld_width, mat_height])
    ax1.set_xticks([])
    ax1.set_yticks([])
    Z1 = fancy_dendrogram(
                            Z,
                            show_contracted=True,
                            leaf_rotation=90.,  # rotates the x axis labels
                            leaf_font_size=8.,  # font size for the x axis labels
                            annotate_above=0.07,
                            orientation="left"
                        )

    # top side dendogram
    ax2 = fig.add_axes([mat_xstart, td_ystart, mat_width, td_height])
    ax2.set_xticks([])
    Z2 = fancy_dendrogram(
                            Z,
                            show_contracted=True,
                            leaf_rotation=90.,  # rotates the x axis labels
                            leaf_font_size=8.,  # font size for the x axis labels
                            annotate_above=0.07,
                            # max_d=max_d,
                        )

    # main heat-map
    orig_cmap    = matplotlib.cm.Greys
    shifted_cmap = orig_cmap
    # shifted_cmap = shiftedColorMap(orig_cmap, start=-0.5, midpoint=0.750, name='shifted')
    axmatrix     = fig.add_axes([mat_xstart, y_start, mat_width, mat_height])
    idx1         = Z1['leaves']
    idx2         = Z2['leaves']
    mat_result   = df.iloc[idx1, :]
    mat_result   = mat_result.iloc[:, idx2]
    im           = axmatrix.matshow(mat_result, aspect='auto', origin='lower', cmap=shifted_cmap, vmin=0)

    axmatrix.set_xticks(range(len(idx1)))
    axmatrix.set_xticklabels(df.index[idx1], minor=False, fontsize=4)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()

    plt.xticks(rotation=-90, fontsize=4)    

    axmatrix.set_yticks(range(len(idx2)))
    axmatrix.set_yticklabels(df.index[idx2], minor=False, fontsize=4)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()

    # colorbar
    cbaxes = fig.add_axes([x_start, y_start, cb_width, mat_height]) 
    cbar   = fig.colorbar(im, cax = cbaxes, ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

    # top side color coding for cluster membership
    [axc_x, axc_y, axc_w, axc_h] = [mat_xstart, ccolor_ystart, mat_width, ccolor_h]
    axc      = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar

    clusters = sch.fcluster(Z, args.max_dist, criterion='distance')
    # clusters = sch.fcluster(Z, 4, depth=depth)
    
    df.loc[:,'cluster'] = clusters
    df                  = df.reset_index()

    cluster_dfs         = df.groupby('cluster')
    
    return clusters,cluster_dfs,cbar,axc,idx2

def create_figure(clusters, idx2, axc, args):
    dc       = np.array(clusters[idx2], dtype=int)
    dc.shape = (1,len(clusters))

    vals = np.linspace(0,1,256)
    np.random.shuffle(vals)
    cmap_c = plt.cm.colors.ListedColormap(plt.cm.hsv(vals))

    im_c     = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
    def find_mid_index(vec, value):
        idx_vec = np.where(vec==value)[0]
        return idx_vec[int(len(idx_vec)/2)]
    labels = dict([(value, find_mid_index(dc[0], value)) for value in set(clusters)])
    for l,x in labels.items():
        axc.text(x, 0, l, horizontalalignment='center', fontsize=12)
    axc.set_xticks([])
    axc.set_yticks([])

    plt.savefig('{}.heatmap.png'.format(args.prefix), dpi=250)

def main(args):
    cols = ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', \
            'tlen', 'tstart', 'tend', 'matches', 'alnlen', 'mapqv', \
            'tp', 'cm', 's1', 'dv', 'rl']
    df = pd.read_csv(args.paf, sep='\t', header=None, usecols=list(range(17)), names=cols)
    
    bin_id = int(args.paf.split('/')[-2])

    if df.shape[0]>0:
        df['longread'] = df[['qlen','tlen']].apply(max, axis=1)

        df['dv'] = df['dv'].apply(lambda x: float(x.split(':')[-1]))
        df['sim'] = 1 - df['dv']
    
        df1 = df.groupby(['qname', 'tname']).apply(lambda x: x.sort_values(['alnlen'], ascending=False)).reset_index(drop=True)
        df1 = df1.groupby(['qname', 'tname']).head(1)

        # Combine top alignment length with its similarity score to get distance metric
        df1['score'] = (df1['alnlen']/10000)**2 * df1['sim']

        df1 = df1.loc[:,['qname','tname','score','qlen','tlen']]
        df2 = df1.rename(columns={'qname':'tname', 'tname':'qname', 'qlen':'tlen', 'tlen':'qlen'})
        df_both = df1.append(df2, ignore_index=True, sort=False)

        nan_value = 0
        df_pivot = df_both.pivot_table(index='qname', columns='tname', values='score').fillna(nan_value)

        clusters,cluster_dfs,cbar,axc,idx2 = run_clustering(df_pivot, args, nan_value)
        
        create_figure(clusters, idx2, axc, args)
        
        rls = dict([(row['qname'],row['qlen']) for idx,row in df_both.loc[:,['qname','qlen']].iterrows()])

        clust_info_out = '{}.info.csv'.format(args.prefix)
        with open(clust_info_out, 'w') as f:
            f.write('bin_id,cluster,read_id,read_len,clust_read_score,frac_max_score\n')
            for cluster,df_c in cluster_dfs:
                keep_cols = df_c['qname'].append(pd.Series(['qname', 'cluster']))
                df_c = df_c.loc[:,keep_cols] # only keep the columns of reads from the cluster
                rls_c = []
                read_means = []
                for idx,row in df_c.iterrows():
                    nonself   = row.replace(1.00, np.nan).dropna()
                    dist_idx   = [i for i in nonself.index if (i!='qname' and i!='cluster')]
                    read_means.append(nonself.loc[dist_idx].mean())
                    rls_c.append(rls[row['qname']])

                df_c.loc[:,'score_mean'] = read_means
                clust_mean_score = df_c.loc[:,'score_mean'].mean()
                clust_mean_rl = np.mean(rls_c)
                max_score = (clust_mean_rl/10000)**2
                max_score_frac = clust_mean_score / max_score
                # print(bin_id, row['cluster'], df_c.shape[0], clust_mean_rl, clust_mean_score, max_score_frac)
                df_c = df_c.set_index('qname')
                if df_c.shape[0]>=args.min_reads and max_score_frac>=args.min_score_frac:
                    for r in df_c.index.values:
                        f.write('{},{},{},{},{},{}\n'.format(bin_id, \
                                                             cluster, \
                                                             r, \
                                                             rls[r], \
                                                             round(df_c.loc[r,'score_mean'],3), \
                                                             round(max_score_frac,3)))

        
    else:
        # There are no alignments to parse
        # print('keep\tread\tall_v_all_aligns')
        # print('{}\t{}\t{}'.format('False', 'none', 0))
        pass

if __name__ == '__main__':
    args = parse_args()

    main(args)