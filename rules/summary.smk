####################################
# Read summary statistics          #
####################################

rule summary_stats:
    input: 
        all_summary=READS_IMPORT_SUMMARY
    output:
        plot=SUMMARY_PLOT,
    run:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        df = pd.read_csv(input.all_summary, quoting=3, sep='\t')
        
        # Sometimes these summary files have quotation marks around each row entry
        if type(df.iloc[0,0])=='str':
            if df.iloc[0,0].find('\"')>-1:
                df.columns    = df.columns.map(lambda x: x.strip('\"'))
                df.iloc[:,0]  = df.iloc[:,0].str.strip('\"')
                df.iloc[:,-1] = df.iloc[:,-1].str.strip('\"').astype('float')

        length = 'sequence_length_template'
        qscore = 'mean_qscore_template'
        if 'sequence_length_2d' in df.columns:
            # This is 1Dsq data, not 1D. Update column ID name for querying
            length = length.replace('template', '2d')
            qscore = qscore.replace('template', '2d')

        fig = plt.figure(figsize=(12,10))
        
        # Prevent outlier lengths from screwing up histogram
        plot_df = df.query('({l} <= 70000)'.format(l=length))

        ax1 = fig.add_subplot(2,1,1)
        ax1.hist(plot_df[qscore], bins=50, linewidth=0, alpha=0.5, color='b')
        ax1.set_xlabel('qscore')
        ax1.set_ylabel('count')
        ax1.set_title('%s reads; mean qscore = %.1f' % (len(df[qscore]), np.mean(df[qscore])))

        ax2 = fig.add_subplot(2,1,2)
        ax2.hist(plot_df[length], bins=50, linewidth=0, alpha=0.5, color='b')
        ax2.set_xlabel('length')
        ax2.set_ylabel('count')
        ax2.set_title('mean length = %i' % np.mean(df[length]))

        fig.savefig(output.plot)
