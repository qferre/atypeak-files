import pybedtools
import pandas as pd
import numpy as np
import sys

from plotnine import ggplot, aes, labs, scale_color_gradient, geom_point, geom_abline, scale_x_log10, scale_y_log10, geom_violin, geom_boxplot, position_dodge, scale_y_log10

"""
Preparation : to run this code, the appropriate files must be placed in the same directory as this script.
They have not been given here to save space.
    - Copy the atypeak result files from /bed/ to /script/atypeak_result/
    - Copy the selected CRM from the main atypeak Git (data/input/crm_id_coord.bed from <https://github.com/qferre/atypeak>) to /script/crm_selection.bed
    - Download the ReMap 2020 files by running the commands below
"""

# Process both cell lines, one at a time
CELL_LINES = ["mcf7", "jurkat"]
for CELL_LINE in CELL_LINES:

    ## Source files
    # Give the paths to atypeak result files. Move them to the appropriate positions.
    atypeak_result_file = "./atypeak_result/"+CELL_LINE+".bed"  # Results are given compressed here in /bed/*.xz
    crm_file = "./crm_selection.bed"                            # The usual CRM file given as input -- found in the atypeak Git directory
    remap_2020_file = "./remap2020_full_"+CELL_LINE+".bed"      # ReMap 2020 too large to be provided here, run the commands below
    # To download ReMAP 2020 files, run these commands:
    # wget http://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/remap2020_all_macs2_hg38_v1_0.bed.gz
    # grep -i {CELL_LINE} remap2020_all_macs2_hg38_v1_0.bed > remap2020_full_jurkat.bed # -i means case insensitive

    FIGURE_DIRECTORY = "./results/"+CELL_LINE+"/"

    print("Beginning...")

    ## ------------------------ Grouping the peaks by CRM ------------------- ##
    crm = pybedtools.BedTool(crm_file)
    total_crm = len(crm)

    ## Group atypeak results by CRM.
    # Use the result of `bedtools intersect -a file.bed -b crm_selected.bed -wa -wb`
    output_bed = pybedtools.BedTool(atypeak_result_file)
    total_peaks = len(output_bed)

    intersected = output_bed.intersect(crm, wa=True, wb=True)
    df = intersected.to_dataframe(names = ['chr','start','end','name','anomaly_score','strand','CRM_chr','CRM_start','CRM_end','CRM_name','CRM_score','CRM_strand','CRM_center','CRM_center+1'])
    df['full_crm_coordinate'] = df["CRM_chr"].map(str)+'.'+df["CRM_start"].map(str)+'.'+df["CRM_end"].map(str)
    df_group = df.groupby(['full_crm_coordinate'])

    # Same for ReMap 2020
    new_remap_bed = pybedtools.BedTool(remap_2020_file)
    intersected_new = new_remap_bed.intersect(crm, wa=True, wb=True)
    df_new = intersected_new.to_dataframe(names = ['chr','start','end','name','score','strand','top_start','top_end','color','CRM_chr','CRM_start','CRM_end','CRM_name','CRM_score','CRM_strand','CRM_center','CRM_center+1'])
    df_new['full_crm_coordinate'] = df_new["CRM_chr"].map(str)+'.'+df_new["CRM_start"].map(str)+'.'+df_new["CRM_end"].map(str)

    df_new_group = df_new.groupby(['full_crm_coordinate'])

    print("Intersections complete.")

    ## ---- Figure production ---- #

    res = pd.DataFrame({'crm_id':[],'nb_peaks_2018':[],'nb_peaks_2020':[],'average_atypeak_score':[], 'max_atypeak_score':[]})

    res_chunks = []

    ## Compare : number of peaks in 2018, in 2020, and atypeak mean score

    ## First integrate atypeak
    i = 0
    for name_crm, group in df_group:

        # Get the corresponding 2020 group
        group_new = df_new_group.get_group(name_crm)

        scores = [float(s) for s in group['anomaly_score']]
        npe = len(scores)
        asn = np.mean(scores)
        msn = np.max(scores)

        new_row = pd.DataFrame({'crm_id':[name_crm],'nb_peaks_2018': [npe], 'average_atypeak_score': [asn], 'max_atypeak_score':[msn] ,
        'nb_peaks_2020': [len(group_new)]})

        res_chunks.append(new_row)

        i += 1

        sys.stdout.write("\r" +"Processed CRMs : " + str(i) + ' / ' + str(total_crm))
        sys.stdout.flush()
   
    print("---")

    res = pd.concat(res_chunks, axis = 0)

    # Now that all CRMs are known, use them as ID for 2020 appending
    res = res.set_index('crm_id')
    res['update_ratio'] = res['nb_peaks_2020'] / res['nb_peaks_2018']


    # -------------- Working by TF

    # Add a TF column to all peaks in df and df_new.
    # Done before the grouping so it can still be grouped by CRM
    tf = [name.split('.')[1].lower() for name in df['name']]
    df['tf'] = tf

    tf_new = [name.split('.')[1].lower() for name in df_new['name']]
    df_new['tf'] = tf_new

    # Now some statistics for each peak
    peakstats = pd.DataFrame({
        'peak_score':[],
        'tf':[],
        'crm_id':[],
        'mean_crm_atypeak_score':[],
        'nb_peaks_same_tf_in_crm_2018':[],
        'nb_peaks_same_tf_in_crm_2020':[],
    })


    peakstats_chunks = []
    i = 0

    for name_crm, peaks_in_this_crm in df_group:

        # Get the corresponding 2020 group
        group_new = df_new_group.get_group(name_crm)

        for index, peak in peaks_in_this_crm.iterrows():

            this_peak_tf = peak['tf']
            same_tf_peaks_2018 = peaks_in_this_crm.loc[peaks_in_this_crm['tf'] == this_peak_tf] 
            same_tf_peaks_2020 = group_new.loc[group_new['tf'] == this_peak_tf] 

            new_row = pd.DataFrame({
                'peak_score':[peak['anomaly_score']],
                'tf':[peak['tf']],
                'crm_id':[name_crm],
                'mean_crm_atypeak_score':[np.mean(peaks_in_this_crm['anomaly_score'])],
                'nb_peaks_same_tf_in_crm_2018':[len(same_tf_peaks_2018)],
                'nb_peaks_same_tf_in_crm_2020':[len(same_tf_peaks_2020)],
                'total_bp_same_tf_in_crm_2018':[
                    np.sum([p['end'] - p['start'] for _, p in same_tf_peaks_2018.iterrows()])
                ],
                'total_bp_same_tf_in_crm_2020':[
                    np.sum([pnew['end'] - pnew['start'] for _, pnew in same_tf_peaks_2020.iterrows()])
                ]
            })

            peakstats_chunks.append(new_row)
            i += 1

        sys.stdout.write("\r" +"Processed peaks : " + str(i) + ' / ' + str(total_peaks))
        sys.stdout.flush()

    peakstats = pd.concat(peakstats_chunks, axis = 0)

    # Since I don't merge, I can do the basepairs update ratio
    peakstats['update_ratio_same_tf'] = peakstats['total_bp_same_tf_in_crm_2020'] / peakstats['total_bp_same_tf_in_crm_2018']
    peakstats['update_ratio_same_tf_peaks_nb'] = peakstats['nb_peaks_same_tf_in_crm_2020'] / peakstats['nb_peaks_same_tf_in_crm_2018']


    # Save
    res.to_csv(FIGURE_DIRECTORY+"crm_res.tsv", sep = '\t')
    peakstats.to_csv(FIGURE_DIRECTORY+"peakstats.tsv", sep = '\t')

    """
    # To reload :
    res = pd.read_csv(FIGURE_DIRECTORY+"crm_res.tsv", sep = '\t')
    peakstats = pd.read_csv(FIGURE_DIRECTORY+"peakstats.tsv", sep = '\t')
    """




    ## --------- For the figures

    p = (ggplot(data = res[0:10000], mapping = aes(x='nb_peaks_2020', y='nb_peaks_2018')) 
        + geom_point(mapping = aes(color='average_atypeak_score')) + scale_x_log10() + scale_y_log10()
        + labs(x="Nb. peaks Remap 2020", y="Nb. peaks Remap 2018", color = "Mean atyPeak score per CRM") + scale_color_gradient(low="red", high="blue"))
    p.save(FIGURE_DIRECTORY + "crm_nb_peaks_update.pdf", verbose = False)


    p = (ggplot(data = res[10000:13000], mapping = aes(x='nb_peaks_2018', y='update_ratio')) 
        + geom_point(mapping = aes(color='average_atypeak_score')) + scale_x_log10() + scale_y_log10()
        + labs(x="Nb. peaks Remap 2018", y="Nb peaks ReMap 2020/2018", color = "Mean atyPeak score per CRM") + scale_color_gradient(low="red", high="blue"))
    p.save(FIGURE_DIRECTORY + "crm_update_ratio.pdf", verbose = False)



    # Update ratio for the same TF (presumably, false positives should not be confirmed by more experiments)

    p3 = (ggplot(data = peakstats[10000:20000], mapping = aes(x='mean_crm_atypeak_score', y='total_bp_same_tf_in_crm_2020')) + scale_y_log10()
            + geom_point(mapping = aes(color='peak_score')) + scale_color_gradient(low="orange", high="blue")
            + labs(x="Mean CRM atypeak score", y="Total bp. in CRM same TF 2020", color="Peak score"))
    p3.save(FIGURE_DIRECTORY + "peak_confirmation_bp_2020.pdf", verbose = False)

    p3 = (ggplot(data = peakstats[10000:13000], mapping = aes(x='mean_crm_atypeak_score',y='update_ratio_same_tf')) + scale_y_log10()
            + geom_point(mapping = aes(color='peak_score')) + scale_color_gradient(low="orange", high="blue")
            + labs(x="Mean CRM atypeak score", y="Total bp. in CRM same TF 2020/2018", color="Peak score"))
    p3.save(FIGURE_DIRECTORY + "peak_confirmation_bp_update_ratio.pdf", verbose = False)


    # Using only crms that are relatively well characterized (above 500) to avoid low-crm biais
    sub = peakstats.loc[(peakstats['mean_crm_atypeak_score'] >= 500) & (peakstats['mean_crm_atypeak_score'] <= 1000)]

    p3 = (ggplot(data = sub[10000:15000], mapping = aes(x='peak_score', y='update_ratio_same_tf_peaks_nb'))
    + geom_point() + scale_y_log10() + labs(x="Peak score", y="Nb. peaks in CRM same TF 2020/2018", title = "CRMs with mean score >=500"))
    p3.save(FIGURE_DIRECTORY + "peak_confirmation_nb_update_ratio_well_characterized_crm.pdf", verbose = False)


    # Do it as a violin plot
    update_ratio_binarized = []
    for index, row in sub.iterrows():
        val = sub.loc[index, 'update_ratio_same_tf_peaks_nb']

        # Different thresholds for the two cell lines

        if CELL_LINE == 'jurkat':
            if val > 0 :   toadd = 'a__0-0.5'
            if val > 0.5 : toadd = 'b__0.5-1'
            if val > 1 :   toadd = 'c__1-2'
            if val > 2 :   toadd = 'd__2-3'
            if val > 3 :   toadd = 'e__3+'

        if CELL_LINE == "mcf7":
            if val > 0 :  toadd = 'a__0-3'
            if val > 3 :  toadd = 'b__3-10'
            if val > 5 :  toadd = 'c__5-10'
            if val > 10 : toadd = 'd__10+'


        update_ratio_binarized += [toadd]

    sub['update_ratio_bin'] = update_ratio_binarized
    sub['sqrt_peak_score'] = np.sqrt(sub['peak_score'])

    
    # Now do violin plot
    p4 = (ggplot(data=sub[0:10000], mapping = aes(x='update_ratio_bin', y='peak_score')) 
        + geom_violin(position = position_dodge(1), width = 1)
        + scale_y_log10()
        + geom_boxplot(position = position_dodge(1), width = 0.25))

    p4.save(FIGURE_DIRECTORY + "peak_confirmation_nb_update_ratio_well_characterized_crm_violin_plot.pdf", verbose = False)









