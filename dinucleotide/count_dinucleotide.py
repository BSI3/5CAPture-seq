import pandas as pd
import scipy.stats as sp
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def count_dinucleotide(file):
    base = ['A', 'T', 'C', 'G']
    dinucleotide = {}
    for i in base:
        for j in base:
            di = i + j
            dinucleotide.update({di: {'-4/-3': 0, '-3/-2': 0, '-2/-1': 0, '-1/1': 0, '1/2': 0, '2/3': 0, '3/4': 0}})
    # print len(dinucleotide)

    for i in open(file):
        if i[0] == '>':
            # print i
            continue
        else:
            # j = i[46:54]
            j = i[0:8]
            # print j
            m = -5
            n = -4
            for a in j:
                if m == -5 and n == -4:
                    b = a
                    m += 1
                    n += 1
                else:
                    # print b+a
                    try:
                        dinucleotide[b + a][str(m) + '/' + str(n)] += 1
                    except:
                        print('cannot find ' + (b + a) + ' in ' + (str(m) + '/' + str(n)))
                    else:
                        b = a
                        m += 1
                        n += 1
                        if m == 0:
                            m += 1
                        if n == 0:
                            n += 1

    # print dinucleotide
    df = pd.DataFrame.from_dict(dinucleotide, orient='index')
    df = df[['-4/-3', '-3/-2', '-2/-1', '-1/1', '1/2', '2/3', '3/4']]
    df.loc['total', :] = df.sum()
    # print df
    return df


def fisher_exact(df, df_random):
    df_t = df.T
    df_t = pd.DataFrame(df_t,
                        columns=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC',
                                 'TG', 'TT', 'total'],
                        index=['-4/-3', '-3/-2', '-2/-1', '-1/1', '1/2', '2/3', '3/4'])
    # print df_t
    df_random_t = df_random.T
    df_random_t = pd.DataFrame(df_random_t,
                               columns=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA',
                                        'TC', 'TG', 'TT', 'total'],
                               index=['-4/-3', '-3/-2', '-2/-1', '-1/1', '1/2', '2/3', '3/4'])

    # oddsratio1, pvalue1 = sp.fisher_exact([[0,9],[10,39]])
    # oddsratio, pvalue = sp.fisher_exact([[df1_t.loc['-2/-1','aa'],df1_t.loc['-2/-1','total']],[df2_t.loc['-2/-1','aa'],df2_t.loc['-2/-1','total']]])
    # print df1_t
    # print pvalue1
    # print pvalue
    # print '\n\n'
    df_pvalue = pd.DataFrame(columns=['-4/-3', '-3/-2', '-2/-1', '-1/1', '1/2', '2/3', '3/4'],
                             index=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC',
                                    'TG', 'TT'])  # calculate left tailed,less
    df_pvalue_log = pd.DataFrame(columns=['-4/-3', '-3/-2', '-2/-1', '-1/1', '1/2', '2/3', '3/4'],
                                 index=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA',
                                        'TC', 'TG', 'TT'])
    df_pvalue_gr = pd.DataFrame(columns=['-4/-3', '-3/-2', '-2/-1', '-1/1', '1/2', '2/3', '3/4'],
                                index=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA',
                                       'TC', 'TG', 'TT'])  # calculate right tailed,greater
    df_pvalue_gr_log = pd.DataFrame(columns=['-4/-3', '-3/-2', '-2/-1', '-1/1', '1/2', '2/3', '3/4'],
                                    index=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA',
                                           'TC', 'TG', 'TT'])

    for (index1, row1), (index2, row2) in zip(df_t.iterrows(), df_random_t.iterrows()):
        ######
        ####calculating right tailed, greater

        oddsratio_aa, pvalue_aa = sp.fisher_exact([[row1['AA'], (row1['total'] - row1['AA'])],
                                                   [row2['AA'], (row2['total'] - row2['AA'])]],
                                                  alternative='less')
        oddsratio_ac, pvalue_ac = sp.fisher_exact([[row1['AC'], (row1['total'] - row1['AC'])],
                                                   [row2['AC'], (row2['total'] - row2['AC'])]],
                                                  alternative='less')
        oddsratio_ag, pvalue_ag = sp.fisher_exact([[row1['AG'], (row1['total'] - row1['AG'])],
                                                   [row2['AG'], (row2['total'] - row2['AG'])]],
                                                  alternative='less')
        oddsratio_at, pvalue_at = sp.fisher_exact([[row1['AT'], (row1['total'] - row1['AT'])],
                                                   [row2['AT'], (row2['total'] - row2['AT'])]],
                                                  alternative='less')
        df_pvalue.loc['AA'][index1] = pvalue_aa
        df_pvalue.loc['AC'][index1] = pvalue_ac
        df_pvalue.loc['AG'][index1] = pvalue_ag
        df_pvalue.loc['AT'][index1] = pvalue_at

        df_pvalue_log.loc['AA'][index1] = -10 * np.log10(pvalue_aa)
        df_pvalue_log.loc['AC'][index1] = -10 * np.log10(pvalue_ac)
        df_pvalue_log.loc['AG'][index1] = -10 * np.log10(pvalue_ag)
        df_pvalue_log.loc['AT'][index1] = -10 * np.log10(pvalue_at)

        oddsratio_ca, pvalue_ca = sp.fisher_exact([[row1['CA'], (row1['total'] - row1['CA'])],
                                                   [row2['CA'], (row2['total'] - row2['CA'])]],
                                                  alternative='less')
        oddsratio_cc, pvalue_cc = sp.fisher_exact([[row1['CC'], (row1['total'] - row1['CC'])],
                                                   [row2['CC'], (row2['total'] - row2['CC'])]],
                                                  alternative='less')
        oddsratio_cg, pvalue_cg = sp.fisher_exact([[row1['CG'], (row1['total'] - row1['CG'])],
                                                   [row2['CG'], (row2['total'] - row2['CG'])]],
                                                  alternative='less')
        oddsratio_ct, pvalue_ct = sp.fisher_exact([[row1['CT'], (row1['total'] - row1['CT'])],
                                                   [row2['CT'], (row2['total'] - row2['CT'])]],
                                                  alternative='less')

        df_pvalue.loc['CA'][index1] = pvalue_ca
        df_pvalue.loc['CC'][index1] = pvalue_cc
        df_pvalue.loc['CG'][index1] = pvalue_cg
        df_pvalue.loc['CT'][index1] = pvalue_ct

        df_pvalue_log.loc['CA'][index1] = -10 * np.log10(pvalue_ca)
        df_pvalue_log.loc['CC'][index1] = -10 * np.log10(pvalue_cc)
        df_pvalue_log.loc['CG'][index1] = -10 * np.log10(pvalue_cg)
        df_pvalue_log.loc['CT'][index1] = -10 * np.log10(pvalue_ct)

        oddsratio_ga, pvalue_ga = sp.fisher_exact([[row1['GA'], (row1['total'] - row1['GA'])],
                                                   [row2['GA'], (row2['total'] - row2['GA'])]],
                                                  alternative='less')
        oddsratio_gc, pvalue_gc = sp.fisher_exact([[row1['GC'], (row1['total'] - row1['GC'])],
                                                   [row2['GC'], (row2['total'] - row2['GC'])]],
                                                  alternative='less')
        oddsratio_gg, pvalue_gg = sp.fisher_exact([[row1['GG'], (row1['total'] - row1['GG'])],
                                                   [row2['GG'], (row2['total'] - row2['GG'])]],
                                                  alternative='less')
        oddsratio_gt, pvalue_gt = sp.fisher_exact([[row1['GT'], (row1['total'] - row1['GT'])],
                                                   [row2['GT'], (row2['total'] - row2['GT'])]],
                                                  alternative='less')

        df_pvalue.loc['GA'][index1] = pvalue_ga
        df_pvalue.loc['GC'][index1] = pvalue_gc
        df_pvalue.loc['GG'][index1] = pvalue_gg
        df_pvalue.loc['GT'][index1] = pvalue_gt

        df_pvalue_log.loc['GA'][index1] = -10 * np.log10(pvalue_ga)
        df_pvalue_log.loc['GC'][index1] = -10 * np.log10(pvalue_gc)
        df_pvalue_log.loc['GG'][index1] = -10 * np.log10(pvalue_gg)
        df_pvalue_log.loc['GT'][index1] = -10 * np.log10(pvalue_gt)

        oddsratio_ta, pvalue_ta = sp.fisher_exact([[row1['TA'], (row1['total'] - row1['TA'])],
                                                   [row2['TA'], (row2['total'] - row2['TA'])]],
                                                  alternative='less')
        oddsratio_tc, pvalue_tc = sp.fisher_exact([[row1['TC'], (row1['total'] - row1['TC'])],
                                                   [row2['TC'], (row2['total'] - row2['TC'])]],
                                                  alternative='less')
        oddsratio_tg, pvalue_tg = sp.fisher_exact([[row1['TG'], (row1['total'] - row1['TG'])],
                                                   [row2['TG'], (row2['total'] - row2['TG'])]],
                                                  alternative='less')
        oddsratio_tt, pvalue_tt = sp.fisher_exact([[row1['TT'], (row1['total'] - row1['TT'])],
                                                   [row2['TT'], (row2['total'] - row2['TT'])]],
                                                  alternative='less')

        df_pvalue.loc['TA'][index1] = pvalue_ta
        df_pvalue.loc['TC'][index1] = pvalue_tc
        df_pvalue.loc['TG'][index1] = pvalue_tg
        df_pvalue.loc['TT'][index1] = pvalue_tt

        df_pvalue_log.loc['TA'][index1] = -10 * np.log10(pvalue_ta)
        df_pvalue_log.loc['TC'][index1] = -10 * np.log10(pvalue_tc)
        df_pvalue_log.loc['TG'][index1] = -10 * np.log10(pvalue_tg)
        df_pvalue_log.loc['TT'][index1] = -10 * np.log10(pvalue_tt)
        ######
        ##########################calculating right tailed, greater,gr
        ######

        oddsratio_aa_gr, pvalue_aa_gr = sp.fisher_exact([[row1['AA'], (row1['total'] - row1['AA'])],
                                                         [row2['AA'], (row2['total'] - row2['AA'])]],
                                                        alternative='greater')
        oddsratio_ac_gr, pvalue_ac_gr = sp.fisher_exact([[row1['AC'], (row1['total'] - row1['AC'])],
                                                         [row2['AC'], (row2['total'] - row2['AC'])]],
                                                        alternative='greater')
        oddsratio_ag_gr, pvalue_ag_gr = sp.fisher_exact([[row1['AG'], (row1['total'] - row1['AG'])],
                                                         [row2['AG'], (row2['total'] - row2['AG'])]],
                                                        alternative='greater')
        oddsratio_at_gr, pvalue_at_gr = sp.fisher_exact([[row1['AT'], (row1['total'] - row1['AT'])],
                                                         [row2['AT'], (row2['total'] - row2['AT'])]],
                                                        alternative='greater')

        df_pvalue_gr.loc['AA'][index1] = pvalue_aa_gr
        df_pvalue_gr.loc['AC'][index1] = pvalue_ac_gr
        df_pvalue_gr.loc['AG'][index1] = pvalue_ag_gr
        df_pvalue_gr.loc['AT'][index1] = pvalue_at_gr

        df_pvalue_gr_log.loc['AA'][index1] = -10 * np.log10(pvalue_aa_gr)
        df_pvalue_gr_log.loc['AC'][index1] = -10 * np.log10(pvalue_ac_gr)
        df_pvalue_gr_log.loc['AG'][index1] = -10 * np.log10(pvalue_ag_gr)
        df_pvalue_gr_log.loc['AT'][index1] = -10 * np.log10(pvalue_at_gr)

        oddsratio_ca_gr, pvalue_ca_gr = sp.fisher_exact([[row1['CA'], (row1['total'] - row1['CA'])],
                                                         [row2['CA'], (row2['total'] - row2['CA'])]],
                                                        alternative='greater')
        oddsratio_cc_gr, pvalue_cc_gr = sp.fisher_exact([[row1['CC'], (row1['total'] - row1['CC'])],
                                                         [row2['CC'], (row2['total'] - row2['CC'])]],
                                                        alternative='greater')
        oddsratio_cg_gr, pvalue_cg_gr = sp.fisher_exact([[row1['CG'], (row1['total'] - row1['CG'])],
                                                         [row2['CG'], (row2['total'] - row2['CG'])]],
                                                        alternative='greater')
        oddsratio_ct_gr, pvalue_ct_gr = sp.fisher_exact([[row1['CT'], (row1['total'] - row1['CT'])],
                                                         [row2['CT'], (row2['total'] - row2['CT'])]],
                                                        alternative='greater')

        df_pvalue_gr.loc['CA'][index1] = pvalue_ca_gr
        df_pvalue_gr.loc['CC'][index1] = pvalue_cc_gr
        df_pvalue_gr.loc['CG'][index1] = pvalue_cg_gr
        df_pvalue_gr.loc['CT'][index1] = pvalue_ct_gr

        df_pvalue_gr_log.loc['CA'][index1] = -10 * np.log10(pvalue_ca_gr)
        df_pvalue_gr_log.loc['CC'][index1] = -10 * np.log10(pvalue_cc_gr)
        df_pvalue_gr_log.loc['CG'][index1] = -10 * np.log10(pvalue_cg_gr)
        df_pvalue_gr_log.loc['CT'][index1] = -10 * np.log10(pvalue_ct_gr)

        oddsratio_ga_gr, pvalue_ga_gr = sp.fisher_exact([[row1['GA'], (row1['total'] - row1['GA'])],
                                                         [row2['GA'], (row2['total'] - row2['GA'])]],
                                                        alternative='greater')
        oddsratio_gc_gr, pvalue_gc_gr = sp.fisher_exact([[row1['GC'], (row1['total'] - row1['GC'])],
                                                         [row2['GC'], (row2['total'] - row2['GC'])]],
                                                        alternative='greater')
        oddsratio_gg_gr, pvalue_gg_gr = sp.fisher_exact([[row1['GG'], (row1['total'] - row1['GG'])],
                                                         [row2['GG'], (row2['total'] - row2['GG'])]],
                                                        alternative='greater')
        oddsratio_gt_gr, pvalue_gt_gr = sp.fisher_exact([[row1['GT'], (row1['total'] - row1['GT'])],
                                                         [row2['GT'], (row2['total'] - row2['GT'])]],
                                                        alternative='greater')

        df_pvalue_gr.loc['GA'][index1] = pvalue_ga_gr
        df_pvalue_gr.loc['GC'][index1] = pvalue_gc_gr
        df_pvalue_gr.loc['GG'][index1] = pvalue_gg_gr
        df_pvalue_gr.loc['GT'][index1] = pvalue_gt_gr

        df_pvalue_gr_log.loc['GA'][index1] = -10 * np.log10(pvalue_ga_gr)
        df_pvalue_gr_log.loc['GC'][index1] = -10 * np.log10(pvalue_gc_gr)
        df_pvalue_gr_log.loc['GG'][index1] = -10 * np.log10(pvalue_gg_gr)
        df_pvalue_gr_log.loc['GT'][index1] = -10 * np.log10(pvalue_gt_gr)

        oddsratio_ta_gr, pvalue_ta_gr = sp.fisher_exact([[row1['TA'], (row1['total'] - row1['TA'])],
                                                         [row2['TA'], (row2['total'] - row2['TA'])]],
                                                        alternative='greater')
        oddsratio_tc_gr, pvalue_tc_gr = sp.fisher_exact([[row1['TC'], (row1['total'] - row1['TC'])],
                                                         [row2['TC'], (row2['total'] - row2['TC'])]],
                                                        alternative='greater')
        oddsratio_tg_gr, pvalue_tg_gr = sp.fisher_exact([[row1['TG'], (row1['total'] - row1['TG'])],
                                                         [row2['TG'], (row2['total'] - row2['TG'])]],
                                                        alternative='greater')
        oddsratio_tt_gr, pvalue_tt_gr = sp.fisher_exact([[row1['TT'], (row1['total'] - row1['TT'])],
                                                         [row2['TT'], (row2['total'] - row2['TT'])]],
                                                        alternative='greater')

        df_pvalue_gr.loc['TA'][index1] = pvalue_ta_gr
        df_pvalue_gr.loc['TC'][index1] = pvalue_tc_gr
        df_pvalue_gr.loc['TG'][index1] = pvalue_tg_gr
        df_pvalue_gr.loc['TT'][index1] = pvalue_tt_gr

        df_pvalue_gr_log.loc['TA'][index1] = -10 * np.log10(pvalue_ta_gr)
        df_pvalue_gr_log.loc['TC'][index1] = -10 * np.log10(pvalue_tc_gr)
        df_pvalue_gr_log.loc['TG'][index1] = -10 * np.log10(pvalue_tg_gr)
        df_pvalue_gr_log.loc['TT'][index1] = -10 * np.log10(pvalue_tt_gr)

    # print df_pvalue
    # print '\n\n'
    # print df_pvalue_log
    # print '\n\n'
    # print df_pvalue_gr
    # print '\n\n'
    # print df_pvalue_gr_log
    # print '\n\n'
    return df_pvalue, df_pvalue_log, df_pvalue_gr, df_pvalue_gr_log


file_output = open(
    '/Users/proton/Desktop/count_dinucleotide_10k/new_clustering_results_editminus.txt',
    'w')
file_output.write('#################Random Sequences#################' + '\n')
#####

random = '/Users/proton/Desktop/random_position_10k/'
random1 = 'cec1_random.txt'
random4 = 'cec2_random.txt'

df_random1 = count_dinucleotide(random + random1)
file_output.write(random1 + '\n')
file_output.write(str(df_random1) + '\n\n')

df_random4 = count_dinucleotide(random + random4)
file_output.write(random4 + '\n')
file_output.write(str(df_random4) + '\n\n')
#####

dirname = '/Users/proton/Desktop/xt_88base_10k/'
file1 = 'cec1_xt8_editminus.txt'
file4 = 'cec2_xt8_editminus.txt'

df1 = count_dinucleotide(dirname + file1)  # cec1
df_pvalue11, df_pvalue12, df_pvalue13, df_pvalue14 = fisher_exact(df1, df_random1)

df4 = count_dinucleotide(dirname + file4)  # cec2
df_pvalue41, df_pvalue42, df_pvalue43, df_pvalue44 = fisher_exact(df4, df_random4)
#####

# df_pvaluex1 = unnormailzed pvalue -less
# df_pvaluex2 = normailzed-log pvalue -less >>plotted this
# df_pvaluex3 = unnormailzed pvalue -greater
# df_pvaluex4 = normailzed-log pvalue -greater >>plotted this

file_output.write('#################Interested Sequences#################' + '\n')
file_output.write(file1 + '\n')
file_output.write(str(df1) + '\n\n')
file_output.write('lesser' + '\n' + str(df_pvalue11) + '\n\n')
file_output.write('lesser-log' + '\n' + str(df_pvalue12) + '\n\n')
file_output.write('greater' + '\n' + str(df_pvalue13) + '\n\n')
file_output.write('greater-log' + '\n' + str(df_pvalue14) + '\n\n')

file_output.write(file4 + '\n')
file_output.write(str(df4) + '\n\n')
file_output.write('lesser' + '\n' + str(df_pvalue41) + '\n\n')
file_output.write('lesser-log' + '\n' + str(df_pvalue42) + '\n\n')
file_output.write('greater' + '\n' + str(df_pvalue43) + '\n\n')
file_output.write('greater-log' + '\n' + str(df_pvalue44) + '\n\n')

file_output.close()

df_pvalue12 = df_pvalue12.astype(float)
df_pvalue42 = df_pvalue42.astype(float)

df_pvalue14 = df_pvalue14.astype(float)
df_pvalue44 = df_pvalue44.astype(float)
#####

plt.figure(figsize=(30, 20))
plt.text(0.5, 0.5, 'Heatmap Relation between Dinucleotide and Position')

plt.subplot2grid((8, 21), (1, 0), rowspan=3, colspan=2)  # less cec1
plt.pcolor(df_pvalue12, cmap='jet', vmin=0, vmax=550)
plt.colorbar(ticks=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550])
plt.yticks(np.arange(0.5, len(df_pvalue12.index), 1), df_pvalue12.index)
plt.xticks(np.arange(0.5, len(df_pvalue12.columns), 1), df_pvalue12.columns, rotation='vertical')
plt.title('CEC1')
plt.ylabel('lesser test\ndinucleotide base', rotation='vertical')
plt.xlabel('position')

plt.subplot2grid((8, 21), (1, 12), rowspan=3, colspan=2)  # less cec2
plt.pcolor(df_pvalue42, cmap='jet', vmin=0, vmax=550)
plt.colorbar(ticks=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550])
plt.yticks(np.arange(0.5, len(df_pvalue42.index), 1), df_pvalue42.index)
plt.xticks(np.arange(0.5, len(df_pvalue42.columns), 1), df_pvalue42.columns, rotation='vertical')
plt.title('CEC2')
plt.ylabel('lesser test\ndinucleotide base', rotation='vertical')
plt.xlabel('position')

plt.subplot2grid((8, 21), (5, 0), rowspan=3, colspan=2)  # greater cec1
plt.pcolor(df_pvalue14, cmap='jet', vmin=0, vmax=550)
plt.colorbar(ticks=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550])
plt.yticks(np.arange(0.5, len(df_pvalue14.index), 1), df_pvalue14.index)
plt.xticks(np.arange(0.5, len(df_pvalue14.columns), 1), df_pvalue14.columns, rotation='vertical')
plt.title('CEC1')
plt.ylabel('greater test\ndinucleotide base', rotation='vertical')
plt.xlabel('position')

plt.subplot2grid((8, 21), (5, 12), rowspan=3, colspan=2)  # greater cec2
plt.pcolor(df_pvalue44, cmap='jet', vmin=0, vmax=550)
plt.colorbar(ticks=[0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550])
plt.yticks(np.arange(0.5, len(df_pvalue44.index), 1), df_pvalue44.index)
plt.xticks(np.arange(0.5, len(df_pvalue44.columns), 1), df_pvalue44.columns, rotation='vertical')
plt.title('CEC2')
plt.ylabel('greater test\ndinucleotide base', rotation='vertical')
plt.xlabel('position')

plt.show()
plt.savefig('heatmap_new_clustering_results_editminus.png')
plt.savefig('heatmap_new_clustering_results_editminus.svg')
