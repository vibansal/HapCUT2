# based on example from:
# http://www.walkingrandomly.com/?p=5215
import matplotlib as mpl
mpl.use('Agg',force=True)
import error_rates
import fileIO
from matplotlib import pyplot as plt
import numpy as np
chroms = list(range(1,23))+['X']
from collections import defaultdict
from estimate_htrans_probs_known_phase import estimate_htrans_probs
import os
import sys
import pickle

bad_paths = ['/opt/scipy/2.7/lib/python2.7/site-packages/d2to1-0.2.11-py2.7.egg', '/opt/scipy/2.7/lib/python2.7/site-packages/stsci.distutils-0.3.7-py2.7.egg', '/opt/scipy/2.7/lib/python2.7/site-packages/pyfits-3.3-py2.7-linux-x86_64.egg','/opt/scipy/2.7/lib/python2.7/site-packages']

for p in bad_paths:
    if p in sys.path:
        sys.path.remove(p)

from scipy.signal import savgol_filter

mpl.rc('legend', fontsize=9)
mpl.rc('xtick', labelsize=9)
mpl.rc('ytick', labelsize=9)
mpl.rc('axes', labelsize=9)
mpl.rc('axes', labelsize=9)
mpl.rcParams.update({'font.size': 9})
mpl.rc('lines', linewidth=1.5)
mpl.rc('mathtext',default='regular')

z1 = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
z2 = [0.1195, 0.5, 1, 2, 3, 4, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160]
z3 = ['sim{}'.format(x) for x in range(0,20)]
z4 = list(range(0,10))
z5 = 'hapcut2_paper/plots/simulation_runtime_comparison.png'

def plot_simulation_runtime_comparison(C_covs=z1,L_lens=z2,sims=z3,reps = z4,runtime_fig=z5):


    # measure the mean reads crossing a SNP from the simulated Hi-C reads
    # we have to do this because the simulation actually parametrized the maximum span
    # but on the final plot we're showing mean reads crossing a SNP
    S_spans_file = 'hapcut2_paper/data/S_spans.p'
    if not os.path.exists(S_spans_file):
        S_spans = []
        for sim in sims:
            cov_dict = defaultdict(int)
            frag_matrix = 'hapcut2_paper/data/sim_vary_span/{}.0'.format(sim)

            with open(frag_matrix,"r") as fm:
                for line in fm:
                    if len(line) < 2:
                        continue

                    el = line.strip().split()

                    num_blks      = int(el[0])

                    call_list  = el[2:(2+2*num_blks)]              # extract base call part of line
                    call_list  = zip(*[iter(call_list)]*2)             # list -> tuple list conversion: credit to http://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
                    call_list  = [(int(a)-1, b) for a,b in call_list]  # convert index to 0-based integer
                    call_list2 = []

                    for ix, blk in call_list:
                        curr_ix = ix
                        for a in blk:
                            call_list2.append((curr_ix, a))
                            curr_ix += 1

                    firstpos = call_list2[0][0]
                    lastpos  = call_list2[-1][0]

                    for i in range(firstpos, lastpos+1):
                        cov_dict[i] += 1

            span_counts = list(cov_dict.values())
            mean_num = sum(span_counts)/len(span_counts)
            S_spans.append(mean_num)

        pickle.dump(S_spans,open(S_spans_file,'wb'))
    else:
        S_spans = pickle.load(open(S_spans_file,'rb'))


    C_hapcut2  = []
    C_probhap  = []
    C_refhap   = []
    C_hapcut1  = []
    #C_dgs      = []
    C_fasthare = []

    L_hapcut2 = []
    L_probhap = []
    L_refhap  = []
    L_hapcut1 = []
    #L_dgs      = []
    L_fasthare = []

    S_hapcut2 = []
    S_probhap = []
    S_refhap  = []
    S_hapcut1 = []
    #S_dgs  = []
    S_fasthare = []

    n_reps = len(reps)

    for lab in sims:

        C_hapcut2.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_coverage/{}.{}/hapcut2.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        C_refhap.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_coverage/{}.{}/refhap.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        C_probhap.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_coverage/{}.{}/probhap.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        C_hapcut1.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_coverage/{}.{}/hapcut1.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        #C_dgs.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_coverage/{}.{}/dgs.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        C_fasthare.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_coverage/{}.{}/fasthare.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])

        L_hapcut2.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_length/{}.{}/hapcut2.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        L_refhap.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_length/{}.{}/refhap.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        L_probhap.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_length/{}.{}/probhap.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        L_hapcut1.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_length/{}.{}/hapcut1.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        #L_dgs.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_length/{}.{}/dgs.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)]) # temp fix, one rep crashed...
        L_fasthare.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_length/{}.{}/fasthare.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])

        S_hapcut2.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_span/{}.{}/hapcut2.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        S_refhap.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_span/{}.{}/refhap.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        S_probhap.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_span/{}.{}/probhap.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        S_hapcut1.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_span/{}.{}/hapcut1.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        #S_dgs.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_span/{}.{}/dgs.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])
        S_fasthare.append([fileIO.parse_runtime_file('hapcut2_paper/experiments/sim_vary_span/{}.{}/fasthare.runtime'.format(lab,rep))/3600 for rep in range(0,n_reps)])

    C_hapcut2  = np.array(C_hapcut2)
    C_probhap  = np.array(C_probhap)
    C_refhap   = np.array(C_refhap)
    C_hapcut1  = np.array(C_hapcut1)
    #C_dgs      = np.array(C_dgs)
    C_fasthare = np.array(C_fasthare)

    L_hapcut2 = np.array(L_hapcut2)
    L_probhap = np.array(L_probhap)
    L_refhap  = np.array(L_refhap)
    L_hapcut1 = np.array(L_hapcut1)
    #L_dgs  = np.array(L_dgs)
    L_fasthare = np.array(L_fasthare)

    S_hapcut2  = np.array(S_hapcut2)
    S_probhap  = np.array(S_probhap)
    S_refhap   = np.array(S_refhap)
    S_hapcut1  = np.array(S_hapcut1)
    #S_dgs      = np.array(S_dgs)
    S_fasthare = np.array(S_fasthare)

    def mn(arr):
        if len(np.shape(arr)) == 2:
            return np.mean(arr,1)
        else:
            return np.mean(arr)

    def up(arr):
        return np.mean(arr,1) + np.std(arr,1)

    def dn(arr):
        return np.mean(arr,1) - np.std(arr,1)

    def myplt(x,y,c,l,xstop=None,mem_err=False):
        if xstop == None:
            xstop = len(x)
        plt.plot(x[:xstop],mn(y[:xstop]), color=c,alpha=1,linestyle='-',label=l)
        plt.fill_between(x[:xstop], up(y[:xstop]), dn(y[:xstop]), alpha=0.25, edgecolor=c, facecolor=c)
        if xstop != len(x):
            if not mem_err:
                plt.scatter([x[xstop-1]],[mn(y[xstop-1])],color=c,marker='x',s=100,lw=3)
            else:
                plt.scatter([x[xstop-1]],[mn(y[xstop-1])],color=c,marker='o',s=100,lw=3)

    def find_runtime_limit(data):
        stop = 0
        for i in data:
            stop += 1
            if np.any(i >= 10.0):
                break

        if np.sum(i >= 10.0) > 8:
            stop -= 1 # the error rates are probably wonky in this case
        return stop

    # doesn't really find memory limit
    # but this is only really used for HapCUT1
    # so we just look for when the runtime plummets
    # because when it memorys out, it just returns immediately
    def find_memory_limit(data):
        stop = 0
        last_i = 0
        for i in data:
            if np.any(i < (last_i / 100)):
                break
            stop += 1
            last_i = i
        return stop

    # find where refhap and probhap hit the runtime limit
    C_refhap_stop = find_runtime_limit(C_refhap)
    C_probhap_stop = find_runtime_limit(C_probhap)
    C_hapcut1_stop = find_memory_limit(C_hapcut1)

    # credit to this post for showing how to use fill-between: http://stackoverflow.com/questions/12957582/matplotlib-plot-yerr-xerr-as-shaded-region-rather-than-error-bars

    fig = plt.figure(figsize=(7.5,4.5))

    ax1=fig.add_subplot(231)

    myplt(C_covs,C_hapcut2,'k','HapCUT2')
    myplt(C_covs,C_hapcut1,'b','HapCUT',C_hapcut1_stop,mem_err=True)
    #myplt(C_covs,C_dgs,'g','DGS')
    myplt(C_covs,C_fasthare,'m','FastHare')
    myplt(C_covs,C_refhap,'y','RefHap',C_refhap_stop)
    myplt(C_covs,C_probhap,'r','ProbHap',C_probhap_stop)

    #plt.title("Completeness")
    plt.xlabel("Coverage Per SNV (d)")
    plt.ylabel("Runtime (hours)")
    plt.yscale("log")
    plt.xlim(5,101)
    plt.ylim(0.001,13)

    plt.grid(True,color='grey')
    ax1.set_yticks([0.01,0.1,1,10])
    ax1.set_yticklabels(['0.01','0.1','1','10'])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax1.text(-0.07, -0.07, 'A', transform=ax1.transAxes,fontsize=10, fontweight='bold', va='top')


    ax2 = fig.add_subplot(232)
    plt.grid(True,color='grey')
    ax2.tick_params(axis=u'both', which=u'both',length=0)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    ax2.text(-0.07, -0.07, 'B', transform=ax2.transAxes,fontsize=10, fontweight='bold', va='top')


    # HapCUT1 experiences memory error rather than runtime errors, find where to mark
    L_hapcut1_stop = find_memory_limit(L_hapcut1)

    #myplt(L_lens,L_dgs,'g','DGS')
    myplt(L_lens,L_fasthare,'m','FastHare')
    myplt(L_lens,L_hapcut1,'b','HapCUT',L_hapcut1_stop,mem_err=True)
    myplt(L_lens,L_refhap,'y','RefHap')
    myplt(L_lens,L_probhap,'r','ProbHap')
    myplt(L_lens,L_hapcut2,'k','HapCUT2')

    #plt.title("Accuracy")
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.set_xlabel("SNVs Per Read (V)")
    plt.ylabel("Runtime (hours)")


    plt.xlim(0.001,200000)
    plt.ylim(0.001,13)
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    xtick1 = [0.5,1,2,4,10,20,40,80,160]
    xticklab1 = [0.5,1,2,4,10,20,40,80,160]
    ax2.set_xlim(0,160)
    ax2.set_xticks(xtick1)
    ax2.set_xticklabels(xticklab1)
    ax2.set_yticks([0.01,0.1,1,10])
    ax2.set_yticklabels(['0.01','0.1','1','10'])

    S_refhap_stop = find_runtime_limit(S_refhap)
    S_probhap_stop = find_runtime_limit(S_probhap)
    #S_dgs_stop = find_runtime_limit(S_dgs)
    S_fasthare_stop = find_runtime_limit(S_fasthare)

    ax4=fig.add_subplot(233)


    #myplt(S_spans,S_dgs,'g','DGS',S_dgs_stop)
    myplt(S_spans,S_fasthare,'m','FastHare',S_fasthare_stop)
    myplt(S_spans,S_hapcut1,'b','HapCUT')
    myplt(S_spans,S_refhap,'y','RefHap',S_refhap_stop)
    myplt(S_spans,S_probhap,'r','ProbHap',S_probhap_stop)
    myplt(S_spans,S_hapcut2,'k','HapCUT2')


    plt.xlabel("Mean Mate-Pairs Crossing an SNV (d')")
    plt.ylabel("Runtime (hours)")
    plt.yscale("log")
    plt.xlim(min(S_spans),max(S_spans))
    plt.ylim(0.001,13)
    plt.grid(True,color='grey')

    ax4.set_yticks([0.01,0.1,1,10])
    ax4.set_yticklabels(['0.01','0.1','1','10'])
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax4.text(-0.07, -0.07, 'C', transform=ax4.transAxes,fontsize=10, fontweight='bold', va='top')

    #plt.subplots_adjust(top=0.9)
    #plt.tight_layout()
    #plt.savefig(runtime_fig,bbox_inches='tight')

    #############################################################################################
    # now make error plot
    ############################################################################################

    def get_switch_mismatch_err(lst_of_lsts):
        new_lst = []
        for lst in lst_of_lsts:
            new_err = []
            for err in lst:
                new_err.append(err.get_switch_mismatch_rate())
            new_lst.append(new_err)
        return new_lst

    def prep(lst_of_lsts):
        return np.array(get_switch_mismatch_err(lst_of_lsts))

    [C_hapcut2, C_hapcut1, C_refhap, C_probhap, C_dgs, C_fasthare] = [prep(z) for z in  pickle.load(open('hapcut2_paper/error_rates/tool_comparisons/sim_vary_coverage.stats.p','rb'))]
    [L_hapcut2, L_hapcut1, L_refhap, L_probhap, L_dgs, L_fasthare] = [prep(z) for z in pickle.load(open('hapcut2_paper/error_rates/tool_comparisons/sim_vary_length.stats.p','rb'))]
    [S_hapcut2, S_hapcut1, S_refhap, S_probhap, S_dgs, S_fasthare] = [prep(z) for z in pickle.load(open('hapcut2_paper/error_rates/tool_comparisons/sim_vary_span.stats.p','rb'))]
    # the dgs variables in the above list should be deleted later

    ax1=fig.add_subplot(234)
    #myplt(C_covs,C_dgs,'g','DGS')
    myplt(C_covs,C_fasthare,'m','FastHare')
    myplt(C_covs,C_hapcut1,'b','HapCUT',C_hapcut1_stop,mem_err=True)
    myplt(C_covs,C_refhap,'y','RefHap',C_refhap_stop)
    myplt(C_covs,C_probhap,'r','ProbHap',C_probhap_stop)
    myplt(C_covs,C_hapcut2,'k','HapCUT2')


    plt.scatter([],[],color='grey',marker='o',s=100,lw=3,label='Memory Limit')
    plt.scatter([],[],color='grey',marker='x',s=100,lw=3,label='Runtime Limit')

    #plt.title("Completeness")
    plt.xlabel("Coverage Per SNV (d)")
    plt.ylabel("Switch+Mismatch Error Rate")
    #plt.yscale("log")
    plt.xlim(5,101)
    plt.ylim(0,0.002)

    plt.grid(True,color='grey')
    #ax1.set_yticks([0.01,0.1,1,10])
    #ax1.set_yticklabels(['0.01','0.1','1','10'])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax1.text(-0.07,-0.07, 'D', transform=ax1.transAxes,fontsize=10, fontweight='bold', va='top')
    #plt.legend(loc='upper right',scatterpoints = 1,ncol=2,columnspacing=0.25)


    ax2 = fig.add_subplot(235)
    plt.grid(True,color='grey')
    ax2.tick_params(axis=u'both', which=u'both',length=0)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    ax2.text(-0.07, -0.07, 'E', transform=ax2.transAxes,fontsize=10, fontweight='bold', va='top')


    #myplt(L_lens,L_dgs,'g','DGS')
    myplt(L_lens,L_fasthare,'m','FastHare')
    myplt(L_lens,L_hapcut1,'b','HapCUT',L_hapcut1_stop,mem_err=True)
    myplt(L_lens,L_refhap,'y','RefHap')
    myplt(L_lens,L_probhap,'r','ProbHap')
    myplt(L_lens,L_hapcut2,'k','HapCUT2')

    #plt.title("Accuracy")
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.set_xlabel("SNVs Per Read (V)")
    plt.ylabel("Switch+Mismatch Error Rate")


    plt.xlim(0.001,200000)
    plt.ylim(0,0.00275)
    ax2.set_xscale("log")
    #ax2.set_yscale("log")
    #ax3.set_xscale("log")
    #ax3.set_yscale("log")

    xtick1 = [0.5,1,2,4,10,20,40,80,160]
    xticklab1 = [0.5,1,2,4,10,20,40,80,160]
    ax2.set_xlim(0,160)
    ax2.set_xticks(xtick1)
    ax2.set_xticklabels(xticklab1)
    #ax2.set_yticks([0.01,0.1,1,10])
    #ax2.set_yticklabels(['0.01','0.1','1','10'])


    ax4=fig.add_subplot(236)

    #myplt(S_spans,S_dgs,'g','DGS',S_dgs_stop)
    myplt(S_spans,S_hapcut2,'k','HapCUT2')
    myplt(S_spans,S_hapcut1,'b','HapCUT')
    #plt.plot([],[],color='g',label='DGS')
    #plt.plot([],[],color='m',label='FastHare')
    myplt(S_spans,S_refhap,'y','RefHap',S_refhap_stop)
    myplt(S_spans,S_probhap,'r','ProbHap',S_probhap_stop)
    myplt(S_spans,S_fasthare,'m','FastHare',S_fasthare_stop)

    plt.xlabel("Mean Mate-Pairs Crossing an SNV (d')")
    plt.ylabel("Switch+Mismatch Error Rate")
    #plt.yscale("log")
    plt.xlim(min(S_spans),max(S_spans))
    plt.ylim(0.0,0.02)
    plt.grid(True,color='grey')

    #ax4.set_yticks([0.01,0.1,1,10])
    #ax4.set_yticklabels(['0.01','0.1','1','10'])
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax4.text(-0.07, -0.07, 'F', transform=ax4.transAxes,fontsize=10, fontweight='bold', va='top')

    plt.scatter([],[],color='grey',marker='o',s=100,lw=3,label='Mem. Limit')
    plt.scatter([],[],color='grey',marker='x',s=100,lw=3,label='Time Limit')
    plt.legend(loc='lower right',scatterpoints = 1,columnspacing=0.15)

    #plt.subplots_adjust(top=0.9)
    plt.tight_layout()
    plt.savefig(runtime_fig,bbox_inches='tight')


if __name__ == '__main__':
    plot_simulation_runtime_comparison()

def plot_hic_completeness(mboI_covs,hindIII_covs,hapblocks_mboI,hapblocks_hindIII,vcf_file,output_fig):

    #mboI_covs = [i*11 for i in range(1,20)]
    #hindIII_covs = [i*11 for i in range(1,20)]
    cov_cut = 0
    for i,x in enumerate(hindIII_covs):
        cov_cut = i+1
        if x >= 200:
            break

    mboI_covs = mboI_covs[:cov_cut]
    hindIII_covs = hindIII_covs[:cov_cut]
    hapblocks_mboI = hapblocks_mboI[:cov_cut]
    hapblocks_hindIII = hapblocks_hindIII[:cov_cut]

    # this isn't the correct frag file but it's only used for getting missing rate and we aren't looking at that
    frag_file = None
    runtime_file = None

    errs_mboI = []
    errs_hindIII = []

    for assembly_file in hapblocks_mboI:
        err = error_rates.hapblock_vcf_error_rate(assembly_file, frag_file, vcf_file, runtime_file,largest_blk_only=True)
        errs_mboI.append(err)

    for assembly_file in hapblocks_hindIII:
        err = error_rates.hapblock_vcf_error_rate(assembly_file, frag_file, vcf_file, runtime_file,largest_blk_only=True)
        errs_hindIII.append(err)


    switch_rates_mboI = []
    mismatch_rates_mboI = []
    switch_mismatch_rates_mboI = []
    max_blks_mboI = []
    AN50_mboI = []
    for err in errs_mboI:
        switch_rates_mboI.append(err.get_switch_rate())
        mismatch_rates_mboI.append(err.get_mismatch_rate())
        switch_mismatch_rates_mboI.append(err.get_switch_mismatch_rate())
        #switch_mismatch_rates_mboI.append(err.get_flat_error_rate())
        max_blks_mboI.append(err.get_max_blk_snp_percent())
        AN50_mboI.append(err.get_AN50())


    switch_rates_hindIII = []
    mismatch_rates_hindIII = []
    max_blks_hindIII = []
    switch_mismatch_rates_hindIII = []

    AN50_hindIII = []
    for err in errs_hindIII:
        switch_rates_hindIII.append(err.get_switch_rate())
        mismatch_rates_hindIII.append(err.get_mismatch_rate())
        switch_mismatch_rates_hindIII.append(err.get_switch_mismatch_rate())
        #switch_mismatch_rates_hindIII.append(err.get_flat_error_rate())
        max_blks_hindIII.append(err.get_max_blk_snp_percent())
        AN50_hindIII.append(err.get_AN50())

    plt.figure(figsize=(5,2.5))

    ax1=plt.subplot(1,2,1)
    plt.plot(mboI_covs,max_blks_mboI, color='k',linestyle='-',label='MboI')
    plt.plot(hindIII_covs,max_blks_hindIII, color='r',linestyle='-',label='HindIII')
    #plt.scatter(cov_selvaraj,[err_selvaraj.max_blk['chr1']],marker='o',facecolors='b', edgecolors='b',s=100,label='HindIII (Selvaraj et al)')
    plt.xlim(0,hindIII_covs[cov_cut-1])
    plt.ylim(0,1)
    plt.xlabel("Read Coverage")
    plt.ylabel("Fraction of SNVs in largest block")
    #plt.legend(loc='lower right',scatterpoints = 1)
    plt.grid(True,color='grey')
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    #ax1.spines["bottom"].set_visible(False)
    #ax1.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax1.text(0, -0.09, 'A', transform=ax1.transAxes,fontsize=10, fontweight='bold', va='top')

    ax2=plt.subplot(1,2,2)
    plt.plot(mboI_covs,switch_mismatch_rates_mboI, color='k',linestyle='-',label='MboI')
    plt.plot(hindIII_covs,switch_mismatch_rates_hindIII, color='r',linestyle='-',label='HindIII')
    #plt.scatter(cov_selvaraj,[1-err_selvaraj.get_switch_rate()],marker='o',facecolors='b', edgecolors='b',s=100,label='HindIII (Selvaraj et al)')
    plt.xlim(0,hindIII_covs[cov_cut-1])
    #plt.ylim(0.975,1)
    plt.xlabel("Read coverage")
    plt.ylabel("Switch + Mismatch Error Rate for largest block")
    plt.tight_layout()
    #plt.legend(loc='lower right',scatterpoints = 1)
    plt.grid(True,color='grey')
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    #ax2.spines["bottom"].set_visible(False)
    #ax2.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax2.text(0, -0.09, 'B', transform=ax2.transAxes,fontsize=10, fontweight='bold', va='top')
    plt.legend(loc='upper right',scatterpoints = 1)

    plt.tight_layout()
    plt.savefig(output_fig,bbox_inches='tight')


# paths to datafiles are hardcoded in and then repeated in the Snakefile
# too many for it to be reasonable pass them in as parameters
def plot_common_snps_error(output_fig):

    from error_rates import hapblock_vcf_error_rate_COMMON_ALL_CHROM as common_error
    chroms = ['chr{}'.format(x) for x in range(1,23)]+['chrX']

    # FOSMID COMMON SNVS ERROR
    tool_list1 = ['hapcut','refhap','hapcut2', 'probhap', 'fasthare']
    A_list1 = [['hapcut2_paper/experiments/fosmid/{}/hapcut1.output'.format(x),'hapcut2_paper/experiments/fosmid/{}/refhap.output'.format(x),'hapcut2_paper/experiments/fosmid/{}/hapcut2.output'.format(x),'hapcut2_paper/experiments/fosmid/{}/probhap.output.uncorrected'.format(x),'hapcut2_paper/experiments/fosmid/{}/fasthare.output'.format(x)] for x in chroms]
    R_list1 = [['hapcut2_paper/experiments/fosmid/{}/hapcut1.runtime'.format(x),'hapcut2_paper/experiments/fosmid/{}/refhap.runtime'.format(x),'hapcut2_paper/experiments/fosmid/{}/hapcut2.runtime'.format(x),'hapcut2_paper/experiments/fosmid/{}/probhap.runtime'.format(x),'hapcut2_paper/experiments/fosmid/{}/fasthare.runtime'.format(x)] for x in chroms]
    frag_files1 = ['hapcut2_paper/data/fosmid/{}'.format(x) for x in chroms]
    vcf_files1 = ['hapcut2_paper/data/NA12878_VCFs_hg18/{}.vcf'.format(x) for x in chroms]

    fosmid_common = common_error(tool_list1,A_list1,R_list1,frag_files1,vcf_files1,dataset_name='fosmid')

    # 10X COMMON SNVS ERROR
    tool_list2 = ['hapcut','hapcut2','fasthare']
    A_list2 = [['hapcut2_paper/experiments/10X/{}/hapcut1.output'.format(x),'hapcut2_paper/experiments/10X/{}/hapcut2.output'.format(x),'hapcut2_paper/experiments/10X/{}/fasthare.output'.format(x)] for x in chroms]
    R_list2 = [['hapcut2_paper/experiments/10X/{}/hapcut1.runtime'.format(x),'hapcut2_paper/experiments/10X/{}/hapcut2.runtime'.format(x),'hapcut2_paper/experiments/10X/{}/fasthare.runtime'.format(x)] for x in chroms]
    frag_files2 = ['hapcut2_paper/data/10X/{}'.format(x) for x in chroms]
    vcf_files2 = ['hapcut2_paper/data/NA12878_VCFs_hg19/{}.vcf'.format(x) for x in chroms]

    tenX_common = common_error(tool_list2,A_list2,R_list2,frag_files2,vcf_files2,dataset_name='10X')

    # PACBIO 44x COMMON SNVS ERROR
    tool_list3 = ['hapcut','refhap','hapcut2','fasthare']
    A_list3 = [['hapcut2_paper/experiments/pacbio44/{}/hapcut1.output'.format(x),'hapcut2_paper/experiments/pacbio44/{}/refhap.output'.format(x),'hapcut2_paper/experiments/pacbio44/{}/hapcut2.output'.format(x),'hapcut2_paper/experiments/pacbio44/{}/fasthare.output'.format(x)] for x in chroms]
    R_list3 = [['hapcut2_paper/experiments/pacbio44/{}/hapcut1.runtime'.format(x),'hapcut2_paper/experiments/pacbio44/{}/refhap.runtime'.format(x),'hapcut2_paper/experiments/pacbio44/{}/hapcut2.runtime'.format(x),'hapcut2_paper/experiments/pacbio44/{}/fasthare.runtime'.format(x)] for x in chroms]
    frag_files3 = ['hapcut2_paper/data/pacbio44/{}'.format(x) for x in chroms]
    vcf_files3 = ['hapcut2_paper/data/NA12878_VCFs_hg19/{}.vcf'.format(x) for x in chroms]

    pacbio44_common = common_error(tool_list3,A_list3,R_list3,frag_files3,vcf_files3,dataset_name='pacbio44')

    # PACBIO 11x COMMON SNVS ERROR
    tool_list4 = ['hapcut','refhap','hapcut2','probhap','fasthare']
    A_list4 = [['hapcut2_paper/experiments/pacbio11/{}/hapcut1.output'.format(x),'hapcut2_paper/experiments/pacbio11/{}/refhap.output'.format(x),'hapcut2_paper/experiments/pacbio11/{}/hapcut2.output'.format(x),'hapcut2_paper/experiments/pacbio11/{}/probhap.output.uncorrected'.format(x),'hapcut2_paper/experiments/pacbio11/{}/fasthare.output'.format(x)] for x in chroms]
    R_list4 = [['hapcut2_paper/experiments/pacbio11/{}/hapcut1.runtime'.format(x),'hapcut2_paper/experiments/pacbio11/{}/refhap.runtime'.format(x),'hapcut2_paper/experiments/pacbio11/{}/hapcut2.runtime'.format(x),'hapcut2_paper/experiments/pacbio11/{}/probhap.runtime'.format(x),'hapcut2_paper/experiments/pacbio11/{}/fasthare.runtime'.format(x)] for x in chroms]
    frag_files4 = ['hapcut2_paper/data/pacbio11/{}'.format(x) for x in chroms]
    vcf_files4 = ['hapcut2_paper/data/NA12878_VCFs_hg19/{}.vcf'.format(x) for x in chroms]

    pacbio11_common = common_error(tool_list4,A_list4,R_list4,frag_files4,vcf_files4,dataset_name='pacbio11')

    # HIC 30x COMMON SNVS ERROR
    tool_list5 = ['hapcut','hapcut2']
    cov = 30
    print("COVERAGE 30")
    A_list5 = [['hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut/{}.output'.format(cov,x),'hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut2/{}.output'.format(cov,x)] for x in chroms]
    R_list5 = [['hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut/{}.runtime'.format(cov,x),'hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut2/{}.runtime'.format(cov,x)] for x in chroms]
    frag_files5 = ['hapcut2_paper/data/hic_mboI_subsamples/cov{}/{}'.format(cov,x) for x in chroms]
    vcf_files5 = ['hapcut2_paper/data/NA12878_VCFs_hg19/{}.vcf'.format(x) for x in chroms]
    hic30_common = common_error(tool_list5,A_list5,R_list5,frag_files5,vcf_files5,dataset_name='hic30')

    # HIC 40x COMMON SNVS ERROR
    tool_list6 = ['hapcut','hapcut2']
    cov= 40
    print("COVERAGE 40")
    A_list6 = [['hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut/{}.output'.format(cov,x),'hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut2/{}.output'.format(cov,x)] for x in chroms]
    R_list6 = [['hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut/{}.runtime'.format(cov,x),'hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut2/{}.runtime'.format(cov,x)] for x in chroms]
    frag_files6 = ['hapcut2_paper/data/hic_mboI_subsamples/cov{}/{}'.format(cov,x) for x in chroms]
    vcf_files6 = ['hapcut2_paper/data/NA12878_VCFs_hg19/{}.vcf'.format(x) for x in chroms]
    hic40_common = common_error(tool_list6,A_list6,R_list6,frag_files6,vcf_files6,dataset_name='hic40')

    # HIC 90x COMMON SNVS ERROR
    tool_list7 = ['hapcut','hapcut2']
    cov= 90
    print("COVERAGE 90")
    A_list7 = [['hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut/{}.output'.format(cov,x),'hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut2/{}.output'.format(cov,x)] for x in chroms]
    R_list7 = [['hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut/{}.runtime'.format(cov,x),'hapcut2_paper/experiments/hic/hapcut_htrans/mboI/cov{}/hapcut2/{}.runtime'.format(cov,x)] for x in chroms]
    frag_files7 = ['hapcut2_paper/data/hic_mboI_subsamples/cov{}/{}'.format(cov,x) for x in chroms]
    vcf_files7 = ['hapcut2_paper/data/NA12878_VCFs_hg19/{}.vcf'.format(x) for x in chroms]
    hic90_common = common_error(tool_list7,A_list7,R_list7,frag_files7,vcf_files7,dataset_name='hic90')

    #[fosmid_common,tenX_common,pacbio44_common,pacbio11_common,hic30_common,hic40_common,hic90_common] = pickle.load(open('temp.p','rb'))
    #pickle.dump(lst,open('temp.p','wb'))

    # adapted from http://matplotlib.org/examples/api/barchart_demo.html
    #!/usr/bin/env python
    # a bar plot with errorbars

    # plot it

    alpha1 = 0.8
    alpha2 = 0.5

    plt.figure(figsize=(7.5,5))

    #####################################################
    # PLOT FOSMID
    #####################################################

    ax = plt.subplot(211)
    IND = np.array([0,1,2,3,4,5])  # the x locations for the groups


    # corrected for same snps
    refhap_errs = (fosmid_common['refhap'].get_switch_rate(), fosmid_common['refhap'].get_mismatch_rate())
    probhap_errs = (fosmid_common['probhap'].get_switch_rate(), fosmid_common['probhap'].get_mismatch_rate())
    hapcut_errs = (fosmid_common['hapcut'].get_switch_rate(), fosmid_common['hapcut'].get_mismatch_rate())
    hapcut2_errs = (fosmid_common['hapcut2'].get_switch_rate(), fosmid_common['hapcut2'].get_mismatch_rate())
    fasthare_errs = (fosmid_common['fasthare'].get_switch_rate(), fosmid_common['fasthare'].get_mismatch_rate())

    ind = IND[0:2]#np.arange(N)  # the x locations for the groups
    width = 0.11


    plt.bar(ind, hapcut2_errs, color='k',
            ecolor='black', # black error bar color
            alpha=alpha1,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT2')
    plt.bar(ind+width, probhap_errs, color='r',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='ProbHap')
    plt.bar(ind+2*width, refhap_errs, color='y',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='RefHap')
    plt.bar(ind+3*width, fasthare_errs, color='m',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='FastHare')
    plt.bar(ind+4*width, hapcut_errs, color='b',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT')



    # add some text for labels, title and axes ticks
    ax.set_ylabel('Error Rate')
    ax.set_xticks(IND+1.5*width)
    ax.set_xticklabels(('Switch\n                                          Fosmid','Mismatch','Switch\n                                     PacBio (11'r'$\times$'')','Mismatch','Switch\n                                      PacBio (44'r'$\times$'')','Mismatch'))

    ax.yaxis.grid(True,color='grey')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    plt.legend(loc='upper left')
    plt.xlim(-0.25,5.75)



    # corrected for same snps
    # groups in order are fosmid, pacbio11, pacbio44
    ax.text(0.015, -0.07, 'A', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')
    ax.text(0.34, -0.07, 'B', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')


    plt.ylim((0,0.011))

    #####################################################
    # PLOT PACBIO11
    #####################################################

    #IND = np.array([-2, 0.5, 1, 1.5, 2.0, 2.5, 3, 5.25]) #np.array([0,1,2,3,4,5])

    ind = IND[2:4] #np.arange(N)  # the x locations for the groups


    #width = 0.25       # the width of the bars
    width = 0.11
    # corrected for same snps
    refhap_errs = (pacbio11_common['refhap'].get_switch_rate(), pacbio11_common['refhap'].get_mismatch_rate())
    probhap_errs = (pacbio11_common['probhap'].get_switch_rate(), pacbio11_common['probhap'].get_mismatch_rate())
    hapcut_errs = (pacbio11_common['hapcut'].get_switch_rate(), pacbio11_common['hapcut'].get_mismatch_rate())
    hapcut2_errs = (pacbio11_common['hapcut2'].get_switch_rate(), pacbio11_common['hapcut2'].get_mismatch_rate())
    fasthare_errs = (pacbio11_common['fasthare'].get_switch_rate(), pacbio11_common['fasthare'].get_mismatch_rate())

    def label(rect, height, value):
        # attach some text labels
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%0.4f' % value,
                ha='center', va='bottom')

    plt.bar(ind, hapcut2_errs, color='k',
            ecolor='black', # black error bar color
            alpha=alpha1,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT2')
    plt.bar(ind+width, probhap_errs, color='r',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='ProbHap')
    plt.bar(ind+2*width, refhap_errs, color='y',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='RefHap')
    plt.bar(ind+3*width, fasthare_errs, color='m',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='FastHare')

    hapcut_errs = (0.0104,hapcut_errs[1]) # make room to label hapcut manually

    rects = plt.bar(ind+4*width, hapcut_errs, color='b',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT2')

    label(rects[0],0.01,hapcut_errs[0])

    #####################################################
    # PLOT PACBIO44
    #####################################################

    ind = IND[4:6]#np.arange(N)  # the x locations for the groups
    #width = 0.12


    refhap_errs = (pacbio44_common['refhap'].get_switch_rate(), pacbio44_common['refhap'].get_mismatch_rate())
    hapcut2_errs = (pacbio44_common['hapcut2'].get_switch_rate(), pacbio44_common['hapcut2'].get_mismatch_rate())
    hapcut_errs = (pacbio44_common['hapcut'].get_switch_rate(), pacbio44_common['hapcut'].get_mismatch_rate())
    fasthare_errs = (pacbio44_common['fasthare'].get_switch_rate(), pacbio44_common['fasthare'].get_mismatch_rate())

    plt.bar(ind, hapcut2_errs, color='k',
            ecolor='black', # black error bar color
            alpha=alpha1,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT2')
    plt.bar(ind+width, refhap_errs, color='y',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='Refhap')
    plt.bar(ind+2*width, fasthare_errs, color='m',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='FastHare')
    plt.bar(ind+3*width, hapcut_errs, color='b',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT')

    #plt.legend(loc='upper left')

    #####################################################
    # PLOT 10X Genomics
    #####################################################

    ax = plt.subplot(212)

    # corrected for same snps
    hapcut_errs = (tenX_common['hapcut'].get_switch_rate(), tenX_common['hapcut'].get_mismatch_rate())
    hapcut2_errs = (tenX_common['hapcut2'].get_switch_rate(), tenX_common['hapcut2'].get_mismatch_rate())
    fasthare_errs = (tenX_common['fasthare'].get_switch_rate(), tenX_common['fasthare'].get_mismatch_rate())

    ind = IND[0:2]#np.arange(N)  # the x locations for the groups
    width = 0.16


    plt.bar(ind+0.5*width, hapcut2_errs, color='k',
            ecolor='black', # black error bar color
            alpha=alpha1,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT2')
    plt.bar(ind+1.5*width, hapcut_errs, color='b',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT')
    plt.bar(ind+2.5*width, fasthare_errs, color='m',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='FastHare')

    # corrected for same snps
    # groups in order are fosmid, pacbio11, pacbio44
    #ax.text(-0.035, -0.035, 'A', transform=ax.transAxes,fontsize=15, fontweight='bold', va='top')

    #####################################################
    # PLOT ALL HIC
    #####################################################

    hapcut2_errs = (hic40_common['hapcut2'].get_switch_rate(), hic40_common['hapcut2'].get_mismatch_rate(),hic90_common['hapcut2'].get_switch_rate(), hic90_common['hapcut2'].get_mismatch_rate())
    hapcut1_errs = (hic40_common['hapcut'].get_switch_rate(), hic40_common['hapcut'].get_mismatch_rate(),hic90_common['hapcut'].get_switch_rate(), hic90_common['hapcut'].get_mismatch_rate())

    ind = IND[2:6] # the x locations for the groups
    width = 0.18
    #width = 0.12       # the width of the bars

    plt.bar(ind+width, hapcut2_errs, color='k',
            ecolor='black', # black error bar color
            alpha=alpha1,      # transparency
            width=width,      # smaller bar width
            align='center',
            label='HapCUT2')
    plt.bar(ind+2*width, hapcut1_errs, color='b',
            ecolor='black', # black error bar color
            alpha=alpha2,      # transparency
            width=width,      # smaller bar width
            align='center',
           label='HapCUT')

    # add some text for labels, title and axes ticks
    #plt.xlim(-0.5,6)

    ax.set_ylabel('Error Rate')
    ax.set_xticks(IND+1.5*width)
    ax.set_xticklabels(('Switch\n                                     10X Genomics','Mismatch','Switch\n                                      Hi-C (40'r'$\times$'')','Mismatch','Switch\n                                          Hi-C (90'r'$\times$'')','Mismatch'))

    ax.yaxis.grid(True,color='grey')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")

    ax.text(0.015, -0.07, 'C', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')
    ax.text(0.34, -0.07, 'D', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')

    plt.xlim(-0.25,5.75)

    plt.tight_layout()
    plt.savefig(output_fig,bbox_inches='tight')


        #plt.legend()
    def percent_improved(old, new):
        return str((old-new)/old*100)[:4]+'%'

    print("HapCUT2 switch improvement over ProbHap on Fosmid Data")
    print(percent_improved(fosmid_common['probhap'].get_switch_rate(),fosmid_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over ProbHap on Fosmid Data")
    print(percent_improved(fosmid_common['probhap'].get_mismatch_rate(),fosmid_common['hapcut2'].get_mismatch_rate()))
    print("")
    print("HapCUT2 switch improvement over ProbHap on PacBio11 Data")
    print(percent_improved(pacbio11_common['probhap'].get_switch_rate(),pacbio11_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over ProbHap on PacBio11 Data")
    print(percent_improved(pacbio11_common['probhap'].get_mismatch_rate(),pacbio11_common['hapcut2'].get_mismatch_rate()))
    print("")
    print("HapCUT2 switch improvement over RefHap on PacBio11 Data")
    print(percent_improved(pacbio11_common['refhap'].get_switch_rate(),pacbio11_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over RefHap on PacBio11 Data")
    print(percent_improved(pacbio11_common['refhap'].get_mismatch_rate(),pacbio11_common['hapcut2'].get_mismatch_rate()))
    print("")
    print("HapCUT2 switch improvement over RefHap on PacBio44 Data")
    print(percent_improved(pacbio44_common['refhap'].get_switch_rate(),pacbio44_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over RefHap on PacBio44 Data")
    print(percent_improved(pacbio44_common['refhap'].get_mismatch_rate(),pacbio44_common['hapcut2'].get_mismatch_rate()))
    print("")

    print("HapCUT2 switch improvement over FastHare on Fosmid Data")
    print(percent_improved(fosmid_common['fasthare'].get_switch_rate(),fosmid_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over FastHare on Fosmid Data")
    print(percent_improved(fosmid_common['fasthare'].get_mismatch_rate(),fosmid_common['hapcut2'].get_mismatch_rate()))
    print("")
    print("HapCUT2 switch improvement over FastHare on PacBio11 Data")
    print(percent_improved(pacbio11_common['fasthare'].get_switch_rate(),pacbio11_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over FastHare on PacBio11 Data")
    print(percent_improved(pacbio11_common['fasthare'].get_mismatch_rate(),pacbio11_common['hapcut2'].get_mismatch_rate()))
    print("")
    print("HapCUT2 switch improvement over FastHare on PacBio44 Data")
    print(percent_improved(pacbio44_common['fasthare'].get_switch_rate(),pacbio44_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over FastHare on PacBio44 Data")
    print(percent_improved(pacbio44_common['fasthare'].get_mismatch_rate(),pacbio44_common['hapcut2'].get_mismatch_rate()))
    print("")
    print("HapCUT2 switch improvement over FastHare on 10X Data")
    print(percent_improved(tenX_common['fasthare'].get_switch_rate(),tenX_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over FastHare on 10X Data")
    print(percent_improved(tenX_common['fasthare'].get_mismatch_rate(),tenX_common['hapcut2'].get_mismatch_rate()))
    print("")
    print("HapCUT2 switch improvement over HapCUT on 10X Data")
    print(percent_improved(tenX_common['hapcut'].get_switch_rate(),tenX_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over HapCUT on 10X Data")
    print(percent_improved(tenX_common['hapcut'].get_mismatch_rate(),tenX_common['hapcut2'].get_mismatch_rate()))
    print("")

    print("HapCUT2 switch improvement over HapCUT on hic30 Data")
    print(percent_improved(hic30_common['hapcut'].get_switch_rate(),hic30_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over HapCUT on hic30 Data")
    print(percent_improved(hic30_common['hapcut'].get_mismatch_rate(),hic30_common['hapcut2'].get_mismatch_rate()))
    print("")
    print("HapCUT2 switch improvement over HapCUT on hic40 Data")
    print(percent_improved(hic40_common['hapcut'].get_switch_rate(),hic40_common['hapcut2'].get_switch_rate()))
    print("HapCUT2 mismatch improvement over HapCUT on HiC40 Data")
    print(percent_improved(hic40_common['hapcut'].get_mismatch_rate(),hic40_common['hapcut2'].get_mismatch_rate()))
    print("")

    print("HapCUT2 minutes of runtime on fosmid data:")
    print(fosmid_common['hapcut2'].get_runtime()/60)
    print("HapCUT2 minutes of runtime on pacbio11 data:")
    print(pacbio11_common['hapcut2'].get_runtime()/60)
    print("HapCUT2 minutes of runtime on pacbio44 data:")
    print(pacbio44_common['hapcut2'].get_runtime()/60)
    print("HapCUT2 minutes of runtime on hic30 data:")
    print(hic30_common['hapcut2'].get_runtime()/60)
    print("HapCUT2 hours of runtime on hic90 data:")
    print(hic90_common['hapcut2'].get_runtime()/3600)
    print("HapCUT2 minutes of runtime on 10X data:")
    print(tenX_common['hapcut2'].get_runtime()/60)
    print("HapCUT hours of runtime on 10X data:")
    print(tenX_common['hapcut'].get_runtime()/3600)
    print("FastHare hours of runtime on 10X data:")
    print(tenX_common['fasthare'].get_runtime()/3600)

def plot_pacbio_vs_hic1(vcf_files19, h_p11_files, h_p44_files, hic40_files, hic90_files, output_fig):

    def per_distance_error_rate_genome(files, vcf_files, FRAC, BINSIZE):

        res = error_rates.PDER_result([],[], binsize=BINSIZE, sample_frac=FRAC,max_dist=0)
        for file,vcf in zip(files,vcf_files):
            res += error_rates.per_distance_error_rate(file, vcf,binsize=BINSIZE, sample_frac=FRAC)
        return res

    # pacbio 11
    h_p11 = per_distance_error_rate_genome(h_p11_files,vcf_files19,0.01,20000)
    ## pacbio 44
    h_p44 = per_distance_error_rate_genome(h_p44_files,vcf_files19,0.01,20000)
    # hic 40
    hic40 = per_distance_error_rate_genome(hic40_files,vcf_files19,0.0001,20000)
    # hic 90
    hic90 = per_distance_error_rate_genome(hic90_files,vcf_files19,0.0001,20000)


    vcf_file = 'hapcut2_paper/data/NA12878_VCFs_hg19/chr1.vcf'
    hapblock_list_pacbio11 = fileIO.parse_hapblock_file(h_p11_files[0])
    num_blks_pacbio11, frac_phased_pacbio11 = error_rates.frac_SNPs_per_num_blks(hapblock_list_pacbio11, vcf_file, use_SNP_index=True)
    hapblock_list_pacbio44= fileIO.parse_hapblock_file(h_p44_files[0])
    num_blks_pacbio44, frac_phased_pacbio44 = error_rates.frac_SNPs_per_num_blks(hapblock_list_pacbio44, vcf_file,use_SNP_index=True)

    hapblock_list_mboI40 = fileIO.parse_hapblock_file(hic40_files[0])
    num_blks_mboI40, frac_phased_mboI40 = error_rates.frac_SNPs_per_num_blks(hapblock_list_mboI40, vcf_file, use_SNP_index=True)
    hapblock_list_mboI90 = fileIO.parse_hapblock_file(hic90_files[0])
    num_blks_mboI90, frac_phased_mboI90 = error_rates.frac_SNPs_per_num_blks(hapblock_list_mboI90, vcf_file, use_SNP_index=True)

    plt.figure(figsize=(5,2.5))

    ############################################################################################

    ax3 = plt.subplot(1,2,1)
    plt.plot(num_blks_mboI90,frac_phased_mboI90, color='k',linestyle='-',label=r'Hi-C, 90$\times$')
    plt.plot(num_blks_mboI40,frac_phased_mboI40, color='k',linestyle='--',label=r'Hi-C, 40$\times$')
    plt.plot(num_blks_pacbio44,frac_phased_pacbio44, color='b',linestyle='-',label=r'PacBio, 44$\times$')
    plt.plot(num_blks_pacbio11,frac_phased_pacbio11, color='b',linestyle='--',label=r'PacBio, 11$\times$')


    plt.ylim(0,1)
    plt.xscale('log')
    plt.xlabel("Number of Largest Blocks")
    plt.ylabel("Fraction of SNVs Phased")
    #plt.legend(loc='lower right')
    plt.grid(True,color='grey')
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    #ax.spines["bottom"].set_visible(False)
    #ax.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax3.text(-0.08, -0.05, 'A', transform=ax3.transAxes,fontsize=10, fontweight='bold', va='top')

    ############################################################################################


    ax4 = plt.subplot(1,2,2)

    # fosmid needs to be cut off at the first point it hits 0.5
    p11_cut = 1
    for stat in h_p11.get_error_stats():
        if stat <= 0.5:
            break
        p11_cut += 1

    plt.plot(hic90.get_bin_starts(),hic90.get_error_stats(), color='k',linestyle='-',label=r'Hi-C, 90$\times$')
    plt.plot(hic40.get_bin_starts(),hic40.get_error_stats(), color='k',linestyle='--',label=r'Hi-C, 40$\times$')
    plt.plot(h_p44.get_bin_starts(),h_p44.get_error_stats(), color='b',linestyle='-',label=r'PacBio, 44$\times$')
    plt.plot(h_p11.get_bin_starts()[0:p11_cut],h_p11.get_error_stats()[0:p11_cut], color='b',linestyle='--',label=r'PacBio, 11$\times$')

    plt.xlim(0,int(1e6))
    plt.ylim(0.5,1)
    plt.xticks([200000,400000,600000,800000,1000000], [200,400,600,800,1000])
    plt.xlabel("Distance Between SNVs (Kb)")
    plt.ylabel("Fraction of SNV Pairs Phased Correctly")
    plt.grid(True,color='grey')
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    #box = ax4.get_position()
    #ax4.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    #ax4.legend(loc='center left', ncol=2,bbox_to_anchor=(0.05, 1.13))
    #ax4.legend
    ax4.text(-0.05, -0.05, 'B', transform=ax4.transAxes,fontsize=10, fontweight='bold', va='top')
    #plt.xlim((0,100000))
    plt.legend(bbox_to_anchor=(0.6, 0.8), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(output_fig,bbox_inches='tight')



def plot_pacbio_vs_hic2(pacbio_hapblocks, hapblocks_mboI, hapblocks_hindIII, vcf_file, pacbio_covs, mboI_covs, hindIII_covs, output_fig):

    ############################################################################################
    plt.figure(figsize=(5,2.5))
    # pacbio error_rates

    # this isn't the correct frag file but it's only used for getting missing rate and we aren't looking at that
    pacbio_frag_file = None # we dont need missing rate
    pacbio_runtime_file = None # we don't need runtime

    pacbio_errs = []

    for assembly_file in pacbio_hapblocks:
        err = error_rates.hapblock_vcf_error_rate(assembly_file, pacbio_frag_file, vcf_file, pacbio_runtime_file)
        pacbio_errs.append(err)

    # this isn't the correct frag file but it's only used for getting missing rate and we aren't looking at that
    frag_file = None
    runtime_file = None

    errs_mboI = []
    errs_hindIII = []

    for assembly_file in hapblocks_mboI:
        err = error_rates.hapblock_vcf_error_rate(assembly_file, frag_file, vcf_file, runtime_file)
        errs_mboI.append(err)

    for assembly_file in hapblocks_hindIII:
        err = error_rates.hapblock_vcf_error_rate(assembly_file, frag_file, vcf_file, runtime_file)
        errs_hindIII.append(err)

    pacbio_switch_rates = []
    pacbio_mismatch_rates = []

    for err in pacbio_errs:
        pacbio_switch_rates.append(err.get_switch_rate())
        pacbio_mismatch_rates.append(err.get_mismatch_rate())

    ax1=plt.subplot(1,2,1)
    plt.plot(pacbio_covs,pacbio_switch_rates,  color='b',linestyle='-',label='PacBio, Switch')
    plt.plot(pacbio_covs,pacbio_mismatch_rates,  color='b',linestyle='--',label='PacBio, Mismatch')

    plt.ylim(0,0.045)
    plt.xlabel("Read coverage")
    plt.ylabel("Error Rate")
    plt.tight_layout()
    #plt.legend(loc='lower right',scatterpoints = 1)
    plt.grid(True,color='grey')
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    #ax2.spines["bottom"].set_visible(False)
    #ax2.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax1.text(-0.05, -0.05, 'A', transform=ax1.transAxes,fontsize=10, fontweight='bold', va='top')
    plt.legend(loc='upper right',scatterpoints = 1)


def plot_pruning_pacbio11(hapcut2_blocks,vcf_file,frag_file,output_fig):

    # import fileIO

    h_err1 = []
    runtime_file = None # we don't care about runtime either
    output_file = None
    expdir = 'hapcut2_paper/experiments/pacbio_error_analysis'
    if not os.path.exists(expdir):
        os.makedirs(expdir)

    #vary SNV confidence for hapcut2
    err = error_rates.hapblock_vcf_error_rate(hapcut2_blocks, frag_file, vcf_file, runtime_file,use_SNP_index=True)
    h_err1.append(err)
    snp_cutoffs = [0.75,0.9,0.95,0.99,0.999,0.9999,0.99999]
    for i,snp_cutoff in enumerate(snp_cutoffs):
        output_file = '{}/hapcut2.prune{}'.format(expdir,i)
        fileIO.prune_hapblock_file(hapcut2_blocks, output_file, snp_cutoff, -1, False)
        err = error_rates.hapblock_vcf_error_rate(output_file, frag_file, vcf_file, runtime_file,use_SNP_index=True)
        h_err1.append(err)

    #vary SNP confidence for probhap
    p_err1 = []
    probhap_blocks = 'hapcut2_paper/experiments/pacbio11/chr1/probhap.output'
    emission_cutoffs = [0,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01]
    for i,emission_cutoff in enumerate(emission_cutoffs):
        output_file = 'hapcut2_paper/experiments/pacbio_error_analysis/probhap.prune{}'.format(i)
        fileIO.prune_probhap_file(probhap_blocks, output_file, emission_cutoff, -1)
        err = error_rates.hapblock_vcf_error_rate(output_file, frag_file, vcf_file, runtime_file,use_SNP_index=True)
        p_err1.append(err)

    h_err2 = []
    split_cutoffs = [-1,0.5,0.75,0.9,0.99,0.999,0.9999,0.99999,0.999999,0.9999999,0.99999999999]
    for i,split_cutoff in enumerate(split_cutoffs):
        output_file = '{}/hapcut2.split{}'.format(expdir,i)
        fileIO.prune_hapblock_file(hapcut2_blocks, output_file, -1, split_cutoff, False)
        err = error_rates.hapblock_vcf_error_rate(output_file, frag_file, vcf_file, runtime_file,use_SNP_index=True)
        h_err2.append(err)

    p_err2 = []
    probhap_blocks = 'hapcut2_paper/experiments/pacbio11/chr1/probhap.output.uncorrected'
    split_cutoffs2 = [-1,0.5,0.6,0.7,0.8,0.9,0.99,0.999,0.9999,0.99999,0.999999]
    for i,split_cutoff in enumerate(split_cutoffs2):
        output_file = 'hapcut2_paper/experiments/pacbio_error_analysis/probhap.split{}'.format(i)
        fileIO.prune_probhap_file(probhap_blocks, output_file, -1, split_cutoff)
        err = error_rates.hapblock_vcf_error_rate(output_file, frag_file, vcf_file, runtime_file,use_SNP_index=True)
        p_err2.append(err)

    refhap_blocks = 'hapcut2_paper/experiments/pacbio11/chr1/refhap.output'
    r_err = error_rates.hapblock_vcf_error_rate(refhap_blocks, frag_file, vcf_file, runtime_file,use_SNP_index=True)

    plt.figure(figsize=(5,2.5))
    ax1 = plt.subplot(121)
    plt.plot([err.get_missing_rate() for err in h_err1],[err.get_mismatch_rate() for err in h_err1],color='k',label='HapCUT2')
    plt.plot([err.get_missing_rate() for err in p_err1],[err.get_mismatch_rate() for err in p_err1],color='r',alpha=0.5,label='ProbHap')
    plt.scatter([r_err.get_missing_rate()],[r_err.get_mismatch_rate()],s=100,color='y',alpha=0.5,label='RefHap')
    #plt.ylim(0.5,1)
    plt.xlim(0,0.175)
    plt.ylim(0,0.0105)

    plt.xlabel("Fraction of Covered SNVs Pruned")
    plt.ylabel("Mismatch Error Rate")
    plt.grid(True,color='grey')
    ax1.set_xticks([0,0.04,0.08,0.12,0.16])
    #ax1.set_xticklabels(['0.01','0.1','1','10'])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    ax1.text(-0.075, -0.05, 'A', transform=ax1.transAxes,fontsize=10, fontweight='bold', va='top')

    plt.legend(scatterpoints = 1)

    ax2 = plt.subplot(122)
    plt.plot([err.get_AN50()/1000 for err in h_err2],[err.get_switch_rate() for err in h_err2],color='k',label='HapCUT2')
    plt.plot([err.get_AN50()/1000 for err in p_err2],[err.get_switch_rate() for err in p_err2],color='r',alpha=0.5,label='ProbHap')
    plt.scatter([r_err.get_AN50()/1000],[r_err.get_switch_rate()],s=100,color='y',alpha=0.5,label='RefHap')

    plt.xlim(0,int(1e2))
    plt.ylim(0,0.008)
    plt.xlabel("AN50 (Kb)")
    plt.ylabel("Switch Error Rate")
    plt.grid(True,color='grey')
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.text(-0.05, -0.05, 'B', transform=ax2.transAxes,fontsize=10, fontweight='bold', va='top')

    plt.tight_layout()
    plt.savefig(output_fig,bbox_inches='tight')

    #plt.tight_layout()
    #plt.legend(loc='upper left',scatterpoints = 1)


def plot_hic_htrans_estimates(chr1_em_file, chr19_em_file, chr1_gt_file, chr19_gt_file, chr1_vcf, chr19_vcf, chr1_frags, chr19_frags, output_fig):

    def read_distfile(f):
        result = []
        with open(f) as infile:
            for line in infile:
                el = line.strip().split()
                result.append(float(el[1]))

        return result

    estimate_htrans_probs(chr1_frags,chr1_vcf,chr1_gt_file,bin_size=1000000)
    estimate_htrans_probs(chr19_frags,chr19_vcf,chr19_gt_file,bin_size=1000000)

    def smooth(l):
        return savgol_filter(l,501,1)

    chr1_gt = read_distfile(chr1_gt_file)
    chr1_em = read_distfile(chr1_em_file)
    chr19_gt = read_distfile(chr19_gt_file)
    chr19_em = read_distfile(chr19_em_file)

    xdata1 = list(range(0,42000000,1000000))
    xdata2 = list(range(0,42000000,5000))

    plt.figure(figsize=(2.5,2.5))
    ax = plt.subplot(111)
    plt.plot(xdata1[:len(chr1_gt)],chr1_gt, color='b',linestyle='-',label='Chr1, Ground truth')
    plt.plot(xdata2[:len(chr1_em)],smooth(chr1_em),  color='b',linestyle='dotted',label='Chr1, HapCUT2')
    plt.plot(xdata1[:len(chr19_gt)],chr19_gt,  color='r',linestyle='-',label='Chr19, Ground truth')
    plt.plot(xdata2[:len(chr19_em)],smooth(chr19_em), color='r',linestyle='dotted',label='Chr19, HapCUT2')
    #plt.title("90x Hi-C H-trans probabilities: Ground truth estimate vs HapCUT2 estimate")
    plt.grid(True,color='grey')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    plt.xlabel("Insert Size, I (Mb)")
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ylabel(r'Estimate of $\tau$(I)')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.ylim(0,0.25)
    plt.xticks(list(range(0,44000000,5000000)),list(range(0,44,5)))
    plt.tight_layout()
    plt.savefig(output_fig,bbox_inches='tight')

def plot_generic(data,labels,skip_list,outname,runtime_as_log,colors,linewidths,alphas,x_list,xticklabels,xlabel,x_lo,x_hi,plot_N50=False,doplots=[1,1,1,1,1,1],cutoffdict=dict(),MB=False,hic=False):

    ############################################################################################
    # READ DATA FROM EXPERIMENT BENCHMARK FILE INTO LISTS
    ############################################################################################
    l = len(labels)
    assert (l == len(data))

    switches   = []
    mismatches = []
    missings   = []
    AN50s      = []
    N50s       = []
    maxblks    = []
    runtimes   = []

    x_len = len(x_list)
    for i in range(0,l):
        switches.append([None]*x_len)
        mismatches.append([None]*x_len)
        missings.append([None]*x_len)
        AN50s.append([None]*x_len)
        N50s.append([None]*x_len)
        maxblks.append([None]*x_len)
        runtimes.append([None]*x_len)

    #masks = [np.zeros(len(x_list)) for i in range(0,l)]

    for i, result_list in enumerate(data):

        if result_list == None:
            continue

        for j, res in enumerate(result_list):

            missing = res.get_missing_rate()
            #if missing == 1.0:
            #    masks[i][j] = 1 # we'll leave this one off the plot
            if j >= x_len:
                continue
            switches[i][j] = res.get_switch_rate()
            mismatches[i][j] = res.get_mismatch_rate()
            missings[i][j] = missing
            AN50s[i][j] = res.get_AN50()
            N50s[i][j] = res.get_N50()
            maxblks[i][j] = res.get_max_blk_snp_percent()
            runtimes[i][j] = res.get_runtime()/60.0

    ############################################################################################
    # TURN LISTS INTO NUMPY ARRAYS
    ############################################################################################

    x_data = np.array(x_list)
    switch_data   = []
    mismatch_data = []
    missing_data  = []
    AN50_data     = []
    N50_data     = []
    maxblk_data   = []
    runtime_data  = []

    for i in range(0,l):
        switch_data.append(np.array(switches[i]))
        mismatch_data.append(np.array(mismatches[i]))
        missing_data.append(np.array(missings[i]))
        AN50_data.append(np.array(AN50s[i]))
        N50_data.append(np.array(N50s[i]))
        maxblk_data.append(np.array(maxblks[i]))
        runtime_data.append(np.array(runtimes[i]))

    ############################################################################################
    # PLOT THE DATA
    ############################################################################################

    def get_axis_limits(ax, xscale=0.05, yscale=.95): # credit to http://stackoverflow.com/questions/24125058/add-label-to-subplot-in-matplotlib
        return ax.get_xlim()[1]*xscale, ax.get_ylim()[1]*yscale


    if sum(doplots) == 3:
        X1 = 1
        X2 = 3
        fs = (5,2.5)
    elif sum(doplots)>4:
        X1 = 3
        X2 = 2
        fs = (5,7.5)
    else:
        X1 = 2
        X2 = 2
        fs = (5,5)

    plt.figure(num=None, figsize=fs)

    letters = 'ABCDEF'

    pltnum = 0
    # Plot Switch Errors
    if doplots[0]:
        pltnum += 1
        ax = plt.subplot(X1,X2,pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            else:
                plt.plot(x_data, switch_data[i], linewidth=linewidths[i], color=colors[i],alpha=alphas[i], label=labels[i])
        plt.xlabel(xlabel)
        plt.ylabel('Switch Error Rate')
        if xticklabels:
            plt.xticks(x_list, xticklabels, rotation='vertical')
        plt.grid(True,color='silver')
        plt.xlim(x_lo, x_hi)

        plt.legend(bbox_to_anchor=(0., 1.05, 1., .105), loc=3,ncol=6,columnspacing=0.6,handletextpad=0.1,borderaxespad=0.)

        # remove axis spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        # hiding axis ticks
        plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
        # add letter
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')

    # Plot Mismatch Rate
    if doplots[1]:
        pltnum += 1
        ax = plt.subplot(X1,X2, pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            else:
                plt.plot(x_data, mismatch_data[i], linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
        plt.xlabel(xlabel)
        plt.ylabel('Mismatch Rate')
        if xticklabels:
            plt.xticks(x_list, xticklabels, rotation='vertical')
        plt.grid(True,color='silver')
        plt.xlim(x_lo, x_hi)
        # remove axis spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        # hiding axis ticks
        plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
        # add letter
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')

    # Plot Missing Rate
    if doplots[2]:
        pltnum += 1

        ax = plt.subplot(X1,X2, pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            else:
                plt.plot(x_data, missing_data[i], linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
        plt.xlabel(xlabel)
        plt.ylabel('Fraction of Covered SNVs Pruned')
        if xticklabels:
            plt.xticks(x_list, xticklabels, rotation='vertical')
        plt.grid(True,color='silver')
        plt.xlim(x_lo, x_hi)
        # remove axis spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        # hiding axis ticks
        plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
        # add letter
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')

    # Plot max blk snp percent
    if doplots[3]:
        pltnum += 1

        ax = plt.subplot(X1,X2,pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            else:
                plt.plot(x_data, maxblk_data[i], linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
        plt.xlabel(xlabel)
        plt.ylabel('Fraction of SNVs phased in largest block')
        if xticklabels:
            plt.xticks(x_list, xticklabels, rotation='vertical')
        plt.grid(True,color='silver')
        plt.xlim(x_lo, x_hi)
        if hic:
            plt.ylim(0,1.1)
        # remove axis spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        # hiding axis ticks
        plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
        # add letter
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')

    # Plot AN50
    if doplots[4]:
        pltnum += 1

        ax = plt.subplot(X1,X2,pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            else:
                if MB:
                    plt.plot(x_data, [j/1000000 for j in AN50_data[i]], linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
                else:
                    plt.plot(x_data, AN50_data[i], linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
        plt.xlabel(xlabel)
        if MB:
            plt.ylabel('AN50 (Mb)')
        else:
            plt.ylabel('AN50')
        if xticklabels:
            plt.xticks(x_list, xticklabels, rotation='vertical')
        plt.grid(True,color='silver')
        plt.xlim(x_lo, x_hi)
        # remove axis spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        # hiding axis ticks
        plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
        # add letter
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')

    # Plot Runtime in minutes
    if doplots[5]:
        pltnum += 1

        ax = plt.subplot(X1,X2, pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            elif plot_N50:
                if MB:
                    plt.plot(x_data, [j/1000000 for j in N50_data[i]],linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
                else:
                    plt.plot(x_data, N50_data[i],linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
            else:
                plt.plot(x_data, runtime_data[i],linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
        plt.xlabel(xlabel)
        if plot_N50:
            if MB:
                plt.ylabel('N50 (Mb)')
            else:
                plt.ylabel('N50')
        else:
            plt.ylabel('Runtime (minutes)')

        if runtime_as_log and not plot_N50:
            plt.yscale('log')
        if xticklabels:
            plt.xticks(x_list, xticklabels, rotation='vertical')
        plt.grid(True,color='silver')
        plt.xlim(x_lo, x_hi)
        # remove axis spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        # hiding axis ticks
        plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
        # add letter
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=10, fontweight='bold', va='top')

    plt.tight_layout()
    plt.subplots_adjust(top=0.96)
    plt.savefig(outname,bbox_inches='tight')

# borrowed publication plot code from
# http://nbviewer.jupyter.org/github/rasbt/matplotlib-gallery/blob/master/ipynb/publication.ipynb

def plot_experiment(data,labels,skip_list,outname,runtime_as_log=False,doplots=[1,1,1,1,1,1],hic=False):

    colors = ['k','b','y','r','m','g']
    linewidths = [3,3,3,3,3,3]
    alphas = [1,0.5,0.5,0.5,0.5,0.5]

    x_list = list(range(1,24))
    xticklabels = [str(i) for i in range(1,23)]
    xticklabels.append('X')
    xlabel = 'Chromosome'
    x_lo = 0.9
    x_hi = 23.1

    plot_generic(data,labels,skip_list,outname,runtime_as_log,colors,linewidths,alphas,x_list,xticklabels,xlabel,x_lo,x_hi,doplots=doplots,hic=hic)
