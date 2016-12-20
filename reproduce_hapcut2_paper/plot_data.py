# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 15:46:12 2015

@author: Peter Edge
"""

# based on example from:
# http://www.walkingrandomly.com/?p=5215
import matplotlib as mpl
mpl.use('Agg',force=True)
from matplotlib import pyplot as plt
import numpy as np

# some minor change

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


def plot_simulation(data,labels,x_list,xlabel,skip_list,outname,runtime_as_log=False,doplots=[1,1,1,1,1,1]):

    colors = ['k','b','y','r','m','deeppink']
    #colors = ['k','y','r']
    linewidths = [2,1.5,1.5,1.5,1.5,1.5]
    alphas = [1,0.5,0.5,0.5,0.5,0.5]

    xticklabels=None
    x_lo = min(x_list)
    x_hi = max(x_list)

    plot_generic(data,labels,skip_list,outname,runtime_as_log,colors,linewidths,alphas,x_list,xticklabels,xlabel,x_lo,x_hi,doplots=doplots)

def plot_experiment_sissor(data,labels,skip_list,outname,runtime_as_log=False,doplots=[1,1,0,0,1,1]):

    colors = ['b','y','r','m','deeppink','k']
    linewidths = [1.5,1.5,1.5,1.5,1.5,1.5]
    alphas = [1,1,1,1,1,1]

    x_list = list(range(1,23))
    xticklabels = [str(i) for i in range(1,23)]
    xlabel = 'Chromosome'
    x_lo = 0.9
    x_hi = 22.1

    plot_generic(data,labels,skip_list,outname,runtime_as_log,colors,linewidths,alphas,x_list,xticklabels,xlabel,x_lo,x_hi,plot_N50=True,doplots=doplots,MB=True)

def plot_experiment_hic(data,labels,skip_list,outname,runtime_as_log=False):

    colors = ['b','y','r','m','deeppink','k']
    linewidths = [3,3,3,3,3,3]
    alphas = [1,1,1,1,1,1]

    x_list = list(range(1,23))
    xticklabels = [str(i) for i in range(1,23)]
    #xticklabels.append('X')
    xlabel = 'Chromosome'
    x_lo = 0.9
    #x_hi = 23.1
    x_hi = 22.1

    plot_generic(data,labels,skip_list,outname,runtime_as_log,colors,linewidths,alphas,x_list,xticklabels,xlabel,x_lo,x_hi)


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

    mpl.rcParams['xtick.labelsize'] = 11
    mpl.rcParams['ytick.labelsize'] = 11

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
        plt.xlabel(xlabel, fontsize=9)
        plt.ylabel('Switch Error Rate', fontsize=9)
        if xticklabels:
            plt.xticks(x_list, xticklabels, rotation='vertical')
        plt.grid(True,color='silver')
        plt.xlim(x_lo, x_hi)

        #plt.legend(fancybox=True,loc='upper left',ncol=2,fontsize=10)
        plt.legend(bbox_to_anchor=(0., 1.05, 1., .105), loc=3,ncol=6,columnspacing=0.6,handletextpad=0.1,fontsize=13,borderaxespad=0.)

        # remove axis spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        # hiding axis ticks
        plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
        # add letter
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=15, fontweight='bold', va='top')

    # Plot Mismatch Rate
    if doplots[1]:
        pltnum += 1
        ax = plt.subplot(X1,X2, pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            else:
                plt.plot(x_data, mismatch_data[i], linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
        plt.xlabel(xlabel, fontsize=13)
        plt.ylabel('Mismatch Rate', fontsize=13)
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
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=15, fontweight='bold', va='top')

    # Plot Missing Rate
    if doplots[2]:
        pltnum += 1

        ax = plt.subplot(X1,X2, pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            else:
                plt.plot(x_data, missing_data[i], linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
        plt.xlabel(xlabel, fontsize=13)
        plt.ylabel('Fraction of Covered SNVs Pruned', fontsize=13)
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
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=15, fontweight='bold', va='top')

    # Plot max blk snp percent
    if doplots[3]:
        pltnum += 1

        ax = plt.subplot(X1,X2,pltnum)
        for i in range(0,l):
            if i in skip_list:
                plt.plot([],[],marker='.',color='w',label=labels[i]+' (N/A)')
            else:
                plt.plot(x_data, maxblk_data[i], linewidth=linewidths[i], color=colors[i],alpha=alphas[i],label=labels[i])
        plt.xlabel(xlabel, fontsize=13)
        plt.ylabel('Fraction of SNVs phased in largest block', fontsize=13)
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
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=15, fontweight='bold', va='top')

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
        plt.xlabel(xlabel, fontsize=13)
        if MB:
            plt.ylabel('AN50 (Mb)', fontsize=13)
        else:
            plt.ylabel('AN50', fontsize=13)
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
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=15, fontweight='bold', va='top')

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
        plt.xlabel(xlabel, fontsize=13)
        if plot_N50:
            if MB:
                plt.ylabel('N50 (Mb)', fontsize=13)
            else:
                plt.ylabel('N50', fontsize=13)
        else:
            plt.ylabel('Runtime (minutes)', fontsize=13)

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
        ax.text(0.05, 0.975, letters[pltnum-1], transform=ax.transAxes,fontsize=15, fontweight='bold', va='top')

    plt.tight_layout()
    plt.subplots_adjust(top=0.96)
    plt.savefig(outname,bbox_inches='tight')
