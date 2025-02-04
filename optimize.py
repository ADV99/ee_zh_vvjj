import ROOT
import os
import numpy as np
#import matplotlib.pyplot as plt

# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()

# ___________________________________________________________________________________________________


def selection_dfs(df, proc, sel, h1s, h2s):

    print("{}: {}".format(sel["name"], sel["formula"]))
    df_sel = df.Filter(sel["formula"], sel["name"])
    counts = df_sel.Count()
    dfs = []
    for h1 in h1s:
        df_sel_h1 = df_sel.Histo1D(
            (
                # "{}_{}_{}".format(h1["name"], sel["name"], proc.name),
                "{}_{}".format(h1["name"], proc.name),
                "{};{};{}".format(h1["name"], h1["xtitle"], h1["ytitle"]),
                h1["nbins"],
                h1["xmin"],
                h1["xmax"],
            ),
            h1["var"],
            "weight_{}".format(proc.name),
        )
        dfs.append(df_sel_h1)

    for h2 in h2s:
        df_sel_h2 = df_sel.Histo2D(
            (
                # "{}_{}".format(h2["name"], sel["name"]),
                "{}_{}".format(h2["name"], proc.name),
                "{};{};{};{}".format(h2["name"], h2["xtitle"], h2["ytitle"], ""),
                h2["nbins_x"],
                h2["xmin"],
                h2["xmax"],
                h2["nbins_y"],
                h2["ymin"],
                h2["ymax"],
            ),
            h2["var_x"],
            h2["var_y"],
            "weight_{}".format(proc.name),
        )
        dfs.append(df_sel_h2)
    dfs.append(counts)

    return dfs, df_sel


def produce_graphs(proc, sels, h1s, h2s):

    df = proc.rdf()
    #df = ROOT.RDataFrame("events", proc.files)
    df = df.Define("weight_{}".format(proc.name), "{}".format(proc.weight))

    df_proc_dict = dict()

    def traverse_tree(df, df_dict, sel):
        df_dict[sel.value["name"]], df_sel = selection_dfs(
            df, proc, sel.value, h1s, h2s
        )
        for child in sel.children:
            traverse_tree(df_sel, df_dict, child)

    traverse_tree(df, df_proc_dict, sels)
    return df_proc_dict

# _____________________________________________________________________________________________
def compute_statistics(df_dict, processes):
    initial_counts = dict()
    weights = {}
    for pr in processes:
        initial_counts[pr.name] = pr.nevents
        weights[pr.name] = pr.weight
        #print(" > " + pr.name + " --> nevents: ", pr.nevents, " --> weight: ", pr.weight)

    stats = dict()
    #print("> Stats: ")
    for sel in list(df_dict.values())[0].keys():
        #print("-> " + sel + ": ")
        stats[sel] = dict()
        counting = 0

        for pr in df_dict.keys():
            stats[sel][pr] = dict()
            #print(df_dict[pr][sel][-1])
            stats[sel][pr]['Counts'] = df_dict[pr][sel][-1].GetValue()
            stats[sel][pr]["Efficiency"] = 1. * stats[sel][pr]['Counts'] / initial_counts[pr]
            stats[sel][pr]["Yields"] = stats[sel][pr]['Counts'] * weights[pr]

        for pr in df_dict.keys():

            stats[sel][pr]["Significance"] = stats[sel][pr]["Yields"] / np.sqrt( np.sum([stats[sel][p]["Yields"] for p in df_dict.keys()] ))

            #print("--> Process: ", pr)
            #print("Counts: ", stats[sel][pr]['Counts'])
            #print("Yields: ", stats[sel][pr]['Yields'])
            #print("Efficiency: ", stats[sel][pr]["Efficiency"])
            #print("Significance: ", stats[sel][pr]["Significance"])

            #print(stats)
            #print(stats.values())

    return stats

def combine_significances(stats, selecs, processes):
    result = []
    #_ = [print(stats[sel][proc]["Significance"]) for sel in selections]
    for proc in processes:
        s = [stats[sel][proc]["Significance"]**2 for sel in selecs]
        result.append(np.sqrt(np.sum(s)))
    return result


def print_table(outputDir, stats, signal = ["Hss"]):

    f = open(outputDir+"/outputTabular.txt","w")

    print('\\begin{table}[H] \n    \\centering \n    \\resizebox{\\textwidth}{!}{ \n    \\begin{tabular}{l',end='',file=f)
    print('c' * ( len(stats.values()) + 2 * len(signal) ) , end='' , file=f)
    print('} \n \\toprule',file=f)
    
    # Print headers
    row = [" "]
    for proc in list(stats.values())[0].keys():
        row.append(proc)
        if proc in signal:
            row.append(" ")
            row.append(" ")
    print(*row, sep = ' & ', end='', file=f)
    
    # Print content
    print('        \\\\ \n \midrule',file=f)
    for sel in stats.keys():
        row = [sel]
        for proc in list(stats.values())[0].keys():
            row.append("{:.2e}".format(stats[sel][proc]["Yields"]))
            if proc in signal:
                row.append("{:.3f}".format(stats[sel][proc]["Efficiency"]))
                row.append("{:.3g}".format(stats[sel][proc]["Significance"]))

        print(*row, sep = ' & ', end='', file=f)
        print(' \\\\ ', file=f)
    print('        \\bottomrule  \n    \\end{tabular}} \n    \\caption{Yields} \n    \\label{tab:my_label} \n\\end{table}', file=f)

    f.close()

# _____________________________________________________________________________________________
from config_optimize import processes, selection_tree, h1s, h2s

## run all selections with RDF producing hists and counts
df_list = []
df_dict = dict()
print(" --> Booking cuts: ")
for proc in processes:
    df_dict[proc.name] = produce_graphs(proc, selection_tree, h1s, h2s)
    for sel in selection_tree.get_all_nodes():
        df_list += df_dict[proc.name][sel.value["name"]]

print(" --> Doing cuts: ")
ROOT.RDF.RunGraphs(df_list)

# _____________________________________________________________________________________________

## Set output directory
#outdir = "/eos/home-a/adelvecc/winter2023/trying_BDT/output/"

## Compute Statistics and print tables
print(" --> Copmute statistics: ")
stats = compute_statistics(df_dict, processes)

interest = dict()
interest['Blike'] = 'Hbb'
interest['Clike'] = 'Hcc'
interest['Slike'] = 'Hss'
interest['Glike'] = 'Hgg'

finals_dict = dict()

#f_opt = open(outdir+"/output_opt.txt","w")
f_opt = open("/output_opt.txt","w")

print("\n --> Computing significances : ")
for c1 in selection_tree.children:
    print(" -> ", c1.value["name"])
    finals_dict[c1.value["name"]] = dict()
    # Compute significances for all limits
    for c2 in c1.children:
        print(" > ", c2.value["name"])
        final_selections = [s.value["name"] for s in c2.children]
        finals_dict[c1.value["name"]][c2.value["name"]] = combine_significances(stats, final_selections,[interest[c1.value['name']]])

    # look for the best significance
    sort_ed = sorted(list(finals_dict[c1.value["name"]].keys()), key = lambda x : finals_dict[c1.value["name"]][x][0])
    #print(sort_ed)
    #print(finals_dict)
    print("Best cuts for " + c1.value["name"] + " : " , sort_ed[-1])
    print("Best significance for " + c1.value["name"] + " : ", finals_dict[c1.value["name"]][sort_ed[-1]])
    print("Best cuts for " + c1.value["name"] + " : " , sort_ed[-1], file = f_opt)
    print("Best significance for " + c1.value["name"] + " : ", finals_dict[c1.value["name"]][sort_ed[-1]], file = f_opt)

f_opt.close()
