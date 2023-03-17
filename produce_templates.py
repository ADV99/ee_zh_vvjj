import ROOT
import os

# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()

# ___________________________________________________________________________________________________


def selection_dfs(df, proc, sel, h1s, h2s):

    print("{}: {}".format(sel["name"], sel["formula"]))
    df_sel = df.Filter(sel["formula"], sel["name"])
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

    return dfs, df_sel


def produce_graphs(proc, sels, h1s, h2s):

    df = proc.rdf()
    df = ROOT.RDataFrame("events", proc.files)
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

from config import processes, selection_tree, h1s, h2s

## run all selections with RDF
df_list = []
df_dict = dict()
for proc in processes:
    df_dict[proc.name] = produce_graphs(proc, selection_tree, h1s, h2s)
    for sel in selection_tree.get_all_nodes():
        df_list += df_dict[proc.name][sel.value["name"]]


ROOT.RDF.RunGraphs(df_list)

## now store histograms in ROOT files
#rdir = "/eos/user/s/selvaggi/analysis/zh_vvjj/templates2D/validation/"
rdir = "/eos/home-a/adelvecc/trying/"
ldir = "output"
os.system("mkdir -p {}".format(ldir))

for sel in selection_tree.children:

    fname = "{}/{}_hist2D.root".format(ldir, sel.value["name"])
    tf = ROOT.TFile.Open(
        fname,
        "RECREATE",
    )
    for proc in processes:
        for df in df_dict[proc.name][sel.value["name"]]:
            df.Write()

    os.system("cp {} {}".format(fname, rdir))
