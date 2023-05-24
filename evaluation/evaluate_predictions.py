import numpy as np
import pandas as pd
import math
from ast import literal_eval
from argparse import ArgumentParser 

#python3 evaluate_predictions.py  -path_results_to_write evaluation.csv -path_results_to_read res4.csv



# compute overlaps between units
def overlap_units(l1, l2) :

    overlap_list = []
    dic_units = {}
    for index, units1 in enumerate(l1) :
        dic_units.setdefault(("unit" + str(index), len(units1)), [])
        list_overlap = []
        for units2 in list(l2) :
            overlap = ([], [], [], False, list(units1), list(units2))
            intersection = units1.intersection(units2)
            if (len(intersection) > 0) :
                intersect = units1.intersection(units2)
                overlap = (list(intersect), list(units1.difference(units2)), list(units2.difference(units1)), True, list(units1), list(units2))
                list_overlap.append(overlap)
            else :
                list_overlap.append(overlap)
        dic_units[("unit" + str(index), len(units1))] = list_overlap
    return dic_units


def check_overlap(dic_units) :
    maximum_overlap_list = []
    matched = []
    for key, values in sorted(dic_units.items(), key=lambda tup : tup[0][1], reverse = True) :
        flag = False
        index = 0
        while flag == False and index < len(values) :
            sorted_tup = sorted(values, key=lambda tup: len(tup[0]), reverse = True)
            if sorted_tup[index][-1] not in matched and len(sorted_tup[index][0]) > 0 :
                flag = True
                matched.append(list(sorted_tup[index][-1]))
                maximum_overlap_list.append(sorted_tup[index])
            else :
                flag = False 
                maximum_overlap_list.append(sorted_tup[index - 1])
            index = index + 1
    if flag == False :
        maximum_overlap_list.append(sorted_tup[index - 1])

    print(maximum_overlap_list)
    TP = []
    FN = []
    FP = []
    FPFN = []
    flag = False
    found = []
    for l in maximum_overlap_list :
        if (l[3] == True) :
            found.extend(l[0])
            TP.extend(l[0])
            for el in l[1] :
                if el not in found :
                    flag = False
                    for units in l1 :
                        if(el in units) :
                            FN.append(el)
                            flag = True
                            break
            for el in l[2] :
                if el not in found :
                    flag = False
                    for units in l2 :
                        if(el in units) :
                            print(el)
                            print("----------")
                            FP.append(el)
                            flag = True
                            break

    for l in maximum_overlap_list :
        if (l[3] == False) :
            for el in l[-2] :
                if el in FP:
                    FPFN.append(el)
                elif el not in FP and el not in FN :
                    FN.append(el)
    for el in list(l2) :
        el = list(el)
        for x in el :
            if x not in found :
                FP.append(x)

    intersection = set(FP).intersection(set(FN))
    FPFN = list(set(FPFN).union(intersection))
    FP = list(set(FP).difference(set(FPFN)))
    FN = list(set(FN).difference(set(FPFN)))
    TN = len(uniprot) - (len(TP) + len(FP) + len(FN) + len(FPFN))
    accuracy = np.round((len(TP) + TN) / (len(TP) + TN + len(FP) + len(FN) + len(FPFN)), 2)
    precision = np.round((len(TP)) / (len(TP) + len(FP) + len(FPFN)), 2)
    recall = np.round((len(TP)) / (len(TP) + len(FN)), 2)
    print("TP : {}".format(TP))
    print("FP : {}".format(FP))
    print("TN : {}".format(TN))
    print("FN : {}".format(FN))
    print("FPFN : {}".format(FPFN))
    #print("TP : {}".format(len(TP)))
    #print("FP : {}".format(len(FP)))
    #print("TN : {}".format(TN))
    #print("FN : {}".format(len(FN)))
    #print("FPFN : {}".format(len(FPFN)))
    #print("accuracy : {}".format(accuracy))
    #print("precision : {}".format(precision))
    #print("recall : {}".format(recall))
    return [accuracy, precision, recall]
    



parser = ArgumentParser()
#parser.add_argument("-sequences", help = "path to the csv file containing sequences for the uniprots related to annotated structures")
parser.add_argument("-model_dataset", help = "path to the csv file containing annotated structures")
parser.add_argument("-path_results_to_write", help = "path where to save results")
parser.add_argument("-path_results_to_read", help = "path where to read results")
args = parser.parse_args()
print(args)
results = pd.read_csv(args.path_results_to_read, sep = ",")
results = pd.DataFrame(results)
print(results)
repeats_dataset = pd.read_csv(args.model_dataset)
repeats_dataset = pd.DataFrame(repeats_dataset)
#sequence_dataset = pd.read_csv(args.sequences, sep = "\t")
#sequence_dataset = pd.DataFrame(sequence_dataset)

list_pred = []

for index, res in results.iterrows() :
        pdb = res[0]
        print ("---PDB : {}".format(pdb))
        uniprot = repeats_dataset[(repeats_dataset["PDB"] == pdb[0:4]) & (repeats_dataset["chain"] == pdb[4])]["sp_primary"].values[0]
        print(uniprot)
        if isinstance(uniprot, float):
            print('No uniprot for this PDB')
            list_pred.append([res['PDB'], "No uniprot for this PDB", "No uniprot for this PDB", "No uniprot for this PDB", "No uniprot for this PDB", "No uniprot for this PDB"])
            continue
        if (pdb[0:4] in list(repeats_dataset["PDB"])) and (pdb[4] in list(repeats_dataset[(repeats_dataset["PDB"] == pdb[0:4]) & (repeats_dataset["chain"] == pdb[4])]["chain"])) and res[1] != "no regions" :
            repeats_dataset_filtered = repeats_dataset[(repeats_dataset["PDB"] == pdb[0:4]) & (repeats_dataset["chain"] == pdb[4])]
            found = False
            for index2, model in repeats_dataset_filtered.iterrows():
                #sequences = str(sequence_dataset.loc[sequence_dataset["Entry"] == uniprot, "Sequence"].values[0])

                seqlen = len(uniprot)
                units_predictions = res[2]
                #units_reference = repeats_dataset[(repeats_dataset["PDB"] == pdb[0:4]) & (repeats_dataset["chain"] == pdb[4])]["units"].values[0]
                units_reference = model[10]
                units_reference = literal_eval(units_reference)
                units_predictions = literal_eval(units_predictions)
                units_predictions = [[int(el) for el in sub] for sub in units_predictions]
                units_reference = [[int(el) for el in sub] for sub in units_reference]
                l1 = []
                l2 = []
                for units1 in units_reference:
                    v = set()
                    r1 = range(units1[0], (units1[1] + 1))
                    for el1 in list(r1) :
                        v.add(el1)
                    l1.append(v)

                for units2 in units_predictions:
                    v = set()
                    r2 = range(units2[0], (units2[1] + 1))
                    for el2 in list(r2) :
                        v.add(el2)
                    l2.append(v)
                #print(l1)
                print(set(range(units_predictions[0][0] + 1, units_predictions[-1][-1])))
                if (len(set(range(units_predictions[0][0] + 1, units_predictions[-1][-1])).intersection(set(range(units_reference[0][0] + 1, units_reference[-1][-1])))) > 1)  :
                    reg_reference = set(range(units_reference[0][0] + 1, units_reference[-1][-1]))
                    reg_predictions = set(range(units_predictions[0][0] + 1, units_predictions[-1][-1]))
                    TP_region = len(reg_reference.intersection(reg_predictions))
                    FP_region = len(reg_predictions.difference(reg_reference))
                    FN_region = len(reg_reference.difference(reg_predictions))
                    reg_accuracy = np.round(TP_region / (TP_region + FP_region + FN_region), 2)
                    reg_precision = np.round(TP_region / (TP_region + FP_region), 2)
                    reg_recall = np.round((TP_region) / (TP_region + FN_region), 2)

                    reference_class = repeats_dataset[(repeats_dataset["PDB"] == pdb[0:4]) & (repeats_dataset["chain"] == pdb[4])]["class"].values[0]
                    reference_topology = repeats_dataset[(repeats_dataset["PDB"] == pdb[0:4]) & (repeats_dataset["chain"] == pdb[4])]["topology"].values[0]
                    reference_fold = repeats_dataset[(repeats_dataset["PDB"] == pdb[0:4]) & (repeats_dataset["chain"] == pdb[4])]["fold"].values[0]
                    reference_clan = repeats_dataset[(repeats_dataset["PDB"] == pdb[0:4]) & (repeats_dataset["chain"] == pdb[4])]["clan"].values[0]
                    if (reference_class == res[4]) :
                        classe = 1
                    else :
                        classe = 0
                    if (reference_topology == res[5]) :
                        topology = 1
                    else :
                        topology = 0
                    if (reference_fold == res[6]) :
                        fold = 1
                    else :
                        fold = 0
                    if (reference_clan == res[7]) :
                        clan = 1
                    else :
                        clan = 0
                    dic_units = overlap_units(l1, l2)
                    list_eval = check_overlap(dic_units)
                    list_eval.insert(0, pdb)
                    list_eval.append(reg_accuracy)
                    list_eval.append(reg_precision)
                    list_eval.append(reg_recall)
                    list_eval.append(units_reference)
                    list_eval.append(units_predictions)
                    list_eval.append(classe)
                    list_eval.append(topology)
                    list_eval.append(fold)
                    list_eval.append(clan)
                    list_eval.append(str(res[8]))
                    list_pred.append(list_eval)
                    found = True
                    break
                else :
                    continue 

        if (found == False) :
            list_pred.append([pdb, "No matches", "No matches", "No matches", "No matches", "No matches", "No matches", "No matches", "No matches", "No matches", "No matches", "No matches", "No matches", res[8]])
print(list_pred)
metrics_results = pd.DataFrame(list_pred, columns = ["PDB", "unit_accuracy", "unit_precision", "unit_recall", "region_accuracy", "region_precision", "region_recall", "reference units", "predicted units", "classe", "topology", "fold", "clan", "time"])
print(metrics_results)
metrics_results.to_csv(args.path_results_to_write, sep = ',') 
