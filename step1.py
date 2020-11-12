#! /usr/bin/env python3
# author : Picolo Floriane

import argparse
import glob

def parser():
    parser = argparse.ArgumentParser(description = "ORTHOLOGY ANALYSIS")
    parser.add_argument("-f", "--fish", metavar = "PATH", type = str, required = True, help = "pathway orthologous files")
    parser.add_argument("-l", "--human", metavar = "PATH", type = str, required = True, help = "specific human list to analysis")
    parser.add_argument("-s", "--specificity", metavar = "STR", type = str, required = True, help = "specific name of gene list")
    return parser

def create_humanlist(pathway):
    """
        Prend un fichier brut de gènes humains spécifiques qu'on met dans une liste
    """
    with open(pathway) as in_file:
        liste = []
        for row in in_file:
            liste.append(row.strip()) # add line without "\n"
    liste = list(set(liste)) # delete double variable 
    return liste

def create_dicoHF(pathway):
    """ 
        Prend un fichier brut sortie de biomart comprenant les orthologues téléostéens des gènes humains qu'on met dans un dico
        Récupère le nom du téléostéen étudié
    """
    dico_HF = {} # key : human gene , values : fish gene
    with open(pathway) as infile: 
        for row in infile:
            if "Gene" in row:
                fish_name = row.strip().split(",")[1][:-15].replace(" ", "_") # a changer quand on aura nos propres fichiers
            else:
                line = row.strip().split(",")
                if line[1] != "": 
                    if line[0] not in dico_HF.keys():
                        dico_HF[line[0]] = line[1]
                    else:
                        dico_HF[line[0]] += ";" + line[1]
    return dico_HF, fish_name

def extract_specificgenes_to_generallist(Hliste, HFdico):
    """
        Extrait du dictionnaire général, les gènes humains étant dans la liste spécifique
        Sépare les gènes retrouvés en duplicat et en singleton
    """
    dicoHF_Dinlist = {} # key = human gene in specific list ; values : fish gene in dupli
    dicoHF_Sinlist = {} # key = human gene in specific list ; values : fish gene in singleton
    for element in Hliste: 
        if element in HFdico.keys():
            if ";" in HFdico[element]:
                dicoHF_Dinlist[element] = dicoHF[element]
            else:
                dicoHF_Sinlist[element] = dicoHF[element]
    return dicoHF_Dinlist, dicoHF_Sinlist

def general_count(HFdico):
    """
        Prend un dico général et compte le nombre de gène humain ayant un orthologue en dupli ou en singleton chez l'espèce choisi
    """
    countS = countD = 0
    for k, v in HFdico.items():
        if ";" in v: 
            countD += 1
        else: 
            countS += 1
    return countS, countD

def write_csvfile_forR(dicocount, spe_name):
    with open (spe_name + "_chi2.csv", "w") as outfile: 
        outfile.write("Species;Nb_single;Nb_dupli;Nb_single_in_list;Nb_dupli_in_list" + "\n")
        for k, v in dicocount.items(): 
            outfile.write(k + ";" + str(v[0]) + ";" + str(v[1]) + ";" + str(v[2]) + ";" + str(v[3]) + "\n")

def write_upSetfile(Hliste, spedico, spe_name, letter):
    """
        Mettre 0 ou 1 en fonction de si le gène humain est présent ou non chez l'espèce
    """
    with open(spe_name + "_" + letter +"_upSet.csv", "w") as outfile:
        outfile.write("Human_gene")
        for specie in spedico.keys():
            outfile.write(";" + specie)
        for humangene in Hliste:
            outfile.write("\n" + humangene)
            for specie in spedico.keys(): 
                if humangene in spedico[specie]:
                    outfile.write(";" + str("1"))
                else: 
                    outfile.write(";" + str("0"))

if __name__ == "__main__":
    
    parser = parser()
    args = parser.parse_args()

    # args
    pathway_fish = args.fish
    pathway_human = args.human
    specificity_name = args.specificity
    
    specificdico_D = {} # key = specie name ; values = dicoHF avec uniquement les gènes humains de la liste
    specificdico_S = {} # key = specie name ; values = dicoHF avec uniquement les gènes humains de la liste
    dico_count = {} # key = specie name ; values = count list [nb S in specie, nb D in specie, nb S in list, nb D in list]

    # programms
    humanlist = create_humanlist(pathway_human)

    for f in glob.glob(pathway_fish + "/*.txt"):
        dicoHF, specie = create_dicoHF(f)
        dicoHFspe_D, dicoHFspe_S = extract_specificgenes_to_generallist(humanlist, dicoHF)
        specificdico_D[specie] = dicoHFspe_D
        specificdico_S[specie] = dicoHFspe_S
        nbS, nbD = general_count(dicoHF)
        dico_count[specie] = [nbS, nbD, len(dicoHFspe_S), len(dicoHFspe_D)]
    
    write_csvfile_forR(dico_count, specificity_name)
    write_upSetfile(humanlist, specificdico_D, specificity_name, "D") 
    write_upSetfile(humanlist, specificdico_S, specificity_name, "S") 



