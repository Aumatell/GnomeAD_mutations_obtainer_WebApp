#/usr/bin/env python3

# 3.7.4 (default, Aug  9 2019, 18:34:13) [MSC v.1915 64 bit (AMD64)]
# Running on winsows 10

from flask import Flask, render_template, request
import os
from Bio import SeqIO
from google.cloud import bigquery
import pandas as pd
import sys
from bs4 import BeautifulSoup
import re
from urllib.request import urlopen
import requests
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mc
##########################################################
# Functions
##########################################################
def mysqlform(gene_id_seq):
    '''
    It returns a MYSQL code to search fot her introduced proteins
    :param gene_id_seq:
    :return: Mysql query ID part
    '''
    gene_of_interest_formated = ""
    for i in range(len(gene_id_seq)):
        gene_of_interest_formated = gene_of_interest_formated + "vep.SWISSPROT = '" + str(gene_id_seq[i]) + "' OR "
    return gene_of_interest_formated.replace("\\", "")[:-3]

def selector(gene_id_seq, Transmemdic, out):
    '''
    Selects the elements inside or outside the dictionary
    :param gene_id_seq: Protein ID
    :param Transmemdic: Dictionary of transmem proteins
    :param out: True = not matching False = matching
    :return: It out == true returns the elements not in the dic , if its False reutns elements of the dic
    '''
    if out == True:
        return [ID for ID in gene_id_seq if ID not in Transmemdic.keys()]
    else:
        return [ID for ID in gene_id_seq if ID in Transmemdic.keys()]


##########################################################
# GENERATE APP
##########################################################
app = Flask(__name__)

@app.route('/')
def main():
    return render_template('transform.html')
##########################################################
# INPUT GENES
##########################################################

@app.route('/transform_results', methods=['POST'])
def transform_results():
    gene_id_seq = str(request.form['input_seq'])
    sequences_doc = False
    mutations_doc = False
    plots_doc = False
    sequences_doc = str(request.form['sequences_doc'])
    mutations_doc = str(request.form['mutations_doc'])
    plots_pergene = str(request.form['per_gene'])
    plots_type1 = str(request.form['plots_type1'])
    plots_type2 = str(request.form['plots_type2'])
    plots_type3 = str(request.form['plots_type3'])
    plots_type4 = str(request.form['plots_type4'])
    plots_type5 = str(request.form['plots_type5'])
    plots_type6 = str(request.form['plots_type6'])

    gene_id_seq = gene_id_seq.replace("\n", "").upper().replace(" ", "")
    gene_id_seq = list(gene_id_seq.split(","))

    ##########################################################
    # UNIPROT
    ##########################################################
    # TRANSMEMBRAIN PROTEIN REVIEWED IDENTIFICATION AND SEQUENCE
    response = requests.get(
        "https://www.uniprot.org/uniprot/?query=annotation:(type:transmem)&fil=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score&format=fasta")
    Transmemdic = {}
    first = True
    for line in response.text.split("\n"):
        if ">" in line:
            if first:
                first = False
            else:
                Transmemdic.update({Entry_key: [Entry_id, Entry_seq]})
            Entry_key = line[4:10]
            Entry_id = line
            Entry_seq = ""
        else:
            Entry_seq = Entry_seq + line

    ##########################################################
    # CHECK IF ALL OF THE INTRODUCED PROTEINS ARE TRANSMEMBRAIN
    ##########################################################
    transmembrain = selector(gene_id_seq, Transmemdic, False)
    notransmembrain = selector(gene_id_seq, Transmemdic, True)

    ##########################################################
    # FASTA DOCUMENT
    ##########################################################
    if sequences_doc:
        file=open("FASTA.fasta",'a')
        sys.stdout = file
        for ID in transmembrain:
            print(str(Transmemdic.get(str(ID))[0]) +"\n"+ re.sub("(.{60})", "\\1\n", str(Transmemdic.get(str(ID))[1]), 0, re.DOTALL)+"\n")

        file.close()
    ##########################################################
    # GENOMAD
    ##########################################################

    # PREPARATION OF THE ID PART OF THE QUERY
    gene_of_interest_query = mysqlform(gene_id_seq)
    # Ubication of json file
    os.environ[
        "GOOGLE_APPLICATION_CREDENTIALS"
    ] = "pymut-90f6558127b3.json"  # c:/Users/Quim/PycharmProjects/pythonProject1/PYTHONIC_APP/pymut-90f6558127b3.json
    df = pd.DataFrame()
    # Query for every chromosome
    client = bigquery.Client()  # Start the BigQuery Client
    for i in (list(map(str, range(1, 22))) + ["X", "Y"]):
        QUERY = ("SELECT names AS ID, reference_name AS CHROM, vep.SWISSPROT AS SWID, start_position AS POS, end_position AS END_POS,"
                 " reference_bases AS REF_B, vep.Consequence AS Consequence, vep.IMPACT AS Impact, vep.SYMBOL AS Symbol,"
                 " vep.Gene AS Gene, vep.CANONICAL AS Canonical, vep.PHENO AS PHENOTYPE, vep.GENE_PHENO AS G_PHENO, "
                 "quality AS QUALITY, filter AS QC_filter, vep.feature AS Feature, vep.allele AS allele, alternate_bases.allele_type AS ALLELE_TYPE, "
                 "vep.BIOTYPE AS BIOTYPE, alternate_bases.AC AS AC, vep.EXON AS EXON, vep.INTRON AS INTRON, vep.Gene AS ENSMBL,"
                 " vep.STRAND AS STRAND, vep.Protein_position AS Protein_pos FROM `bigquery-public-data.gnomAD.v2_1_1_exomes__chr" + str(
            i) + "` AS main_table,"
                 " main_table.alternate_bases AS alternate_bases, alternate_bases.vep AS vep WHERE "
                 + str(gene_of_interest_query)
                 + "LIMIT 100;"
                 )

        query_job = client.query(QUERY)  # Start Query API Request
        try:
            query_result = query_job.result()  # Get Query Result
        except:
            pass
        dfn = query_result.to_dataframe()  # Save the Query to Dataframe
        df = pd.concat([df, dfn])
    ####################################################
    # MUTATIONS FILE
    ####################################################
    if mutations_doc:
        df.to_csv("mutations_GenomAD.csv", index=False)
    ##########################################################
    # RETURN TEXT
    ##########################################################

    if len(notransmembrain) > 0:
        output_text = str("The following proteins are not found on the uniprot reviewed sequences of transmembrain proteins of humans: \n" + ",".join(notransmembrain) )
    else:
        output_text = str("All the proteins introduced are found in Swissprot database")
    input_text = str("The proteins introduced are: \n" + ",".join(gene_id_seq))

# PLOTS

    if plots_type1 or plots_type2 or plots_type3 or plots_type4 or plots_type5 or plots_type6:
        if plots_pergene == "True":
            for i in gene_id_seq: #For each gene
                subset=df[df.SWID== str(i)]
                if plots_type1:
                    fig = plt.figure(figsize=(15,10))
                    column ="Impact"  #column name
                    Values=[]
                    Labels=[]
                    Positions=[]
                    Color = np.array([])
                    counter=float(0)
                    for category in subset[column].unique():
                        Values.append(np.count_nonzero(np.array(subset[column]==str(category))))
                        Labels.append(str(category)+"\n"+str(Values[int(counter)]))
                        Positions.append(counter)
                        counter += 1

                    script_dir = os.path.dirname(__file__) #DIRECTORY
                    results_dir = os.path.join(script_dir, str(i)+'/')

                    if not os.path.isdir(results_dir): #Chech if directory exists
                        os.makedirs(results_dir)

                    clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                             (0.7, "green"), (0.75, "blue"), (1, "blue")] # COLORS
                    rvb = mc.LinearSegmentedColormap.from_list("", clist)
                    nuberelements = int(len(Positions))
                    x = np.arange(nuberelements).astype(float)
                    ##PLOT
                    plt.bar(Positions, Values, color=rvb(x/nuberelements)) # Columns
                    plt.xlabel(column) # label x
                    plt.ylabel("counts") #label y
                    plt.xticks(Positions,list(subset[column].unique()), rotation = 45) # column labels
                    plt.savefig(results_dir+str(i)+str(column)+".jpg")  # Generate file
                    ##
                if plots_type2:
                    fig = plt.figure(figsize=(15, 10))
                    column = "Consequence"
                    Values = []
                    Labels = []
                    Positions = []
                    Color = np.array([])
                    counter = float(0)
                    for category in subset[column].unique():
                        Values.append(np.count_nonzero(np.array(subset[column] == str(category))))
                        Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                        Positions.append(counter)
                        counter += 1

                    script_dir = os.path.dirname(__file__)
                    results_dir = os.path.join(script_dir, str(i) + '/')

                    if not os.path.isdir(results_dir):
                        os.makedirs(results_dir)

                    clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                             (0.7, "green"), (0.75, "blue"), (1, "blue")]
                    rvb = mc.LinearSegmentedColormap.from_list("", clist)
                    nuberelements = int(len(Positions))
                    x = np.arange(nuberelements).astype(float)
                    plt.bar(Positions, Values, color=rvb(x / nuberelements))
                    plt.xlabel(column)
                    plt.ylabel("counts")
                    plt.xticks(Positions, list(subset[column].unique()), rotation = 45)
                    plt.savefig(results_dir + str(i) + str(column) + ".jpg")
                if plots_type3:
                    fig = plt.figure(figsize=(15, 10))
                    column = "Feature"
                    Values = []
                    Labels = []
                    Positions = []
                    Color = np.array([])
                    counter = float(0)
                    for category in subset[column].unique():
                        Values.append(np.count_nonzero(np.array(subset[column] == str(category))))
                        Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                        Positions.append(counter)
                        counter += 1

                    script_dir = os.path.dirname(__file__)
                    results_dir = os.path.join(script_dir, str(i) + '/')

                    if not os.path.isdir(results_dir):
                        os.makedirs(results_dir)

                    clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                             (0.7, "green"), (0.75, "blue"), (1, "blue")]
                    rvb = mc.LinearSegmentedColormap.from_list("", clist)
                    nuberelements = int(len(Positions))
                    x = np.arange(nuberelements).astype(float)
                    plt.bar(Positions, Values, color=rvb(x / nuberelements))
                    plt.xlabel(column)
                    plt.ylabel("counts")
                    plt.xticks(Positions, list(subset[column].unique()), rotation = 45)
                    plt.savefig(results_dir + str(i) + str(column) + ".jpg")
                if plots_type4:
                    fig = plt.figure(figsize=(15, 5))
                    column = "EXON"
                    Values = []
                    Labels = []
                    Positions = []
                    Color = np.array([])
                    counter = float(0)
                    for category in subset[column].unique():
                        Values.append(np.count_nonzero(np.array(subset[column] == str(category))))
                        Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                        Positions.append(counter)
                        counter += 1

                    script_dir = os.path.dirname(__file__)
                    results_dir = os.path.join(script_dir, str(i) + '/')

                    if not os.path.isdir(results_dir):
                        os.makedirs(results_dir)

                    clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                             (0.7, "green"), (0.75, "blue"), (1, "blue")]
                    rvb = mc.LinearSegmentedColormap.from_list("", clist)
                    nuberelements = int(len(Positions))
                    x = np.arange(nuberelements).astype(float)
                    plt.bar(Positions, Values, color=rvb(x / nuberelements))
                    plt.xlabel(column)
                    plt.ylabel("counts")
                    plt.xticks(Positions, list(subset[column].unique()), rotation = 45)
                    plt.savefig(results_dir + str(i) + str(column) + ".jpg")
                if plots_type5:
                    pass
                    fig = plt.figure(figsize=(15, 10))
                    column = "INTRON"
                    Values = []
                    Labels = []
                    Positions = []
                    Color = np.array([])
                    counter = float(0)
                    for category in subset[column].unique():
                        Values.append(np.count_nonzero(np.array(subset[column] == str(category))))
                        Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                        Positions.append(counter)
                        counter += 1

                    script_dir = os.path.dirname(__file__)
                    results_dir = os.path.join(script_dir, str(i) + '/')

                    if not os.path.isdir(results_dir):
                        os.makedirs(results_dir)

                    clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                             (0.7, "green"), (0.75, "blue"), (1, "blue")]
                    rvb = mc.LinearSegmentedColormap.from_list("", clist)
                    nuberelements = int(len(Positions))
                    x = np.arange(nuberelements).astype(float)
                    plt.bar(Positions, Values, color=rvb(x / nuberelements))
                    plt.xlabel(column)
                    plt.ylabel("counts")
                    plt.xticks(Positions, list(subset[column].unique()), rotation = 45)
                    plt.savefig(results_dir + str(i) + str(column) + ".jpg")
                if plots_type6:
                    pass
                    fig = plt.figure(figsize=(15, 10))
                    column = "Protein_pos"
                    Values = []
                    Labels = []
                    Positions = []
                    Color = np.array([])
                    counter = float(0)
                    for category in subset[column].unique():
                        Values.append(np.count_nonzero(np.array(subset[column] == str(category))))
                        Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                        Positions.append(counter)
                        counter += 1

                    script_dir = os.path.dirname(__file__)
                    results_dir = os.path.join(script_dir, str(i) + '/')

                    if not os.path.isdir(results_dir):
                        os.makedirs(results_dir)

                    clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                             (0.7, "green"), (0.75, "blue"), (1, "blue")]
                    rvb = mc.LinearSegmentedColormap.from_list("", clist)
                    nuberelements = int(len(Positions))
                    x = np.arange(nuberelements).astype(float)
                    plt.bar(Positions, Values, color=rvb(x / nuberelements))
                    plt.xlabel(column)
                    plt.ylabel("counts")
                    plt.xticks(Positions, list(subset[column].unique()), rotation = 45)
                    plt.savefig(results_dir + str(i) + str(column) + ".jpg")
        else:

            if plots_type1:
                fig = plt.figure(figsize=(15, 10))
                column = "Impact"
                Values = []
                Labels = []
                Positions = []
                Color = np.array([])
                counter = float(0)
                for category in df[column].unique():
                    Values.append(np.count_nonzero(np.array(df[column] == str(category))))
                    Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                    Positions.append(counter)
                    counter += 1

                script_dir = os.path.dirname(__file__)
                results_dir = os.path.join(script_dir, "ALL/")

                if not os.path.isdir(results_dir):
                    os.makedirs(results_dir)

                clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                         (0.7, "green"), (0.75, "blue"), (1, "blue")]
                rvb = mc.LinearSegmentedColormap.from_list("", clist)
                nuberelements = int(len(Positions))
                x = np.arange(nuberelements).astype(float)
                plt.bar(Positions, Values, color=rvb(x / nuberelements))
                plt.xlabel(column)
                plt.ylabel("counts")
                plt.xticks(Positions, list(df[column].unique()), rotation = 45)
                plt.savefig(results_dir + "ALL" + str(column) + ".jpg")
            if plots_type2:
                fig = plt.figure(figsize=(15, 5))
                column = "Consequence"
                Values = []
                Labels = []
                Positions = []
                Color = np.array([])
                counter = float(0)
                for category in df[column].unique():
                    Values.append(np.count_nonzero(np.array(df[column] == str(category))))
                    Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                    Positions.append(counter)
                    counter += 1

                script_dir = os.path.dirname(__file__)
                results_dir = os.path.join(script_dir, "ALL" + '/')

                if not os.path.isdir(results_dir):
                    os.makedirs(results_dir)

                clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                         (0.7, "green"), (0.75, "blue"), (1, "blue")]
                rvb = mc.LinearSegmentedColormap.from_list("", clist)
                nuberelements = int(len(Positions))
                x = np.arange(nuberelements).astype(float)
                plt.bar(Positions, Values, color=rvb(x / nuberelements))
                plt.xlabel(column)
                plt.ylabel("counts")
                plt.xticks(Positions, list(df[column].unique()),rotation = 45)
                plt.savefig(results_dir + "ALL" + str(column) + ".jpg")
            if plots_type3:
                fig = plt.figure(figsize=(15, 5))
                column = "Feature"
                Values = []
                Labels = []
                Positions = []
                Color = np.array([])
                counter = float(0)
                for category in df[column].unique():
                    Values.append(np.count_nonzero(np.array(df[column] == str(category))))
                    Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                    Positions.append(counter)
                    counter += 1

                script_dir = os.path.dirname(__file__)
                results_dir = os.path.join(script_dir, "ALL" + '/')

                if not os.path.isdir(results_dir):
                    os.makedirs(results_dir)

                clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                         (0.7, "green"), (0.75, "blue"), (1, "blue")]
                rvb = mc.LinearSegmentedColormap.from_list("", clist)
                nuberelements = int(len(Positions))
                x = np.arange(nuberelements).astype(float)
                plt.bar(Positions, Values, color=rvb(x / nuberelements))
                plt.xlabel(column)
                plt.ylabel("counts")
                plt.xticks(Positions, list(df[column].unique()),rotation = 45)
                plt.savefig(results_dir + "ALL" + str(column) + ".jpg")
            if plots_type4:
                fig = plt.figure(figsize=(15, 5))
                column = "EXON"
                Values = []
                Labels = []
                Positions = []
                Color = np.array([])
                counter = float(0)
                for category in df[column].unique():
                    Values.append(np.count_nonzero(np.array(df[column] == str(category))))
                    Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                    Positions.append(counter)
                    counter += 1

                script_dir = os.path.dirname(__file__)
                results_dir = os.path.join(script_dir, "ALL" + '/')

                if not os.path.isdir(results_dir):
                    os.makedirs(results_dir)

                clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                         (0.7, "green"), (0.75, "blue"), (1, "blue")]
                rvb = mc.LinearSegmentedColormap.from_list("", clist)
                nuberelements = int(len(Positions))
                x = np.arange(nuberelements).astype(float)
                plt.bar(Positions, Values, color=rvb(x / nuberelements))
                plt.xlabel(column)
                plt.ylabel("counts")
                plt.xticks(Positions, list(df[column].unique()), rotation = 45)
                plt.savefig(results_dir + "ALL" + str(column) + ".jpg")
            if plots_type5:
                pass
                fig = plt.figure(figsize=(15, 5))
                column = "INTRON"
                Values = []
                Labels = []
                Positions = []
                Color = np.array([])
                counter = float(0)
                for category in df[column].unique():
                    Values.append(np.count_nonzero(np.array(df[column] == str(category))))
                    Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                    Positions.append(counter)
                    counter += 1

                script_dir = os.path.dirname(__file__)
                results_dir = os.path.join(script_dir, "ALL" + '/')

                if not os.path.isdir(results_dir):
                    os.makedirs(results_dir)

                clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                         (0.7, "green"), (0.75, "blue"), (1, "blue")]
                rvb = mc.LinearSegmentedColormap.from_list("", clist)
                nuberelements = int(len(Positions))
                x = np.arange(nuberelements).astype(float)
                plt.bar(Positions, Values, color=rvb(x / nuberelements))
                plt.xlabel(column)
                plt.ylabel("counts")
                plt.xticks(Positions, list(df[column].unique()), rotation = 45)
                plt.savefig(results_dir + "ALL" + str(column) + ".jpg")
            if plots_type6:
                pass
                fig = plt.figure(figsize=(15, 5))
                column = "Protein_pos"
                Values = []
                Labels = []
                Positions = []
                Color = np.array([])
                counter = float(0)
                for category in df[column].unique():
                    Values.append(np.count_nonzero(np.array(df[column] == str(category))))
                    Labels.append(str(category) + "\n" + str(Values[int(counter)]))
                    Positions.append(counter)
                    counter += 1

                script_dir = os.path.dirname(__file__)
                results_dir = os.path.join(script_dir, "ALL" + '/')

                if not os.path.isdir(results_dir):
                    os.makedirs(results_dir)

                clist = [(0, "red"), (0.125, "red"), (0.25, "orange"), (0.5, "green"),
                         (0.7, "green"), (0.75, "blue"), (1, "blue")]
                rvb = mc.LinearSegmentedColormap.from_list("", clist)
                nuberelements = int(len(Positions))
                x = np.arange(nuberelements).astype(float)
                plt.bar(Positions, Values, color=rvb(x / nuberelements))
                plt.xlabel(column)
                plt.ylabel("counts")
                plt.xticks(Positions, list(df[column].unique()), rotation = 45)
                plt.savefig(results_dir + "ALL" + str(column) + ".jpg")
                
    return render_template('transform_results.html', **locals())




app.run(debug=True)