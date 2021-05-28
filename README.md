from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from ete3 import Tree, faces, TreeStyle
import hashlib
import csv
import os
import time
import codecs
import subprocess
import pyfastx

counter = 0
#input filepath and split it to get the filepath and filename(-extension)
filepath = input("Enter your Filepath (including filename):")
path = "/".join(filepath.split("/")[:-1]) + "/"
filename_exe = os.path.basename(filepath)
filename = os.path.splitext(filename_exe)[0]

cleanpath = os.path.join(path, "cleaned")
filteredpath = os.path.join(path, "filtered")
pangolinpath = os.path.join(path, "pangolin_output")
alignedpath = os.path.join(path, "aligned")
trimmedpath = os.path.join(path, "trimmed")
iqtreepath = os.path.join(path, "iqtree_output")
visualizedpath = os.path.join(path, "visualized")

#see if intermediate files exist. If they do move on, so you dont make files again.
if not os.path.exists(cleanpath):
    os.makedirs(cleanpath)
if not os.path.exists(filteredpath):
    os.makedirs(filteredpath)
if not os.path.exists(pangolinpath):
    os.makedirs(pangolinpath)
if not os.path.exists(alignedpath):
    os.makedirs(alignedpath)
if not os.path.exists(trimmedpath):
    os.makedirs(trimmedpath)
if not os.path.exists(iqtreepath):
    os.makedirs(iqtreepath)
if not os.path.exists(visualizedpath):
    os.makedirs(visualizedpath)

def cleanstep():
    ### name cleaner ###
    #making the sequence header just the seqid
    num = 0
    with open(filepath) as data:
        records = list(SeqIO.parse(data, "fasta"))
        number_of_elements = len(records)
        wf = open(cleanpath + "/" + "c_" + filename + "_" + str(counter) + ".fasta", "w")
        data1 = open(filepath, "r")
        while num < number_of_elements:
            for line in data1:
                if line.startswith(">"):
                    wf.write(">" + records[num].id + "\n")
                    num += 1
                else:
                    wf.write(line)
    wf.close()

    #changing any / to _ in header (forwardslashes in headers mess with the tree builder as it reads it as a filepath)
    with open(cleanpath + "/" + "cleaned_" + filename + "_" + str(counter) + ".fasta", "w+") as cleanfile:
        name = codecs.open(cleanpath + "/" + "c_" + filename + "_" + str(counter) + ".fasta").read()
        clean_name = name.replace("/", "_")
        cleanfile.write(clean_name)
    #removing the 1st round cleaning file once the fully cleaned file is made.
    if os.path.exists(cleanpath + "/" + "cleaned_" + filename + "_" + str(counter) + ".fasta"):
        os.remove(cleanpath + "/" + "c_" + filename + "_" + str(counter) + ".fasta")

    ### filter ###
    def filter(fastafile, output1, output2):
        filtered50 = 0
        filtered90 = 0
        num = 0
        g = 0; c = 0; a = 0; t = 0
        unknown = 0
        seq_list = []
        seq_id = None
        seq_str = None

        #once it reaches a ">" in the file it will print the id and seq if base content is higher than threshold
        for line in fastafile:
            if line.startswith(">"):
                num += 1
                if (g + c + a + t) > 0:
                    total = g + c + a + t
                    percentage = round((total / (total + unknown)) * 100)
                    if percentage > 90:
                        output1.write(seq_id + "\n"); output1.write(seq_str + "\n")
                        output2.write(seq_id + "\n"); output2.write(seq_str + "\n")
                        filtered90 += 1
                        filtered50 += 1
                    elif percentage > 50:
                        output2.write(seq_id + "\n"); output2.write(seq_str + "\n")
                        filtered50 += 1
                    seq_id = line.strip()

                #reset count
                g = 0; c = 0; a = 0; t = 0
                unknown = 0
                seq_list = []
            #it starts counting the number of bases after the first ">"
            #this means that it only prints the first sequence id and seq when it reaches the 2nd ">"
            else:
                nuc_str = list(line.strip())
                seq_list.append(nuc_str)
                seq_list = [item for elem in seq_list for item in elem]
                seq_str = "".join(seq_list)
                for n in nuc_str:
                    if n == "G" or n == "g":
                        g += 1
                    elif n == "C" or n == "c":
                        c += 1
                    elif n == "A" or n == "a":
                        a += 1
                    elif n == "T" or n == "t":
                        t += 1
                    else:
                        unknown += 1

        #it needs a count for the last sequence in file
        total = g + c + a + t + unknown
        percentage = round((total / (total + unknown)) * 100)
        if percentage > 90:
            filtered90 += 1
            filtered50 += 1
            output1.write(seq_id + "\n"); output1.write(seq_str + "\n")
            output2.write(seq_id + "\n"); output2.write(seq_str + "\n")
        elif percentage > 50:
            filtered50 += 1
            output2.write(seq_id + "\n"); output2.write(seq_str + "\n")

        print("filter90 (f90) represents " + str(filtered90) + "/" + str(num) + " sequences")
        print("filter50 (f50) represents " + str(filtered50) + "/" + str(num) + " sequences")

        output1.close()
        output2.close()

    #if there is already a filtered file, then take its contents and make it the base of a new filtered file, then append new contents to this file
    if os.path.exists(filteredpath + "/" + "f90_" + filename + "_" + str(counter-1) + ".fasta"):
        fastafile = open(cleanpath + "/" + "cleaned_" + filename + "_" + str(counter) + ".fasta", "r")
        oldseqs1 = open(filteredpath + "/" + "f90_" + filename + "_" + str(counter-1) + ".fasta", "r")
        oldseqs2 = open(filteredpath + "/" + "f50_" + filename + "_" + str(counter-1) + ".fasta", "r")
        output1 = open(filteredpath + "/" + "f90_" + filename + "_" + str(counter) + ".fasta", "a")
        output2 = open(filteredpath + "/" + "f50_" + filename + "_" + str(counter) + ".fasta", "a")
        for line in oldseqs1:
            output1.write(line)
        for line in oldseqs2:
            output2.write(line)

        filter(fastafile, output1, output2)

    else:
        fastafile = open(cleanpath + "/" + "cleaned_" + filename + "_" + str(counter) + ".fasta", "r")
        output1 = open(filteredpath + "/" + "f90_" + filename + "_" + str(counter) + ".fasta", "w")
        output2 = open(filteredpath + "/" + "f50_" + filename + "_" + str(counter) + ".fasta", "w")

        filter(fastafile, output1, output2)

    #backingup uncleaned, unfiltered sequences to backup file
    with open(filepath) as f:
        with open(path + "/" + "backup_seqs_" + filename, "a") as f1:
            for line in f:
                f1.write(line)   # I am using "a" so "write" will append to the file

    #empyting datafile for new sequences
    empty = open(filepath, "r+")
    empty.truncate(0)
    empty.close()

def buildstep():
    ### getting pangolin lineages for the seqs ###
    subprocess.run(["pangolin", filteredpath + "/" + "f90_" + filename + "_" + str(counter) + ".fasta", "--outdir", pangolinpath, "--outfile", "f90_" + filename + "_" + str(counter) + ".csv"])
    subprocess.run(["pangolin", filteredpath + "/" + "f50_" + filename + "_" + str(counter) + ".fasta", "--outdir", pangolinpath, "--outfile", "f50_" + filename + "_" + str(counter) + ".csv"])

    def align(filtereddata, alignedoutfile):
        ### alignment ###
        mafft_exe = "/usr/bin/mafft"
        reference = path + "MN908947.3.fasta"
        mafft_cline = subprocess.run([mafft_exe, "--6merpair", "--addfragments", filtereddata, reference, ">", alignedoutfile])
        stdout, stderr = mafft_cline
        with open(alignedoutfile, "w") as align:
            align.write(stdout)

    align(filteredpath + "/" + "f90_" + filename + "_" + str(counter) + ".fasta", alignedpath + "/" + "aligned_f90_" + filename + "_" + str(counter) + ".fasta")
    align(filteredpath + "/" + "f50_" + filename + "_" + str(counter) + ".fasta", alignedpath + "/" + "aligned_f50_" + filename + "_" + str(counter) + ".fasta")

    def triming(trimpath, trimoutput):
        ### trimming aligned seqs ###
        fa = pyfastx.Fasta(trimpath, full_name=True)
        list = []
        num = 0
        with open(trimpath, "r") as data:
            trimmed = open(trimoutput, "w+")
            for line in data:
                if line.startswith(">"):
                    line = line.rstrip("\n")
                    line = line.strip(">")
                    list.append(line)

            number_of_elements = len(list)
            while num < number_of_elements:
                trimmed.write(">" + list[num] + "\n")
                s1 = fa[list[num]]
                s1 = str(s1)
                trimmed.write(s1[129:-50] + "\n")
                num += 1

    triming(alignedpath + "/" + "aligned_f90_" + filename + "_" + str(counter) + ".fasta", trimmedpath + "/" + "t90_" + filename + "_" + str(counter) + ".fasta")
    triming(alignedpath + "/" + "aligned_f50_" + filename + "_" + str(counter) + ".fasta", trimmedpath + "/" + "t50_" + filename + "_" + str(counter) + ".fasta")

    ### tree making ###
    subprocess.run(["/media/sf_4th_year_project_files/iqtree-2.1.2-Linux/bin/iqtree2", "-s", alignedpath + "/" + "aligned_f90_" + filename + "_" + str(counter) + ".fasta", "-keep-ident", "-m", "GTR+R", "--prefix", iqtreepath + "/" + "90_" + filename + "_" + str(counter)])
    subprocess.run(["/media/sf_4th_year_project_files/iqtree-2.1.2-Linux/bin/iqtree2", "-s", alignedpath + "/" + "aligned_f50_" + filename + "_" + str(counter) + ".fasta", "-keep-ident", "-m", "GTR+R", "--prefix", iqtreepath + "/" + "50_" + filename + "_" + str(counter)])

    def visualizingtree(pangopath, treepath, pngoutput):
        ### visualising tree as png###
        #reading lineages into dictionary
        with open(pangopath, "r") as csvfile:
            reader = csv.reader(csvfile)
            lineage = {rows[0]: rows[1] for rows in reader}

        #turning newick tree into ete tree object
        nw = treepath
        t = Tree(nw, format=1, quoted_node_names=True)
        #formating
        nameface = faces.AttrFace("name", fsize=20, fgcolor="#000000")

        def mylayout(node):
            #formating
            if node.is_leaf():
                faces.add_face_to_node(nameface, node, column=0)
                lineageface = faces.TextFace(lineage[node.name], fsize=20)
                faces.add_face_to_node(lineageface, node, column=0)
                node.img_style["size"] = 12
                node.img_style["shape"] = "circle"
            else:
                node.img_style["size"] = 3
                node.img_style["shape"] = "circle"
                node.img_style["fgcolor"] = "#000000"

        ts = TreeStyle()
        ts.layout_fn = mylayout
        ts.show_leaf_name = False
        ts.show_branch_length = True
        ts.branch_vertical_margin = 20
        ts.scale = 500
        t.render(pngoutput, w=1000, tree_style=ts)

    visualizingtree(pangolinpath + "/" + "f90_" + filename + "_" + str(counter) + ".csv", iqtreepath + "/" + "90_" + filename + "_" + str(counter) + ".treefile", visualizedpath + "/" + "90_" + filename + "tree" + "_" + str(counter) + ".png")
    visualizingtree(pangolinpath + "/" + "f50_" + filename + "_" + str(counter) + ".csv", iqtreepath + "/" + "50_" + filename + "_" + str(counter) + ".treefile", visualizedpath + "/" + "50_" + filename + "tree" + "_" + str(counter) + ".png")


def on_modified():
    print(f"new data found")

def sha1(filename):
    BUF_SIZE = 65536
    sha1 = hashlib.sha1()
    with open(filename, "rb") as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            sha1.update(data)
    return sha1.hexdigest()

def check_for_change():
    #get hashes
    hash_new = sha1(filepath)
    #compare hashes
    if hash_current != hash_new:
        on_modified()
        global counter
        counter += 1
        cleanstep()
        buildstep()
    else:
        print("no new data")

#running cleaning and building on initial data
cleanstep()
buildstep()

global hash_current
hash_current = sha1(filepath)

#check every 30 minutes for new data. if new data found then run code on it.
try:
    while True:
        print("waiting for more data")
        check_for_change()
        time.sleep(1800)
except KeyboardInterrupt:
    print("finished")

