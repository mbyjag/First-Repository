from Bio import SeqIO
with open("file.fastq", "r") as data, open("100q.fastq", "w") as file100q, open("200q.fastq", "w") as file200q, open("500q.fastq", "w") as file500q,  open("100a.fasta", "w") as file100a, open("200a.fasta", "w") as file200a, open("500a.fasta", "w") as file500a:
    for record in SeqIO.parse(data, "fastq"):
        count = 100
        sequencesq = ("@" + record.title + "\n" + record.seq[:count] + "\n" + "+" + "\n" + record.letter_annotations[:count] + "\n")
        sequencesa = ("@" + record.title + "\n" + record.seq[:count])
        while count is 100 or 200 or 300:
            count += 100
            SeqIO.write(sequencesq, file100q, "fastq")
            SeqIO.write(sequencesa, file100a, "fasta")
