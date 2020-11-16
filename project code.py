from Bio.SeqIO.QualityIO import FastqGeneralIterator
count = 0
with open("file.fastq", "r") as data:
    with open("100q.fastq", "w") as file100q:
        with open("200q.fastq", "w") as file200q:
            with open("500q.fastq", "w") as file500q:
                with open("100a.fasta", "w") as file100a:
                    with open("200a.fasta", "w") as file200a:
                        with open("500a.fasta", "w") as file500a:
                            for title, seq, qual in FastqGeneralIterator(data):
                                count += 1
                                file100q.write("@" + title + "\n" + seq[:100] + "\n" + "+" + "\n" + qual[:100] + "\n")
                                file200q.write("@" + title + "\n" + seq[:200] + "\n" + "+" + "\n" + qual[:200] + "\n")
                                file500q.write("@" + title + "\n" + seq[:500] + "\n" + "+" + "\n" + qual[:500] + "\n")
                                file100a.write("@" + title + "\n" + seq[:100] + "\n")
                                file200a.write("@" + title + "\n" + seq[:200] + "\n")
                                file500a.write("@" + title + "\n" + seq[:500] + "\n")
