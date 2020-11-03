filename = "FAL00432_af307028575a34db3268b4df36412691cb25b558_0.fastq"
file_data = open(filename, "r")
file_100 = open("file100.txt", 'w')
file_200 = open("file200.txt", 'w')
file_500 = open("file500.txt", 'w')
counter = 0
for line in file_data:
    counter += 1
    if (counter - 2) % 4 == 0:
        file_100.write(line[:100]+"\n")
        file_200.write(line[:200]+"\n")
        file_500.write(line[:500]+"\n")

file_data.close()
file_100.close()
file_200.close()
file_500.close()
