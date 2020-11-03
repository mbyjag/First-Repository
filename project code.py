filedata = input('enter file name')
file_data = open(filedata, 'r')
file_100 = open("file100", 'w')
file_200 = open("file200", 'w')
file_500 = open("file500", 'w')
counter = 0
for line in file_data:
    counter += 1
    if (counter - 2) % 4 == 0:
        file_100.write(line[:100])
        file_200.write(line[:200])
        file_500.write(line[:500])


# maybe need a try/except?
# if there isnt 500chr will it still write what there is?
# do we need seperate try/excepts for each file?

# saw stuff about
with open('ERR4425679_1.fastq', 'r') as file_data:
    with open('file100', 'w') as file_100:
        with open('file200', 'w') as file_200:
            with open('file500', 'w') as file_500:
# instead of 'file_date = open' so it closes itself?
# saw about using x instead of w if file doesnt already exist?


file_data = open("ERR4425679_1.fastq", 'r')
file_100 = open("file100", 'w')
file_200 = open("file200", 'w')
file_500 = open("file500", 'w')
counter = 0
for line in file_data:
    counter += 1
    if (counter - 2) % 4 == 0:
        file_100.write(line[:100])
        file_200.write(line[:200])
        file_500.write(line[:500])