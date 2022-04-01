import csv
import needlemanBase
import sys

# Receives Input File from User
filename = str(sys.argv[1])

f = open(filename, 'r')   # Opens file path for designated .csv file
reader = csv.reader(f)

# Lists for storing data
align = []
score = []
seq_1 = []
seq_2 = []

# Adds data from input file to separate lists
for row in reader:
    try:
        seq_1.append(row[0])
        seq_2.append(row[1])
    except:
        pass

# Removes Sequence Headers (sequence1, sequence2) located on first row
seq_1.pop(0)
seq_2.pop(0)

# Adds (Aligned Results) and (Score) obtained from Needleman-Wunsch algorithm to their respective lists
for i in range(len(seq_1)):
    tempSave = needlemanBase.needle(str(seq_1[i]), str(seq_2[i]))
    align.append(tempSave.pop())
    score.append(tempSave.pop())

    # Test Prints for viewing results in terminal
    # print("Aligned Result is:", align[i])
    # print("Score is:", score[i])


# Opens a new .csv file for result output
p = open("result.csv", 'a', newline='')

# Adds Headers for columns to .csv file
row = ("sequence1", "sequence2", "alignment text", "alignment score")
writer = csv.writer(p)
writer.writerow(row)

# Adds the data stored in each list to the .csv file
for x in range(len(seq_1)):
    row = (str(seq_1[x]), str(seq_2[x]), str(align[x]), score[x])
    writer = csv.writer(p)
    writer.writerow(row)

# Closes .csv file
p.close()