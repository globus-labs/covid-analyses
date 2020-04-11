import pickle
import csv
import sys

pickle_list = pickle.load(open(sys.argv[1], 'rb'))

csv_results = []

for s, d in pickle_list.items():
    print (s)
    data = [s, ''.join(d[0])] + list(d[1])
    data = ['' if str(i) == 'nan' else i for i in data]
    csv_results.append(data)


with open(sys.argv[2], 'w', newline='') as o_file:
    writer = csv.writer(o_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    writer.writerows(csv_results)
