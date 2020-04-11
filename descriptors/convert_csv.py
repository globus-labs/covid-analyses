import pickle
import csv
import sys

pickle_list = pickle.load(open(sys.argv[1], 'rb'))

csv_results = []

count = 0

for s, d in pickle_list.items():
    
    if len(d[0]) == 0 or d[0] is None or d[0][0] == '':
        data =  [s, 'tmp_id%s' % count] + list(d[1])
    else: 
        data = [s, ''.join(d[0])] + list(d[1])
    data = ['' if str(i) == 'nan' else i for i in data]
    csv_results.append(data)
    count += 1

with open(sys.argv[2], 'w', newline='') as o_file:
    writer = csv.writer(o_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    writer.writerows(csv_results)
