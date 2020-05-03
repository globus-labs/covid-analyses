import os
import pickle
p = "/projects/CVD_Research/datasets/release/descriptors/"

for d in os.listdir(p):
    f = os.listdir(os.path.join(p, d))[0]

    if f.endswith('.pkl'):
        print (os.path.join(p,d,f))
        o = pickle.load(open(os.path.join(p, d, f), 'rb'))
        k = list(o.keys())[0]
        print ("%s: %s" % (k, o[k][0]))
    else: 
        f2 = os.listdir(os.path.join(p, d, f))[0]
        print(os.path.join(p,d,f,f2))
        o = pickle.load(open(os.path.join(p, d, f, f2), 'rb'))
        k = list(o.keys())[0]
        print ("%s: %s" % (k, o[k][0]))
