total_dist = 250000000

with open("test_htrans_probs.txt",'w') as outfile:
    for dist in range(0,total_dist,50000):

        prob = dist/30000000

        if prob > 0.35:
            prob = 0.35

        print("{}\t{}".format(dist, prob),file = outfile)
