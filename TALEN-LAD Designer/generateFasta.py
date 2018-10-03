
import sys
import random

# consensus taken from Figure 5 of:

# DNA sequence-dependent compartmentalization and silencing of chromatin at the nuclear lamina.
# Zullo JM, Demarco IA, Pique-Regi R, Gaffney DJ, Epstein CB, Spooner CJ, Luperchio TR, Bernstein BE, Pritchard JK, Reddy KL, Singh H.

# Weights applied by guess from size of characters

consensus = [
             'a'*95+'t'*5,
             'g'*95+'c'*5,
             'a'*100,
             'g'*80+'t'*20,
             'a'*50+'g'*50,
             'g'*80+'t'*20,
             'a'*50+'g'*50,
             'c'*50+'g'*50,
             'a'*100,
             'g'*50+'c'*25+'a'*25,
             'a'*60+'t'*40,
             'g'*100,
             'a'*60+'g'*40,
             'a'*50+'g'*50,
             'a'*50+'c'*25+'t'*25,
             'g'*50+'t'*25+'a'*25,
             'a'*50+'t'*25+'c'*25,
             'a'*25+'c'*25+'t'*25+'g'*25,
             'a'*60+'g'*40,
             'g'*80+'t'*20,
             'g'*50+'c'*25+'a'*25,
             'g'*80+'c'*20,
             'a'*80+'g'*20,
             'a'*70+'t'*30,
             'g'*50+'t'*25+'a'*25,
             'a'*80+'g'*20,
             'g'*80+'t'*20,
             'a'*50+'c'*25+'t'*25,
             'g'*80+'a'*20,
             'a'*50+'t'*25+'g'*25,
             'g'*80+'a'*20,
             'a'*50+'g'*25+'t'*25,
             'a'*40+'g'*30+'c'*30,
             'a'*80+'g'*20,
             'g'*100
             ]

seqs = set()

for i in range(1000000):
    seq = ""
    for base in consensus:
        seq = seq + random.choice(base)
    
    if i % 100000 == 0:
        print >> sys.stderr, i
    
    for j in range(5):
        
        length = random.randrange(13,len(seq))
        
        start = random.randrange(len(seq)-length)
        
        seqs.add(seq[start:start+length])
    
for seq in seqs:
    print ">"+seq
    print seq