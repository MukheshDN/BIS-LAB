# Parallel Cellular Algorithm – HP protein folding (NO GRAPH VERSION)

import random
import math
import copy

def neighbors(p):
    x,y = p
    return [(x+1,y),(x-1,y),(x,y+1),(x,y-1)]

class Protein:
    def __init__(self,seq):
        self.seq = seq.upper()
        self.n = len(seq)
        self.coords=[(i,0) for i in range(self.n)]
        self.occ=set(self.coords)

    def clone(self):
        p=Protein(self.seq)
        p.coords=list(self.coords)
        p.occ=set(self.coords)
        return p

    def set_coords(self,c):
        self.coords=c
        self.occ=set(c)

    def energy(self):
        pos2i={p:i for i,p in enumerate(self.coords)}
        E=0
        for i,p in enumerate(self.coords):
            if self.seq[i]=='H':
                for nb in neighbors(p):
                    j=pos2i.get(nb,None)
                    if j!=None and abs(i-j)>1 and self.seq[j]=='H':
                        if i<j:
                            E-=1
        return E

    def pivot(self):
        k=random.randint(1,self.n-2)
        d=random.choice([1,-1])
        px,py=self.coords[k]
        nc=self.coords[:k+1]
        for (x,y) in self.coords[k+1:]:
            dx,dy=x-px,y-py
            if d==1: rx,ry=-dy,dx
            else:    rx,ry=dy,-dx
            nc.append((px+rx,py+ry))
        if len(nc)==len(set(nc)): return nc
        return None

    def end(self):
        for end in [0,self.n-1]:
            mid=1 if end==0 else self.n-2
            for nb in neighbors(self.coords[mid]):
                if nb not in self.occ and nb!=self.coords[end]:
                    if end==0: nc=[nb]+self.coords[1:]
                    else: nc=self.coords[:-1]+[nb]
                    if len(nc)==len(set(nc)): return nc
        return None

def propose(p):
    if random.random()<0.25: return p.end()
    return p.pivot()

def accept(dE,T):
    if dE<=0: return True
    return random.random()<math.exp(-dE/T)

def parallel_iter(p,T,par=8):
    E0=p.energy()
    cands=[]
    for _ in range(par):
        nc=propose(p)
        if nc!=None:
            tmp=p.clone()
            tmp.set_coords(nc)
            cands.append((nc,tmp.energy()))
    random.shuffle(cands)
    for nc,E in cands:
        if accept(E-E0,T):
            p.set_coords(nc)
            return True,E
    return False,E0

def run(seq,iters=2000,temp=1.2,dec=0.999,par=10):
    p=Protein(seq)
    best=p.clone();bestE=best.energy()
    T=temp
    for _ in range(iters):
        ok,E=parallel_iter(p,T,par)
        if ok and E<bestE:
            best=p.clone();bestE=E
        T*=dec
    return bestE, p.coords

# ---------------- RUN HERE ------------------

seq="HHPHPPHHPH"  # <<< YOUR INPUT SEQUENCE
bestE,coords = run(seq)

print("Sequence:",seq)
print("BEST ENERGY FOUND:",bestE)
print("BEST COORDINATES:",coords)
