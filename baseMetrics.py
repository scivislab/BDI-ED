import numpy as np
import math

def cost_wasserstein_split(n1,n2):
    if(n2==None):
        if(n1[0]!=3):
            return 0.0
        b1,d1 = n1[1]
        b2 = (b1+d1)*0.5
        d2 = (b1+d1)*0.5
        return math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    if(n1==None):
        if(n2[0]!=3):
            return 0.0
        b1,d1 = n2[1]
        b2 = (b1+d1)*0.5
        d2 = (b1+d1)*0.5
        return math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    if((n1[0]==3) and (n2[0]!=3)):
        return float("inf")
    if((n2[0]==3) and (n1[0]!=3)):
        return float("inf")
    if((n1[0]!=3) and (n2[0]!=3)):
        return 0.0
    b1,d1 = n1[1]
    b2,d2 = n2[1]
    d = math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    return d
def cost_wasserstein_split_squared(n1,n2):
    return cost_wasserstein_split(n1,n2)*cost_wasserstein_split(n1,n2)
    
def cost_wasserstein_join(n1,n2):
    if(n2==None):
        if(n1[0]!=0):
            return 0.0
        b1,d1 = n1[1]
        b2 = (b1+d1)*0.5
        d2 = (b1+d1)*0.5
        return math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    if(n1==None):
        if(n2[0]!=0):
            return 0.0
        b1,d1 = n2[1]
        b2 = (b1+d1)*0.5
        d2 = (b1+d1)*0.5
        return math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    if((n1[0]==0) and (n2[0]!=0)):
        return float("inf")
    if((n2[0]==0) and (n1[0]!=0)):
        return float("inf")
    if((n1[0]!=0) and (n2[0]!=0)):
        return 0.0
    b1,d1 = n1[1]
    b2,d2 = n2[1]
    d = math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    return d
def cost_wasserstein_join_squared(n1,n2):
    return cost_wasserstein_join(n1,n2)*cost_wasserstein_join(n1,n2)
    
def cost_wasserstein_branch(n1,p1,n2,p2):
    if(n1==None):
        b1 = n2
        d1 = p2
        b2 = (b1+d1)*0.5
        d2 = (b1+d1)*0.5
        return math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    if(n2==None):
        b1 = n1
        d1 = p1
        b2 = (b1+d1)*0.5
        d2 = (b1+d1)*0.5
        return math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    b1 = n1
    d1 = p1
    b2 = n2
    d2 = p2
    d = math.sqrt(abs(b1-b2)*abs(b1-b2)+abs(d1-d2)*abs(d1-d2))
    return d
def cost_wasserstein_branch_squared(n1,p1,n2,p2):
    return cost_wasserstein_branch(n1,p1,n2,p2)*cost_wasserstein_branch(n1,p1,n2,p2)
    
def cost_Linf_split(n1,n2):
    if(n2==None):
        if(n1[0]!=3):
            return 0.0
        b1,d1 = n1[1]
        return abs(d1-b1)*0.5
    if(n1==None):
        if(n2[0]!=3):
            return 0.0
        b1,d1 = n2[1]
        return abs(d1-b1)*0.5
    if((n1[0]==3) and (n2[0]!=3)):
        return float("inf")
    if((n2[0]==3) and (n1[0]!=3)):
        return float("inf")
    if((n1[0]!=3) and (n2[0]!=3)):
        return 0.0
    b1,d1 = n1[1]
    b2,d2 = n2[1]
    d = min(max(abs(b1-b2),abs(d1-d2)),
            (abs(d1-b1)+abs(d2-b2))*0.5)
    return d
def cost_Linf_split_squared(n1,n2):
    return cost_Linf_split(n1,n2)*cost_Linf_split(n1,n2)
    
def cost_Linf_join(n1,n2):
    if(n2==None):
        if(n1[0]!=0):
            return 0.0
        b1,d1 = n1[1]
        return abs(d1-b1)*0.5
    if(n1==None):
        if(n2[0]!=0):
            return 0.0
        b1,d1 = n2[1]
        return abs(d1-b1)*0.5
    if((n1[0]==0) and (n2[0]!=0)):
        return float("inf")
    if((n2[0]==0) and (n1[0]!=0)):
        return float("inf")
    if((n1[0]!=0) and (n2[0]!=0)):
        return 0.0
    b1,d1 = n1[1]
    b2,d2 = n2[1]
    d = min(max(abs(b1-b2),abs(d1-d2)),
            (abs(d1-b1)+abs(d2-b2))*0.5)
    return d
def cost_Linf_join_squared(n1,n2):
    return cost_Linf_join(n1,n2)*cost_Linf_join(n1,n2)
    
def cost_Linf_branch(n1,p1,n2,p2):
    if(n1==None):
        b1 = n2
        d1 = p2
        return abs(d1-b1)*0.5
    if(n2==None):
        b1 = n1
        d1 = p1
        return abs(d1-b1)*0.5
    b1 = n1
    d1 = p1
    b2 = n2
    d2 = p2
    d = min(max(abs(b1-b2),abs(d1-d2)),
            (abs(d1-b1)+abs(d2-b2))*0.5)
    return d
def cost_Linf_branch_squared(n1,p1,n2,p2):
    return cost_Linf_branch(n1,p1,n2,p2)*cost_Linf_branch(n1,p1,n2,p2)
    
def cost_overhang_split(n1,n2):
    if(n2==None):
        if(n1[0]!=3):
            return 0.0
        b1,d1 = n1[1]
        return abs(d1-b1)
    if(n1==None):
        if(n2[0]!=3):
            return 0.0
        b1,d1 = n2[1]
        return abs(d1-b1)
    if((n1[0]==3) and (n2[0]!=3)):
        return float("inf")
    if((n2[0]==3) and (n1[0]!=3)):
        return float("inf")
    if((n1[0]!=3) and (n2[0]!=3)):
        return 0.0
    b1,d1 = n1[1]
    b2,d2 = n2[1]
    d = min(abs(b1-b2)+abs(d1-d2),
            abs(b1-d1)+abs(b2-d2))
    return d
def cost_overhang_split_squared(n1,n2):
    return cost_overhang_split(n1,n2)*cost_overhang_split(n1,n2)
    
def cost_overhang_join(n1,n2):
    if(n2==None):
        if(n1[0]!=0):
            return 0.0
        b1,d1 = n1[1]
        return abs(d1-b1)
    if(n1==None):
        if(n2[0]!=0):
            return 0.0
        b1,d1 = n2[1]
        return abs(d1-b1)
    if((n1[0]==0) and (n2[0]!=0)):
        return float("inf")
    if((n2[0]==0) and (n1[0]!=0)):
        return float("inf")
    if((n1[0]!=0) and (n2[0]!=0)):
        return 0.0
    b1,d1 = n1[1]
    b2,d2 = n2[1]
    d = min(abs(b1-b2)+abs(d1-d2),
            abs(b1-d1)+abs(b2-d2))
    return d
def cost_overhang_join_squared(n1,n2):
    return cost_overhang_join(n1,n2)*cost_overhang_join(n1,n2)
    
def cost_overhang_branch(n1,p1,n2,p2):
    if(n1==None):
        b1 = n2
        d1 = p2
        return abs(d1-b1)
    if(n2==None):
        b1 = n1
        d1 = p1
        return abs(d1-b1)
    b1 = n1
    d1 = p1
    b2 = n2
    d2 = p2
    d = min(abs(b1-b2)+abs(d1-d2),
            abs(b1-d1)+abs(b2-d2))
    return d
def cost_overhang_branch_squared(n1,p1,n2,p2):
    return cost_overhang_branch(n1,p1,n2,p2)*cost_overhang_branch(n1,p1,n2,p2)
   
