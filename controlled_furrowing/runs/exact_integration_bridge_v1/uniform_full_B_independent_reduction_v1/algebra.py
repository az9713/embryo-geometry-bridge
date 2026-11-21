"""Exact interval endpoint algebra for independent uniform-packet reduction."""
from fractions import Fraction as Q
def interval(v):
    a,b=Q(v['lower_exact']),Q(v['upper_exact'])
    if a>b:raise ValueError('ordered exact interval required')
    return a,b
def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(k,a):
    k=Q(k);return(k*a[0],k*a[1])if k>=0 else(k*a[1],k*a[0])
def abs_upper(a):return max(map(abs,a))
def abs_lower(a):return Q(0)if a[0]<=0<=a[1]else min(map(abs,a))
def contains(a,b):return a[0]<=b[0]<=b[1]<=a[1]
def overlap(a,b):return a[0]<=b[1]and b[0]<=a[1]
def scatter(ids,force,upper,total_force,total_upper):
    n=len(ids)
    if len(set(ids))!=n or any(type(i)is not int or not 0<=i<len(total_force)for i in ids)or len(force)!=n:raise ValueError('valid unique local/global free map required')
    if set(upper)!={(i,j)for i in range(n)for j in range(i,n)}:raise ValueError('complete unique local symmetric upper triangle required')
    for j,value in zip(ids,force):total_force[j]=add(total_force[j],value)
    for(i,j),value in upper.items():
        key=tuple(sorted((ids[i],ids[j])));total_upper[key]=add(total_upper.get(key,(Q(0),Q(0))),value)
def h_rows(upper,n):
    rows=[Q(0)]*n
    for(i,j),v in upper.items():
        if not 0<=i<=j<n:raise ValueError('canonical upper matrix index required')
        rows[i]+=abs_upper(v)
        if i!=j:rows[j]+=abs_upper(v)
    return rows
def signed_matvec(R,vector):
    if any(len(row)!=len(vector)for row in R):raise ValueError('matching matrix/vector dimensions required')
    out=[]
    for row in R:
        value=(Q(0),Q(0))
        for a,b in zip(row,vector):value=add(value,scale(a,b))
        out.append(value)
    return out
def absolute_row_product(R,rows):
    if any(len(row)!=len(rows)for row in R):raise ValueError('matching matrix/row dimensions required')
    return[sum((abs(Q(a))*b for a,b in zip(row,rows)),Q(0))for row in R]
