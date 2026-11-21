"""Synthetic tests only: cannot accept unfinished uniform producer results."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
from flint import arb,arb_mat,ctx,fmpq
import algebra as A
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(q):return arb(fmpq(q.numerator,q.denominator))
def main():
    checks={};checks['negative_endpoint_scaling']=A.scale(-2,(Q(-1),Q(3)))==(Q(-6),Q(2));checks['zero_crossing_abs_lower']=A.abs_lower((Q(-1),Q(3)))==0
    F=[(Q(0),Q(0))for _ in range(3)];H={};A.scatter([2,0],[(Q(1),Q(2)),(Q(-3),Q(-2))],{(0,0):(Q(-2),Q(1)),(0,1):(Q(3),Q(4)),(1,1):(Q(-5),Q(-4))},F,H)
    checks['unsorted_symmetric_scatter']=H=={(2,2):(Q(-2),Q(1)),(0,2):(Q(3),Q(4)),(0,0):(Q(-5),Q(-4))};checks['signed_force_scatter']=F==[(Q(-3),Q(-2)),(Q(0),Q(0)),(Q(1),Q(2))]
    rows=A.h_rows(H,3);checks['offdiag_counts_both_rows_once']=rows==[Q(9),Q(0),Q(6)]
    R=[[Q(1),Q(-2),Q(3)],[Q(-4),Q(1),Q(-2)],[Q(0),Q(1),Q(1)]];product=A.signed_matvec(R,F);checks['signed_RF_exact_endpoint_sum']=product==[(Q(0),Q(4)),(Q(4),Q(10)),(Q(1),Q(2))]
    major=A.absolute_row_product(R,rows);checks['absolute_R_times_Hrows_exact']=major==[Q(27),Q(48),Q(6)]
    ctx.prec=256;RM=arb_mat(3,3,[ball(x)for row in R for x in row]);HM=arb_mat(3,3)
    for(i,j),v in H.items():HM[i,j]=HM[j,i]=ball(v[0]).union(ball(v[1]))
    RH=RM*HM
    for i in range(3):
        for j in range(3):
            exact=(Q(0),Q(0))
            for k in range(3):exact=A.add(exact,A.scale(R[i][k],H.get(tuple(sorted((k,j))),(Q(0),Q(0)))))
            checks[f'signed_Arb_matrix_encloses_exact_{i}_{j}']=Q(str(RH[i,j].lower().fmpq()))<=exact[0]<=exact[1]<=Q(str(RH[i,j].upper().fmpq()))
    checks['signed_product_rows_below_absolute_R_majorants']=all(sum(Q(str(abs(RH[i,j]).abs_upper().fmpq()))for j in range(3))<=major[i]for i in range(3))
    try:A.scatter([0,0],[],{},F,H);checks['duplicate_map_rejected']=False
    except ValueError:checks['duplicate_map_rejected']=True
    try:A.interval({'lower_exact':'1','upper_exact':'0'});checks['reversed_interval_rejected']=False
    except ValueError:checks['reversed_interval_rejected']=True
    assert all(v is True for v in checks.values()),checks
    out={'status':'PASS_SYNTHETIC_REDUCER_ONLY','checks':checks,'source_sha256':sha(Path(__file__)),'algebra_sha256':sha(HERE/'algebra.py'),'actual_uniform_producer_results_accepted':False,'actual_integrals_or_physical_solves_run':False};p=HERE/'selfchecks.json';assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'checks':len(checks),'sha256':sha(p)}))
if __name__=='__main__':main()
