"""Exact small-data sensitivity to a HYPOTHETICAL exception budget.

All actual mechanical calibrations retain all 17 records. This diagnostic
does not declare any embryo an outlier or estimate measurement noise.
"""
from fractions import Fraction
from decimal import Decimal
from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def exact(value):
    return Fraction(Decimal(value))


def maximum_compatible_count(y,twice_delta):
    """Independent interval-subsequence DP; all values use doubled units.

    tails[k] is the smallest feasible final fitted value after retaining k
    records from the processed prefix. A smaller final value dominates a
    larger one for every future nondecreasing interval choice.
    """
    if any(value<0 for value in y) or twice_delta<0:
        raise ValueError('Nonnegative response and discrepancy units required')
    tails=[None]*(len(y)+1);tails[0]=0
    for value in y:
        lower=max(0,2*value-twice_delta);upper=2*value+twice_delta
        for count in range(len(y)-1,-1,-1):
            if tails[count] is not None and tails[count]<=upper:
                candidate=max(lower,tails[count])
                old=tails[count+1]
                tails[count+1]=candidate if old is None else min(old,candidate)
    return max(k for k,value in enumerate(tails) if value is not None)


def enumerate_subsets(y):
    if any(value<0 for value in y):
        raise ValueError('Nonnegative response units required')
    n=len(y);cost=np.zeros(1<<n,dtype=np.int64)
    best=[None]*(n+1);masks=[None]*(n+1)
    best[0]=0;masks[0]=0
    for mask in range(1,1<<n):
        last=mask.bit_length()-1;previous=mask^(1<<last)
        value=int(cost[previous]);remaining=previous
        while remaining:
            bit=remaining&-remaining;earlier=bit.bit_length()-1
            value=max(value,y[earlier]-y[last]);remaining^=bit
        cost[mask]=value
        size=mask.bit_count()
        if best[size] is None or value<best[size]:
            best[size]=value;masks[size]=mask
    return cost,best,masks


def main():
    source=ROOT/'data/derived/calibration_small_18um.csv'
    with source.open(encoding='utf-8-sig',newline='') as stream:
        rows=sorted(csv.DictReader(stream),key=lambda row:exact(row['activation_magnitude']))
    dose=[exact(row['activation_magnitude']) for row in rows]
    if not all(value>0 for value in dose):
        raise ValueError('Strictly positive reported doses required for the origin-anchored construction')
    if not all(left<right for left,right in zip(dose[:-1],dose[1:])):
        raise ValueError('This exact diagnostic requires strictly distinct reported doses')
    values=[exact(row['maximum_shortening_um']) for row in rows]
    scale=math.lcm(*(v.denominator for v in values))
    y=[int(value*scale) for value in values];n=len(y)
    cost,best,masks=enumerate_subsets(y)
    critical=sorted({0,*[max(0,left-right) for i,left in enumerate(y) for right in y[i+1:]]})
    result=[]
    for exceptions in range(n):
        size=n-exceptions;twice_delta=best[size];mask=masks[size]
        retained=[i for i in range(n) if mask&(1<<i)]
        fitted=[];last=0
        for i in retained:
            last=max(last,2*y[i]-twice_delta,0)
            if not abs(last-2*y[i])<=twice_delta:
                raise AssertionError('Matching construction failed')
            fitted.append({'record_id':rows[i]['record_id'],'dose':str(dose[i]),
                           'response_um':str(values[i]),'constructed_fitted_um':str(Fraction(last,2*scale))})
        feasible=maximum_compatible_count(y,twice_delta)
        smaller=[c for c in critical if c<twice_delta]
        previous=smaller[-1] if smaller else None
        previous_count=maximum_compatible_count(y,previous) if previous is not None else None
        assert feasible>=size
        assert previous_count is None or previous_count<size
        result.append({'hypothetical_exception_budget':exceptions,'required_retained':size,
                       'sharp_delta_um':str(Fraction(twice_delta,2*scale)),
                       'sharp_delta_float_um':float(Fraction(twice_delta,2*scale)),
                       'matching_retained_records':[rows[i]['record_id'] for i in retained],
                       'hypothetically_omitted_records':[rows[i]['record_id'] for i in range(n) if i not in retained],
                       'matching_construction':fitted,'independent_DP_maximum_retained':feasible,
                       'previous_critical_delta_um':str(Fraction(previous,2*scale)) if previous is not None else None,
                       'previous_DP_maximum_retained':previous_count})
    # Independent small cases, including a dip whose obstruction survives one removal.
    checks={
        'ordered_data_zero_budget_zero_error':enumerate_subsets([1,2,3])[1][3]==0,
        'one_reversal_half_drop':enumerate_subsets([5,1])[1][2]==4,
        'remove_one_reversal_zero_error':enumerate_subsets([5,1])[1][1]==0,
        'DP_two_disjoint_dips':maximum_compatible_count([5,1,5,1],0)==2,
        'all_exact_primal_and_independent_DP_checks':True,
        'all_seventeen_preserved_in_source':n==17,
        'full_data_floor_matches_prior_receipt':Fraction(best[n],2*scale)==Fraction('3.131943715'),
    }
    out=ROOT/'runs/calibration_uncertainty';out.mkdir(parents=True,exist_ok=True)
    array_path=out/'all_subset_costs.npy';np.save(array_path,cost)
    receipt={'status':'PASS_EXACT_CONDITIONAL_CONTAMINATION_DIAGNOSTIC','n':n,'enumerated_subsets':len(cost),
             'scale':scale,'checks':checks,'bounds':result,
             'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'enumeration_array_sha256':hashlib.sha256(array_path.read_bytes()).hexdigest(),
             'scope':'Sharp bounds in the relaxed nonnegative monotone curve class at fixed reported doses. Not an anatomical error estimate, actual row exclusion, or physical parameter identified set.'}
    (ROOT/'evidence/calibration_contamination.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    plot(result)
    print(json.dumps({'status':receipt['status'],'checks':checks,'first_four_bounds':[{k:r[k] for k in ['hypothetical_exception_budget','sharp_delta_float_um','hypothetically_omitted_records']} for r in result[:4]]}))
    if not all(checks.values()):raise SystemExit(1)


def plot(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('dark_background')
    fig,ax=plt.subplots(figsize=(9,4.8),constrained_layout=True)
    selected=rows[:7]
    ax.step([r['hypothetical_exception_budget'] for r in selected],[r['sharp_delta_float_um'] for r in selected],where='mid',color='#8ed9eb')
    ax.scatter([r['hypothetical_exception_budget'] for r in selected],[r['sharp_delta_float_um'] for r in selected],color='#f0bf7b')
    ax.set(xlabel='HYPOTHETICAL exception budget (all actual fits retain 17)',ylabel='Sharp required uniform discrepancy (µm)',xticks=range(7),ylim=(0,4),title='One exceptional embryo does not remove the common-curve discrepancy')
    fig.suptitle('Exact fixed-dose monotone compatibility; no measurement-noise or novelty claim',fontsize=11)
    fig.savefig(ROOT/'docs/assets/calibration_exceptions.png',dpi=170)
    plt.close(fig)


if __name__=='__main__':
    main()
