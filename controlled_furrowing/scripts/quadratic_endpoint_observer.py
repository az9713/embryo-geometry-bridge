"""Conditional four-endpoint observers on current-radius quadratic graphs.

Exact rationals represent serialized node values. Eulerian line/disk means
use polynomial antiderivatives, with certified quadratic-inverse brackets at
ROI cuts. Source aggregation/reference/registration ambiguities are explicit.
No solver or physical-model source is modified or optimized by this module.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path
import argparse,hashlib,heapq,json
from annulus_observer import Interval,rational,serial,intersection,_quadratic_controls,_root_bracket,_poly

ROOT=Path(__file__).resolve().parents[1]

def interval(lo,hi=None):
    lo=rational(lo);hi=lo if hi is None else rational(hi)
    if lo>hi:raise ValueError('ordered interval required')
    return Interval(lo,hi)

def add(a,b):return interval(a.lo+b.lo,a.hi+b.hi)
def scale(a,c):
    c=rational(c);return interval(min(c*a.lo,c*a.hi),max(c*a.lo,c*a.hi))
def subtract(a,b):return add(a,scale(b,-1))
def multiply(a,b):
    values=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi];return interval(min(values),max(values))
def polynomial_interval(coeff,t):
    result=interval(0)
    for c in reversed(coeff):result=add(multiply(result,t),interval(c))
    return result
def polynomial_product(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def roi_pair(roi):
    lo,hi=map(rational,roi)
    if not 0<=lo<hi:raise ValueError('positive-width nonnegative current-radius ROI required')
    return lo,hi

class Curve:
    """A continuous radially monotone piecewise quadratic physical graph."""
    def __init__(self,controls,root_steps=80):
        if not isinstance(root_steps,int) or root_steps<1:raise ValueError('positive integer root steps required')
        self.curves=_quadratic_controls(controls);self.root_steps=root_steps
        self.domain=(self.curves[0][0][0][0],self.curves[-1][0][2][0])

    def guard_roi(self,roi):
        lo,hi=map(rational,roi)
        if lo>hi or lo<self.domain[0] or hi>self.domain[1]:raise ValueError('current-radius ROI not covered by graph')
        return lo,hi

    def cuts(self,roi):
        lo,hi=self.guard_roi(roi)
        for q,c0,c1,c2 in self.curves:
            cut=intersection(interval(lo,hi),interval(q[0][0],q[2][0]))
            if cut is None:continue
            left=_root_bracket(c0[0],c1[0],c2[0],cut.lo,self.root_steps)
            right=_root_bracket(c0[0],c1[0],c2[0],cut.hi,self.root_steps)
            yield q,c0,c1,c2,left,right

    def range(self,lo,hi):
        values=[]
        for q,c0,c1,c2,left,right in self.cuts((lo,hi)):
            a,b=left[0],right[1]
            values.extend((_poly(c0[1],c1[1],c2[1],a),_poly(c0[1],c1[1],c2[1],b)))
            if c2[1]:
                critical=-c1[1]/(2*c2[1])
                if a<=critical<=b:values.append(_poly(c0[1],c1[1],c2[1],critical))
        return interval(min(values),max(values))

    def knots(self):
        values=[]
        for q,c0,c1,c2 in self.curves:
            values.extend((q[0][0],q[2][0]))
            if c2[1]:
                critical=-c1[1]/(2*c2[1])
                if 0<critical<1:values.append(_poly(c0[0],c1[0],c2[0],critical))
        return values

    def integral(self,roi,weighting):
        if weighting not in ('current_radial_line','current_disk'):raise ValueError('explicit current_radial_line or current_disk weighting required')
        total=interval(0)
        for q,c0,c1,c2,left,right in self.cuts(roi):
            R=[c0[0],c1[0],c2[0]];Z=[c0[1],c1[1],c2[1]];dR=[c1[0],2*c2[0]]
            p=polynomial_product(Z,dR)
            if weighting=='current_disk':p=polynomial_product(p,R)
            primitive=[Q(0)]+[v/Q(i+1) for i,v in enumerate(p)]
            lower=polynomial_interval(primitive,interval(*left));upper=polynomial_interval(primitive,interval(*right))
            total=add(total,subtract(upper,lower))
        return total

def mass(roi,weighting):
    lo,hi=roi_pair(roi)
    if weighting=='current_radial_line':return hi-lo
    if weighting=='current_disk':return (hi*hi-lo*lo)/2
    raise ValueError('unknown spatial weighting')

class Field:
    def __init__(self,constant,terms):self.constant=rational(constant);self.terms=[(rational(c),g) for c,g in terms]
    def range(self,lo,hi):
        value=interval(self.constant)
        for c,g in self.terms:value=add(value,scale(g.range(lo,hi),c))
        return value
    def knots(self):return [r for c,g in self.terms for r in g.knots()]
    def integral(self,roi,weighting):
        value=interval(self.constant*mass(roi,weighting))
        for c,g in self.terms:value=add(value,scale(g.integral(roi,weighting),c))
        return value

def transformed_graph(graph,constant,coefficient):
    return Curve([[(r,constant+coefficient*z) for r,z in q] for q,c0,c1,c2 in graph.curves],graph.root_steps)

def shared_radius_difference(apical,basal,reference_height):
    if len(apical.curves)!=len(basal.curves):return None
    controls=[]
    for (a,*unused),(b,*unused2) in zip(apical.curves,basal.curves):
        if [p[0] for p in a]!=[p[0] for p in b]:return None
        controls.append([(a[i][0],reference_height-a[i][1]+b[i][1]) for i in range(3)])
    return Curve(controls,apical.root_steps)

def graph_maximum(graph,roi):
    lo,hi=roi_pair(roi);bound=graph.range(lo,hi)
    points=sorted(set([lo,hi]+[p for p in graph.knots() if lo<=p<=hi]))
    lower=max(graph.range(p,p).lo for p in points)
    return interval(lower,bound.hi)

def field_supremum(field,roi,tolerance=Q('0.001'),maximum_cells=8192):
    lo,hi=roi_pair(roi);tolerance=rational(tolerance)
    if tolerance<=0 or maximum_cells<1:raise ValueError('positive reporting tolerance and cell budget required')
    knots=sorted(set([lo,hi]+[p for p in field.knots() if lo<p<hi]))
    lower=None;heap=[];counter=0
    def cell(left,right):
        nonlocal lower,counter
        upper=field.range(left,right).hi
        for p in (left,(left+right)/2,right):
            value=field.range(p,p).lo;lower=value if lower is None else max(lower,value)
        counter+=1;return (-upper,counter,left,right)
    for left,right in zip(knots[:-1],knots[1:]):heapq.heappush(heap,cell(left,right))
    while -heap[0][0]-lower>tolerance and len(heap)<maximum_cells:
        _,_,a,b=heapq.heappop(heap);mid=(a+b)/2
        heapq.heappush(heap,cell(a,mid));heapq.heappush(heap,cell(mid,b))
    result=interval(lower,-heap[0][0])
    return {'interval':result,'requested_width_um':tolerance,'requested_width_achieved':result.hi-result.lo<=tolerance,
            'active_cells':len(heap),'maximum_cells':maximum_cells,'method':'Exact rational adaptive signed supremum, not a sampled-point estimate.'}

def positive_thickness(apical,basal,roi,tolerance,maximum_cells):
    thickness=Field(0,[(1,apical),(-1,basal)]);b=thickness.range(*roi)
    if b.lo>0:return b.lo
    negative=Field(0,[(-1,apical),(1,basal)]);test=field_supremum(negative,roi,tolerance,maximum_cells)['interval']
    if test.lo>=0:raise ValueError('nonpositive current Eulerian tissue thickness at a certified point')
    if test.hi>=0:raise ValueError('positive current Eulerian thickness not certified; do not classify this as a proven inversion')
    return -test.hi

def observe_graphs(apical_controls,basal_controls,*,roi_um,reference_apical_um,reference_basal_um=0,
                   weighting='current_disk',component_aggregation='spatial_max',tolerance_um='0.001',maximum_cells=8192,root_steps=80):
    """One static state; references/ROI/weighting are explicit assumptions.

    Components remain signed. No positive-part clipping or temporal peak is
    inferred. Both source-plausible component aggregation variants are retained.
    """
    roi=roi_pair(roi_um);a0=rational(reference_apical_um);b0=rational(reference_basal_um)
    if a0<=b0:raise ValueError('positive flat reference thickness required')
    if weighting not in ('current_radial_line','current_disk'):raise ValueError('declared current-radius spatial weighting required')
    if component_aggregation not in ('spatial_max','regional_mean'):raise ValueError('declared component aggregation required')
    a=Curve(apical_controls,root_steps);b=Curve(basal_controls,root_steps);a.guard_roi(roi);b.guard_roi(roi)
    thickness_lower=positive_thickness(a,b,roi,rational(tolerance_um),maximum_cells)
    apical=transformed_graph(a,a0,-1);basal=transformed_graph(b,-b0,1)
    A=graph_maximum(apical,roi);B=graph_maximum(basal,roi)
    shortening=Field(a0-b0,[(-1,a),(1,b)])
    shared=shared_radius_difference(a,b,a0-b0)
    if shared is None:maximum=field_supremum(shortening,roi,rational(tolerance_um),maximum_cells)
    else:
        value=graph_maximum(shared,roi)
        maximum={'interval':value,'requested_width_um':rational(tolerance_um),'requested_width_achieved':value.hi-value.lo<=rational(tolerance_um),
                 'active_cells':len(shared.curves),'method':'Exact shared-radial-parameter quadratic extrema; no independent-radius subtraction.'}
    alternatives={}
    for weight in ('current_radial_line','current_disk'):
        denominator=mass(roi,weight)
        amean=add(interval(a0),scale(a.integral(roi,weight),-1/denominator))
        bmean=add(interval(-b0),scale(b.integral(roi,weight),1/denominator))
        alternatives[weight]={'apical_regional_mean_downward':amean,'basal_regional_mean_upward':bmean,'height_regional_mean_shortening':add(amean,bmean),
            'method':'Polynomial change-of-variable antiderivative with exact rational quadratic-inverse endpoint brackets.'}
    means=alternatives[weighting]
    components=(A,B) if component_aggregation=='spatial_max' else (means['apical_regional_mean_downward'],means['basal_regional_mean_upward'])
    return {'status':'CERTIFIED_STATIC_INTERVALS','four_endpoints':{'apical_downward':components[0],'basal_upward':components[1],
                'maximum_Eulerian_height_shortening':maximum['interval'],'regional_mean_Eulerian_height_shortening':means['height_regional_mean_shortening']},
            'spatial_max_components':{'apical_downward':A,'basal_upward':B},'mean_alternatives':alternatives,'maximum_shortening_certificate':maximum,
            'positive_thickness_lower_um':thickness_lower,
            'contract':{'current_radius_ROI_um':roi,'reference_apical_um':a0,'reference_basal_um':b0,'weighting':weighting,'component_aggregation':component_aggregation,
                'root_bracketing_steps':root_steps,'state_scope':'One supplied static state only; no time/path extremum implied.',
                'source_alignment':'Conditional references and registered current-radius ROI; source does not uniquely fix component reference split, component aggregation, ROI transport or line/disk weighting.',
                'signed':'No positive-part clipping. Negative local/component/mean shortening is permitted.',
                'not_proved':'Continuum mechanics/biological validity, experimental temporal coverage, statistical error or canonical source observation operator.'}}

def _state_controls(state,length_unit_um):
    p=state['parameters']
    if p.get('reference_curvature',0)!=0:raise ValueError('flat reference required by this observer')
    L=rational(length_unit_um)
    if L<=0:raise ValueError('positive declared physical length unit required')
    ref=[[rational(v) for v in q] for q in state['reference_nodes']];cur=[[rational(v) for v in q] for q in state['current_nodes']]
    tri=state['triangles']
    if any(len(q)!=2 for q in ref+cur):raise ValueError('two-dimensional meridional node coordinates required')
    if len(ref)!=len(cur) or any(len(t)!=6 for t in tri):raise ValueError('matching P2 reference/current nodes and six-node triangles required')
    n=max(i for t in tri for i in t[:3])+1
    if set(i for t in tri for i in t[:3])!=set(range(n)):raise ValueError('dense vertex prefix required')
    H=rational(p['height']);top=max(q[1] for q in ref[:n]);bottom=min(q[1] for q in ref[:n])
    if bottom!=0 or top!=H:raise ValueError('declared flat reference height disagrees with saved nodes')
    edges={}
    for a,b,c,d,e,f in tri:
        for i,j,m in ((a,b,d),(b,c,e),(c,a,f)):
            if min(i,j,m)<0 or max(i,j,m)>=len(ref):raise ValueError('invalid P2 connectivity')
            key=tuple(sorted((i,j)))
            if key in edges and edges[key]!=m:raise ValueError('nonconforming quadratic surface trace')
            edges[key]=m
    output={}
    for name,z in (('top',H),('bottom',Q(0))):
        boundary=[]
        for (i,j),m in edges.items():
            if ref[i][1]==ref[j][1]==z:
                if ref[i][0]>ref[j][0]:i,j=j,i
                boundary.append((ref[i][0],[cur[i],cur[m],cur[j]]))
        controls=[v for r,v in sorted(boundary)]
        saved=state.get('boundary_'+name+'_quadratic_controls')
        if not controls:raise ValueError('missing quadratic boundary')
        if saved is not None and controls!=[[[rational(v) for v in q] for q in e] for e in saved]:raise ValueError('saved boundary controls disagree with current-node/connectivity field')
        output[name]=[[(L*r,L*z) for r,z in e] for e in controls]
    return output,H*L

def observe_state(state,*,roi_um=None,length_unit_um='33.5',reference_apical_um=None,reference_basal_um=0,
                  weighting='current_disk',component_aggregation='spatial_max',**kwargs):
    if not state.get('converged',False):raise ValueError('saved state must be a verified equilibrium in its declared rule')
    controls,H=_state_controls(state,length_unit_um)
    if roi_um is None:roi_um=(0,rational(state['parameters']['patch_radius'])*rational(length_unit_um))
    result=observe_graphs(controls['top'],controls['bottom'],roi_um=roi_um,
        reference_apical_um=H if reference_apical_um is None else reference_apical_um,reference_basal_um=reference_basal_um,
        weighting=weighting,component_aggregation=component_aggregation,**kwargs)
    result['contract']['reference_choice']='Flat model initial reference unless explicitly overridden; source counterfactual growth correction is not inferred from target observations.'
    result['contract']['default_ROI_choice']='Fixed registered current-radius ROI equal to initial patch radius; a conditional assumption, not a recovered PIV transport rule.'
    result['activation_strength']=state.get('activation_strength')
    return result

def finite_path_mean_of_local_peaks(fields,roi,weighting,tolerance=Q('0.001'),maximum_cells=8192):
    if not fields or rational(tolerance)<=0 or not isinstance(maximum_cells,int) or maximum_cells<1:raise ValueError('nonempty fields and positive numerical budgets required')
    denominator=mass(roi,weighting);lo,hi=roi_pair(roi);heap=[];counter=0;low=Q(0);high=Q(0)
    knots=sorted(set([lo,hi]+[r for f in fields for r in f.knots() if lo<r<hi]))
    def cell(a,b):
        nonlocal counter
        ranges=[f.range(a,b) for f in fields];winner=None
        for i,bound in enumerate(ranges):
            if all(i==j or bound.lo>=other.hi for j,other in enumerate(ranges)):winner=i;break
        if winner is not None:value=fields[winner].integral((a,b),weighting)
        else:value=scale(interval(max(v.lo for v in ranges),max(v.hi for v in ranges)),mass((a,b),weighting))
        counter+=1;return (-(value.hi-value.lo),counter,a,b,value)
    for a,b in zip(knots[:-1],knots[1:]):
        c=cell(a,b);heapq.heappush(heap,c);low+=c[-1].lo;high+=c[-1].hi
    tolerance=rational(tolerance)
    while (high-low)/denominator>tolerance and len(heap)<maximum_cells:
        old=heapq.heappop(heap);_,_,a,b,value=old;low-=value.lo;high-=value.hi;mid=(a+b)/2
        for left,right in ((a,mid),(mid,b)):
            c=cell(left,right);heapq.heappush(heap,c);low+=c[-1].lo;high+=c[-1].hi
    value=interval(low/denominator,high/denominator)
    return {'interval':value,'requested_width_um':tolerance,'requested_width_achieved':value.hi-value.lo<=tolerance,
            'active_cells':len(heap),'scope':'Mean of pointwise maxima over ONLY the supplied finite state collection, not an unsampled time/path supremum.'}

def observe_graph_path(states,*,roi_um,reference_apical_um,reference_basal_um=0,weighting='current_disk',component_aggregation='spatial_max',**kwargs):
    """Finite explicit graph-state collection sharing registered ROI/reference."""
    if not states:raise ValueError('nonempty explicit state collection required')
    results=[observe_graphs(a,b,roi_um=roi_um,reference_apical_um=reference_apical_um,reference_basal_um=reference_basal_um,weighting=weighting,component_aggregation=component_aggregation,**kwargs) for a,b in states]
    keys=tuple(results[0]['four_endpoints'])
    peaks={k:interval(max(v['four_endpoints'][k].lo for v in results),max(v['four_endpoints'][k].hi for v in results)) for k in keys}
    H=rational(reference_apical_um)-rational(reference_basal_um)
    fields=[Field(H,[(-1,Curve(a,kwargs.get('root_steps',80))),(1,Curve(b,kwargs.get('root_steps',80)))]) for a,b in states]
    local=finite_path_mean_of_local_peaks(fields,roi_um,weighting,rational(kwargs.get('tolerance_um','0.001')),kwargs.get('maximum_cells',8192))
    return {'sampled_state_peak_endpoints':peaks,'mean_of_pointwise_sampled_state_peaks':local,'static_states':results,
            'scope':'Only the provided finite states. No invented temporal coverage or continuous quasistatic-path maximum; peak-of-mean and mean-of-peaks remain distinct.'}

def observe_state_path(states,*,roi_um=None,length_unit_um='33.5',reference_apical_um=None,reference_basal_um=0,
                       weighting='current_disk',component_aggregation='spatial_max',**kwargs):
    """A finite supplied same-model collection; not a continuous-time path."""
    if not states:raise ValueError('nonempty explicit saved-state collection required')
    p=states[0]['parameters'];graphs=[];height=None
    for state in states:
        if not state.get('converged') or state['parameters']!=p:raise ValueError('finite path requires verified states with identical model/mesh parameters')
        controls,H=_state_controls(state,length_unit_um)
        if height is not None and H!=height:raise ValueError('fixed registered reference height required')
        height=H;graphs.append((controls['top'],controls['bottom']))
    if roi_um is None:roi_um=(0,rational(p['patch_radius'])*rational(length_unit_um))
    result=observe_graph_path(graphs,roi_um=roi_um,reference_apical_um=height if reference_apical_um is None else reference_apical_um,
        reference_basal_um=reference_basal_um,weighting=weighting,component_aggregation=component_aggregation,**kwargs)
    result['supplied_activation_strengths']=[state.get('activation_strength') for state in states]
    result['coverage']='Finite supplied states only; no unsampled strength or experimental-time coverage is certified.'
    return result

def _tests():
    checks={}
    def graph(fn,Rfn=lambda t:t):return [[(Rfn(t),fn(t)) for t in (Q(0),Q(1,2),Q(1))]]
    def contains(value,x):return value.lo<=x<=value.hi
    def run(a,b,H=10,roi=(0,1),**kwargs):return observe_graphs(a,b,roi_um=roi,reference_apical_um=H,tolerance_um='0.0001',**kwargs)
    constant=run(graph(lambda t:8),graph(lambda t:1))
    checks['constant_four_exact']=all(v.lo==v.hi==q for v,q in zip(constant['four_endpoints'].values(),(Q(2),Q(1),Q(3),Q(3))))
    tilted=run(graph(lambda t:10-t),graph(lambda t:t/2))
    checks['tilted_component_maxima']=contains(tilted['four_endpoints']['apical_downward'],Q(1)) and contains(tilted['four_endpoints']['basal_upward'],Q(1,2))
    checks['tilted_line_mean_exact']=tilted['mean_alternatives']['current_radial_line']['height_regional_mean_shortening'].lo==tilted['mean_alternatives']['current_radial_line']['height_regional_mean_shortening'].hi==Q(3,4)
    checks['tilted_disk_mean_exact']=tilted['mean_alternatives']['current_disk']['height_regional_mean_shortening'].lo==tilted['mean_alternatives']['current_disk']['height_regional_mean_shortening'].hi==1
    quad=run(graph(lambda t:10-t*t),graph(lambda t:t*t/5))
    checks['quadratic_line_mean_exact']=quad['mean_alternatives']['current_radial_line']['height_regional_mean_shortening'].lo==quad['mean_alternatives']['current_radial_line']['height_regional_mean_shortening'].hi==Q(2,5)
    checks['quadratic_disk_mean_exact']=quad['mean_alternatives']['current_disk']['height_regional_mean_shortening'].lo==quad['mean_alternatives']['current_disk']['height_regional_mean_shortening'].hi==Q(3,5)
    split=run(graph(lambda t:10-t),graph(lambda t:1-t))
    checks['sum_of_component_maxima_not_height_max']=split['four_endpoints']['apical_downward'].lo+split['four_endpoints']['basal_upward'].lo==2 and split['four_endpoints']['maximum_Eulerian_height_shortening'].hi==1
    moving=run(graph(lambda t:10-t,lambda t:Q(3,2)*t-t*t/2),graph(lambda t:t),roi=(0,Q(5,8)))
    checks['Eulerian_height_excludes_material_column_proxy']=contains(moving['four_endpoints']['maximum_Eulerian_height_shortening'],Q(9,8)) and moving['four_endpoints']['maximum_Eulerian_height_shortening'].hi<Q(5,4)
    negative=run(graph(lambda t:11),graph(lambda t:-1))
    checks['signed_values_not_clipped']=all(v.lo==v.hi==q for v,q in zip(negative['four_endpoints'].values(),(Q(-1),Q(-1),Q(-2),Q(-2))))
    regional=run(graph(lambda t:10-t),graph(lambda t:t/2),component_aggregation='regional_mean')
    checks['component_mean_variant_explicit']=regional['four_endpoints']['apical_downward'].lo==Q(2,3) and regional['four_endpoints']['basal_upward'].lo==Q(1,3)
    path=[(graph(lambda t:3-t),graph(lambda t:0)),(graph(lambda t:Q(14,5)),graph(lambda t:0))]
    path_result=observe_graph_path(path,roi_um=(0,1),reference_apical_um=3,weighting='current_radial_line',tolerance_um='0.00001')
    checks['finite_path_peak_of_mean']=contains(path_result['sampled_state_peak_endpoints']['regional_mean_Eulerian_height_shortening'],Q(1,2))
    checks['finite_path_mean_of_peaks_differs']=contains(path_result['mean_of_pointwise_sampled_state_peaks']['interval'],Q(13,25)) and path_result['mean_of_pointwise_sampled_state_peaks']['interval'].lo>Q(1,2)
    guards=[lambda:run([[(0,10),(Q(4,5),10),(1,10)]],graph(lambda t:0)),lambda:run(graph(lambda t:10),graph(lambda t:0),roi=(0,2)),lambda:run(graph(lambda t:0),graph(lambda t:1)),lambda:run(graph(lambda t:10),graph(lambda t:0),component_aggregation='unspecified')]
    for i,call in enumerate(guards):
        try:call();checks[f'guard_{i}']=False
        except ValueError:checks[f'guard_{i}']=True
    identity={'converged':True,'activation_strength':0.,'parameters':{'height':10.,'patch_radius':.5,'reference_curvature':0.},
              'reference_nodes':[[0,0],[1,0],[1,10],[0,10],[.5,0],[1,5],[.5,5],[.5,10],[0,5]],
              'triangles':[[0,1,2,4,5,6],[0,2,3,6,7,8]]}
    identity['current_nodes']=[p[:] for p in identity['reference_nodes']]
    out=observe_state(identity,length_unit_um=1)
    checks['saved_flat_P2_state_parser_identity']=all(v.lo==v.hi==0 for v in out['four_endpoints'].values())
    import copy
    bad=copy.deepcopy(identity);bad['parameters']['height']=11.
    mismatch=copy.deepcopy(identity);mismatch['boundary_top_quadratic_controls']=[[[0,10],[.5,9],[1,10]]]
    bad_dimension=copy.deepcopy(identity);bad_dimension['current_nodes'][0]=[0,0,0]
    moreguards=[lambda:observe_state(bad,length_unit_um=1),lambda:observe_state(mismatch,length_unit_um=1),lambda:run(graph(lambda t:10),graph(lambda t:0),weighting='unspecified'),lambda:observe_state_path([identity,bad],length_unit_um=1),lambda:observe_state(bad_dimension,length_unit_um=1),lambda:observe_state_path([],length_unit_um=1)]
    for i,call in enumerate(moreguards):
        try:call();checks[f'input_guard_{i}']=False
        except ValueError:checks[f'input_guard_{i}']=True
    receipt={'status':'PASS_EXACT_OBSERVER_SYNTHETICS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'examples':{'constant':constant,'tilted':tilted,'quadratic':quad,'separate_component_peaks':split,'radial_motion':moving,'negative_signed':negative,'finite_path_ordering':path_result},
             'scope':'Analytical graph tests, not biological validity or continuum solver accuracy.'}
    (ROOT/'evidence/quadratic_endpoint_observer_synthetic.json').write_text(json.dumps(serial(receipt),indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'checks':len(checks),'failures':[k for k,v in checks.items() if not v]}),flush=True)
    if not all(checks.values()):raise RuntimeError('observer synthetic checks failed')

def _diagnostics():
    paths=[ROOT/'runs/activated_quadrature'/f'lambda_{lam}_strength_{strength}'/'reequilibrated.json' for lam in (500.,5000.) for strength in (1.,5.,10.)]
    paths += [ROOT/'runs/joint_refinement/final.json',ROOT/'runs/sparse_quadratic/24x12_lambda5000/gain_1.0/final.json']
    rows=[]
    for path in paths:
        state=json.loads(path.read_text());out=observe_state(state)
        row={'file':path.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'mesh':[state['parameters']['nr'],state['parameters']['nz']],
             'lambda':state['parameters']['bulk'],'strength':state['activation_strength'],'observer':out,
             'legacy_material_column_shortening_um':state['maximum_activated_column_shortening']*33.5,
             'legacy_scope':'Diagnostic comparison only; this material-column proxy is NOT substituted for current-radius height shortening.'}
        rows.append(row)
        (ROOT/'evidence/quadratic_endpoint_observer_diagnostics.json').write_text(json.dumps(serial({'status':'COMPLETE' if len(rows)==len(paths) else 'RUNNING','rows':rows,
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Existing verified static fields only; no solver calls, source-data fits or invented temporal extrema.'}),indent=2)+'\n')
        print(json.dumps(serial({'file':row['file'],'four_endpoints':out['four_endpoints'],'max_resolution_met':out['maximum_shortening_certificate']['requested_width_achieved']})),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--test',action='store_true');parser.add_argument('--diagnostics',action='store_true');args=parser.parse_args()
    if args.test:_tests()
    if args.diagnostics:_diagnostics()
