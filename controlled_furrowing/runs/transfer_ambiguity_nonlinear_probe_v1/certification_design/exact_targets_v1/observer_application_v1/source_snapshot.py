"""Transfer four independently certified spatial FE roots to frozen observers.

This is scalar exact-rational postprocessing. It imports no mechanics kernel.
The older cached-root application supplies the reviewed stored-curve geometry,
not the root error used here. The new error comes only from the exact-spatial
independent certificates. No calibration, observation or model is changed.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(__file__).resolve()
DESIGN = ROOT/'runs/transfer_ambiguity_nonlinear_probe_v1/certification_design'
OLD = DESIGN/'cached_pilots_v1/observer_application_v1'
EXACT = DESIGN/'exact_targets_v1'
OUT = EXACT/'observer_application_v1'
PINS = {
    OLD/'results.json':'21dc32c0b8bf073ffa9bf7a3c163cca107698e0f847304d53483caf071db1a40',
    OLD/'independent_readback/receipt.json':'939212f0690f4ccdf27ede99a5b7e7be6eff4a080adfd2dda78c0352c2e9c2e5',
    EXACT/'bridge_reduction_v1/independent_actual_preparation/all_four_readback.json':'859f507ad3762cfef44d53672bbe15a43592ce54d6e08ba8828bf50a925ec340',
}
REVIEW_STATUS = 'PASS_INDEPENDENT_EXACT_FOUR_OBSERVER_SOURCE'
LENGTH = Q(67,2)
OBS = ('maximum_shortening_um','region_averaged_shortening_um',
       'apical_max_shortening_um','basal_max_shortening_um')

def require(value, message):
    if not value:
        raise ValueError(message)

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def rel(path):
    return path.relative_to(ROOT).as_posix()

def read(path, expected, links):
    path=path.resolve()
    require(path.is_relative_to(ROOT), 'outside study')
    require(digest(path)==expected, 'changed authority '+str(path))
    require(links.setdefault(rel(path),expected)==expected, 'conflicting authority')
    return json.loads(path.read_bytes())

def save(path, value):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write('\n')
    return digest(path)

def interval(value):
    lo,hi=map(Q,value)
    require(lo<=hi,'reversed interval')
    return lo,hi

def inflate(value, error):
    lo,hi=interval(value)
    require(error>=0,'negative error')
    return [str(lo-error),str(hi+error)]

def classify(value, observed, delta):
    lo,hi=interval(value)
    require(delta>=0,'negative allowance')
    lo,hi=lo-observed,hi-observed
    lower=max(Q(0),lo,-hi)
    upper=max(abs(lo),abs(hi))
    return ('IN' if upper<=delta else 'OUT' if lower>delta else 'UNRESOLVED'),[str(lower),str(upper)]

def fixtures():
    require(classify(['-2','-1'],Q(0),Q(2))[0]=='IN','closed IN boundary')
    require(classify(['2','3'],Q(0),Q(2))[0]=='UNRESOLVED','strict OUT boundary')
    require(classify(['9/4','3'],Q(0),Q(2))[0]=='OUT','strict OUT')
    require(inflate(['1','2'],Q(1,4))==['3/4','9/4'],'outward inflation')
    # A graph with slope cap 3 has coefficient 5 and radial margin a-8 eps.
    eps0,eps1=Q(1,32),Q(1,16)
    require((5*eps0)*eps1/eps0==5*eps1 and Q(1)-8*eps1>0,
            'linear graph error and enlarged-error domain must be checked')
    return {'closed_IN_strict_OUT_and_unresolved_boundary':True,
            'outward_exact_inflation':True,'linear_graph_error_and_new_radial_domain':True}

def authorities():
    links={}
    docs={p:read(p,h,links) for p,h in PINS.items()}
    old=docs[OLD/'results.json']
    reviewed=docs[OLD/'independent_readback/receipt.json']
    new=docs[EXACT/'bridge_reduction_v1/independent_actual_preparation/all_four_readback.json']
    require(new['status']=='PASS_FOUR_INDEPENDENT_FIXED_INPUT_EXACT_SPATIAL_ROOTS_AND_MINIMA'
            and new['independent_exact_spatial_root_count']==4
            and new['full_five_gate_goal_achieved'] is False,'four independent spatial roots required')
    # Rehash complete independent spatial proof closure and old graph inputs.
    for name,h in new['current_file_sha256'].items():
        p=(ROOT/name).resolve()
        require(p.is_relative_to(ROOT) and digest(p)==h,'changed spatial proof '+name)
        require(links.setdefault(name,h)==h,'conflicting spatial input')
    for name,h in old['inputs_sha256'].items():
        p=(ROOT/name).resolve()
        require(p.is_relative_to(ROOT) and digest(p)==h,'changed graph premise '+name)
        require(links.setdefault(name,h)==h,'conflicting graph input')
    require(PINS[OLD/'results.json'] in json.dumps(reviewed), 'independent old application binding missing')
    require(len(old['four_pilots'])==4 and old['PDE_or_source_fit_claimed'] is False,
            'reviewed same four stored observers required')
    return old,new,links

def plan():
    old,new,links=authorities()
    require({r['pilot_id'] for r in old['four_pilots']}=={r['pilot_id'] for r in new['cases']},'same targets')
    return {'status':'FROZEN_FOUR_EXACT_SPATIAL_OBSERVER_TRANSFER','source_sha256':digest(SOURCE),
            'inputs_sha256':links,'length_um_exact':str(LENGTH),
            'allowance_um_exact':old['frozen_model_discrepancy_allowance_um_exact'],
            'scope':'Four fixed-input spatial finite-element roots; two small L states (32 alternative comparisons) and two large contrast values. No other calibration states, source fit, temporal path, PDE or novelty conclusion.'}

def prepare():
    p=plan(); checks=fixtures()
    OUT.mkdir(exist_ok=False)
    ph=save(OUT/'protocol.json',p)
    pf=save(OUT/'preflight.json',{'status':'PASS','source_sha256':digest(SOURCE),
               'protocol_sha256':ph,'checks':checks,'actual_observer_applications':0})
    print(json.dumps({'source_sha256':digest(SOURCE),'protocol_sha256':ph,'preflight_sha256':pf}))

def run(args):
    links={}
    p=read(OUT/'protocol.json',args.protocol_sha256,links)
    pf=read(OUT/'preflight.json',args.preflight_sha256,links)
    rv=read(args.review,args.review_sha256,links)
    binding={'source_sha256':digest(SOURCE),'protocol_sha256':args.protocol_sha256,'preflight_sha256':args.preflight_sha256}
    require(rv['status']==REVIEW_STATUS and rv['checks'] and all(v is True for v in rv['checks'].values())
            and all(rv[k]==v for k,v in binding.items()),'independent frozen source review required')
    require(pf['status']=='PASS' and pf['checks'] and all(v is True for v in pf['checks'].values())
            and type(pf['actual_observer_applications']) is int and pf['actual_observer_applications']==0
            and all(pf[k]==binding[k] for k in ('source_sha256','protocol_sha256')) and plan()==p,'cold phase changed')
    old,new,closure=authorities(); links.update(closure)
    delta=Q(p['allowance_um_exact']); rows=[]
    for prior in old['four_pilots']:
        case,=[r for r in new['cases'] if r['pilot_id']==prior['pilot_id']]
        ident=case['target_identity']; sm=case['selfmap']
        # Exact-spatial and cached operators have the identical evaluation center,
        # free map, field, drive and preconditioner; only integration differs.
        contract,=[read(ROOT/name,h,links) for name,h in old['inputs_sha256'].items()
                   if name.endswith('/model_contract.json') and prior['pilot_id']+'/' in name]
        require(ident['field_sha256']==contract['source_field_sha256']
                and ident['cache_packet_sha256']==contract['point_cache_sha256']
                and ident['R_binary_sha256']==contract['binary_preconditioner_sha256']
                and ident['model_contract_sha256']==digest(ROOT/next(name for name in old['inputs_sha256']
                    if name.endswith('/model_contract.json') and prior['pilot_id']+'/' in name))
                and Q(ident['dose_exact'])==Q(prior['exact_dose'])
                and Q(ident['gain_exact'])==Q(prior['exact_gain'])
                and Q(ident['drive_exact'])==Q(prior['exact_drive']),'same original point and physical input required')
        eta,z,r=Q(sm['eta_upper_exact']),Q(sm['z_upper_exact']),Q(case['radius_exact'])
        require(eta>=0 and 0<=z<1 and eta+z*r<r and sm['strict_selfmap'] is True,'strict spatial root required')
        error=eta/(1-z)
        require(error==Q(sm['posterior_error_upper_exact'])<r,'posterior mismatch')
        xi=Q(prior['serialized_X_plus_u_gap_exact']); eps=LENGTH*(error+xi)
        old_eps=Q(prior['physical_nodal_error_um_exact'])
        require(old_eps==LENGTH*(Q(prior['posterior_radius_exact'])+xi)>0,'old addition-gap premise')
        side_error={}; margins={}
        for side in ('apical','basal'):
            face=prior['graph_transfer'][side]
            slope=Q(face['base_graph_Lipschitz_upper_exact'])
            coeff=Q(5,4)*(1+slope)
            require(slope>=0 and Q(face['Eulerian_vertical_error_um_exact'])==coeff*old_eps,'reviewed graph coefficient')
            base=Q(face['radial_derivative_margin_exact'])+8*old_eps
            require(base-8*eps>0,'new error leaves monotone graph domain')
            margins[side]=str(base-8*eps); side_error[side]=coeff*eps
        ea,eb=side_error['apical'],side_error['basal']
        row={'pilot_id':prior['pilot_id'],'role':prior['role'],'endpoint_id':prior['endpoint_id'],
             'spatial_root_proof_sha256':case['independent_proof_sha256'],'posterior_error_exact':str(error),
             'xi_exact':str(xi),'physical_nodal_error_um_exact':str(eps),'new_radial_margins_exact':margins}
        if prior['role']=='large_final':
            eq=eps+ea
            row.update(stored_Q_um_exact=prior['stored_Q_um_exact'],Q_error_um_exact=str(eq),
                       exact_spatial_root_Q_um_exact=inflate(prior['stored_Q_um_exact'],eq))
        else:
            require(prior['record_id']=='optoGEF2_small_col_L','only L calibration covered')
            errors=dict(zip(OBS,(ea+eb,ea+eb,ea,eb),strict=True)); comparisons=[]
            for saved in prior['root_observer_comparisons']:
                iv=inflate(saved['saved_um_exact'],errors[saved['observable']])
                status,res=classify(iv,Q(saved['recorded_value_um_exact']),delta)
                comparisons.append({'variant':saved['variant'],'observable':saved['observable'],
                    'root_interval_um_exact':iv,'absolute_residual_um_exact':res,'status':status})
            require(len(comparisons)==16,'all L variants required'); row['comparisons']=comparisons
        rows.append(row)
    large={r['endpoint_id']:interval(r['exact_spatial_root_Q_um_exact']) for r in rows if r['role']=='large_final'}
    a,b=large['max_contrast_change'],large['max_negative_contrast_change']
    gap=[str(a[0]-b[1]),str(a[1]-b[0])]
    result={'status':'APPLIED_FOUR_INDEPENDENT_SPATIAL_ROOTS_AWAITING_READBACK','bindings':binding,
            'inputs_sha256':links,'cases':rows,'large_signed_gap_um_exact':gap,'strict_positive_gap':Q(gap[0])>0,
            'small_comparison_counts':dict(Counter(c['status'] for r in rows for c in r.get('comparisons',[]))),
            'full136_calibration_claimed':False,'source_fit_PDE_loading_or_novelty_claimed':False,'scope':p['scope']}
    rh=save(OUT/'results.json',result)
    print(json.dumps({'results_sha256':rh,'large_gap_approx':list(map(lambda v:float(Q(v)),gap)),
                      'counts':result['small_comparison_counts']}))

if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__); sub=cli.add_subparsers(dest='mode',required=True)
    sub.add_parser('prepare'); actual=sub.add_parser('run')
    for k in ('protocol-sha256','preflight-sha256','review-sha256'): actual.add_argument('--'+k,required=True)
    actual.add_argument('--review',type=Path,required=True)
    args=cli.parse_args()
    prepare() if args.mode=='prepare' else run(args)
