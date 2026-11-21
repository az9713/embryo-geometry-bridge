"""Separate v2 tension revision generated from retained v1; see generation receipt."""
"""Calibration-only nonlinear spatial experiment, with checkpointed evidence.

This does not silently pass the original goal's unresolved identified-set,
specimen geometry, transient measurement, or continuum-accuracy requirements.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import torch
from axisymmetric_prestress import Parameters,build_model
from spatial_observer import observe

ROOT=Path(__file__).resolve().parents[1]
torch.set_num_threads(1)


def load(relative):
    return json.loads((ROOT/relative).read_text(encoding='utf-8'))


def save(relative,value):
    path=ROOT/relative;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def sha(relative):
    return hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()


def csv_rows(relative):
    with (ROOT/relative).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))


def manifest():
    files=['SPATIAL-TENSION-PROTOCOL.json','scripts/spatial_tension_transfer.py','scripts/axisymmetric_mechanics.py','scripts/axisymmetric_prestress.py','scripts/spatial_observer.py','data/derived/calibration_small_18um.csv','data/derived/spatial_inputs_only.json']
    return {name:sha(name) for name in files}


def prepare():
    rows=csv_rows('data/derived/evaluation_larger_bending.csv')
    data=[{'record_id':r['record_id'],'reported_diameter_um':float(r['reported_diameter_um']),'activation_magnitude':float(r['activation_magnitude'])} for r in rows]
    save('data/derived/spatial_inputs_only.json',data)
    print(json.dumps({'inputs_only_records':len(data),'held_out_response_fields':False}))


def build(kb,diameter,protocol):
    h=protocol['units']['height_um'];geometry=protocol['geometry'];material=protocol['material_scenarios'];mesh=protocol['discretization']
    return build_model(Parameters(nr=mesh['nr'],nz=mesh['nz'],height=1.,outer_radius=geometry['outer_radius_over_height'],patch_radius=diameter/(2*h),activation_z_min=(h-25)/h,activation_z_max=(h-5)/h,mu=material['mu'],bulk=material['logJ_penalty_coefficient_over_mu'],basal_area_modulus=material['basal_area_modulus_over_mu_height'],basal_surface_tension=kb,apical_surface_tension=0.,apical_area_modulus=material['apical_area_modulus_over_mu_height'],reference_curvature=0.,quadrature=mesh['quadrature']))


def branch_id(kb):
    return f'tension_{kb:g}'


def checked_relax(model,strength,initial):
    result=model.relax(float(strength),initial_free_vector=initial,maxiter=100)
    if not result['converged']:
        # One fixed, target-independent retry from the same state with more
        # L-BFGS iterations. Preserve the first attempt, never hide a failure.
        first=result
        result=model.relax(float(strength),initial_free_vector=initial,maxiter=1200)
        result['retained_first_attempt']={k:first[k] for k in ['converged','scaled_free_force_linf','minimum_J','continuation']}
    return result


def table(kb):
    protocol=load('SPATIAL-TENSION-PROTOCOL.json');hashes=manifest();branch=branch_id(kb)
    if kb not in protocol['material_scenarios']['basal_surface_tension_over_mu_height']:
        raise ValueError('Undeclared material branch')
    model=build(kb,protocol['geometry']['small_diameter_um'],protocol)
    nodes=protocol['calibration']['initial_response_strength_nodes']
    initial=None;records=[];started=time.monotonic()
    for index,strength in enumerate(nodes):
        path=f'runs/spatial_tension/{branch}/table_{index:03d}.json'
        if (ROOT/path).exists():
            record=load(path)
            if record['hashes']!=hashes:raise RuntimeError('Model/protocol changed: retain old run and use a versioned experiment, not stale checkpoints')
            result=record['result']
        else:
            result=checked_relax(model,strength,initial)
            record={'hashes':hashes,'strength':strength,'result':result}
            save(path,record)
        if not result['converged']:
            save(f'evidence/spatial_tension_{branch}_failure.json',{'status':'FAIL','stage':'calibration_table','index':index,'strength':strength,'hashes':hashes,'result':result})
            raise RuntimeError('Nonstationary or nonphysical table point; branch not silently dropped')
        initial=np.array(result['free_displacement'])
        records.append({'strength':strength,'shortening_um':protocol['units']['height_um']*result['maximum_activated_column_shortening'],'path':path})
        print(json.dumps({'branch':branch,'table_index':index,'strength':strength,'shortening_um':records[-1]['shortening_um'],'force':result['scaled_free_force_linf'],'seconds':time.monotonic()-started}),flush=True)
    if hashes!=manifest():raise RuntimeError('Inputs changed during table run')
    save(f'runs/spatial_tension/{branch}/table.json',{'hashes':hashes,'records':records,'scope':'Conditional discrete-response table, not a continuum enclosure'})


def calibrate(kb):
    protocol=load('SPATIAL-TENSION-PROTOCOL.json');branch=branch_id(kb);tab=load(f'runs/spatial_tension/{branch}/table.json');hashes=manifest()
    if tab['hashes']!=hashes:raise RuntimeError('Table provenance differs')
    rows=csv_rows(protocol['calibration']['input'])
    dose=np.array([float(r['activation_magnitude']) for r in rows]);response=np.array([float(r['maximum_shortening_um']) for r in rows])
    strengths=np.array([r['strength'] for r in tab['records']]);values=np.array([r['shortening_um'] for r in tab['records']])
    grid=protocol['calibration']['gain_grid']
    gains=np.linspace(grid['minimum'],grid['maximum'],round((grid['maximum']-grid['minimum'])/grid['step'])+1)
    if np.max(gains[:,None]*dose)>max(strengths):raise RuntimeError('Surrogate would extrapolate')
    predicted=np.interp(gains[:,None]*dose,strengths,values)
    rss=np.sum((predicted-response)**2,axis=1);selected=int(np.argmin(rss));gain=float(gains[selected])
    model=build(kb,protocol['geometry']['small_diameter_um'],protocol)
    initial=None;actual=np.zeros(len(rows));diagnostics=[]
    for index in np.argsort(dose):
        strength=gain*dose[index]
        result=checked_relax(model,strength,initial)
        path=f'runs/spatial_tension/{branch}/calibration_direct_{index:02d}.json'
        save(path,{'hashes':hashes,'record_id':rows[index]['record_id'],'strength':strength,'result':result})
        if not result['converged']:raise RuntimeError('Direct calibration solve failed; outcome retained')
        actual[index]=protocol['units']['height_um']*result['maximum_activated_column_shortening']
        initial=np.array(result['free_displacement']);diagnostics.append({'index':int(index),'path':path,'force':result['scaled_free_force_linf']})
    error=float(np.max(abs(actual-predicted[selected])))
    result={'status':'CALIBRATION_CHECKED' if error<=protocol['calibration']['direct_response_surrogate_tolerance_um'] else 'SURROGATE_TOLERANCE_FAIL',
            'hashes':hashes,'branch':branch,'basal_surface_tension':kb,'gain':gain,'gain_at_boundary':selected in [0,len(gains)-1],
            'surrogate_calibration_RSS_um2':float(rss[selected]),'direct_calibration_RSS_um2':float(np.sum((actual-response)**2)),
            'maximum_surrogate_error_um':error,'direct_predictions_um':actual.tolist(),'record_ids':[r['record_id'] for r in rows],
            'gain_grid':gains.tolist(),'surrogate_gain_RSS_um2':rss.tolist(),'checks':diagnostics,
            'scope':'One approximate finite-grid gain fit per declared material scenario; no material identified-set or confidence claim'}
    save(f'runs/spatial_tension/{branch}/calibration.json',result)
    print(json.dumps({k:result[k] for k in ['status','branch','gain','gain_at_boundary','direct_calibration_RSS_um2','maximum_surrogate_error_um']}))
    if result['status']!='CALIBRATION_CHECKED':raise RuntimeError('Do not proceed to target predictions until calibration approximation is resolved')


def predict(kb):
    protocol=load('SPATIAL-TENSION-PROTOCOL.json');branch=branch_id(kb);calibration=load(f'runs/spatial_tension/{branch}/calibration.json');hashes=manifest()
    if calibration['hashes']!=hashes or calibration['status']!='CALIBRATION_CHECKED':raise RuntimeError('Unverified calibration')
    inputs=load('data/derived/spatial_inputs_only.json');results=[]
    for diameter in sorted({r['reported_diameter_um'] for r in inputs}):
        model=build(kb,diameter,protocol);initial=None
        for row in sorted([r for r in inputs if r['reported_diameter_um']==diameter],key=lambda r:r['activation_magnitude']):
            path=f'runs/spatial_tension/{branch}/prediction_{row["record_id"]}.json'
            if (ROOT/path).exists():
                record=load(path)
                if record['hashes']!=hashes:raise RuntimeError('Prediction provenance differs')
                result=record['result']
            else:
                strength=calibration['gain']*row['activation_magnitude']
                result=checked_relax(model,strength,initial)
                record={'hashes':hashes,'input':row,'gain':calibration['gain'],'result':result}
                if result['single_valued_surface_graphs']:
                    top=np.array(result['top_profile']);bottom=np.array(result['bottom_profile'])
                    try: record['slice_observer']=observe(top[:,0],top[:,1],bottom[:,0],bottom[:,1])
                    except ValueError as error: record['slice_observer']={'status':'UNEVALUABLE','reason':str(error),'possible_labels':[0,1]}
                else:record['slice_observer']={'status':'UNEVALUABLE','reason':'nonmonotone surfaces','possible_labels':[0,1]}
                save(path,record)
            if not result['converged']:raise RuntimeError('Target solve failed; retained, never silently discarded')
            initial=np.array(result['free_displacement'])
            results.append({'record_id':row['record_id'],'path':path,'sha256':sha(path)})
            print(json.dumps({'branch':branch,'prediction':row['record_id'],'diameter_um':diameter,'force':result['scaled_free_force_linf']}),flush=True)
    if hashes!=manifest():raise RuntimeError('Inputs changed during prediction')
    save(f'runs/spatial_tension/{branch}/frozen_predictions.json',{'status':'FROZEN_CONDITIONAL_DISCRETE_GEOMETRY','hashes':hashes,'gain':calibration['gain'],'records':results,'no_heldout_response_loads':True,'scope':'Finite-scenario, finite-mesh quasistatic prediction; no specimen/continuum/transient validation claim'})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','table','calibrate','predict']);parser.add_argument('--tension',type=float,default=2.)
    args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:{'table':table,'calibrate':calibrate,'predict':predict}[args.stage](args.tension)
