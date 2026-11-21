"""Score already frozen conditional spatial predictions; never refit them."""
import csv
import hashlib
import json
from pathlib import Path
import math
import numpy as np
from spatial_observer import observe

ROOT=Path(__file__).resolve().parents[1]


def load(relative):
    return json.loads((ROOT/relative).read_text(encoding='utf-8'))


def sha(relative):
    return hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()


def csv_rows(relative):
    with (ROOT/relative).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    protocol=load('SPATIAL-PROTOCOL.json');height=protocol['units']['height_um']
    frozen=[]
    for kb in protocol['material_scenarios']['basal_area_modulus_over_mu_height']:
        branch=f'basal_{kb:g}';path=f'runs/spatial_transfer/{branch}/frozen_predictions.json'
        prediction=load(path)
        for source,expected in prediction['hashes'].items():
            if sha(source)!=expected:raise RuntimeError('Frozen source changed: '+source)
        records=[]
        for entry in prediction['records']:
            if sha(entry['path'])!=entry['sha256']:raise RuntimeError('Frozen result changed')
            records.append(load(entry['path']))
        if len(records)!=31:raise RuntimeError('Incomplete larger-region prediction set')
        frozen.append((branch,kb,path,prediction,records))
    # Larger responses are loaded only after every branch is frozen and checked.
    evaluation=csv_rows('data/derived/evaluation_larger_bending.csv')
    labels={r['record_id']:int(r['bends']) for r in evaluation}
    scalar=csv_rows('data/derived/evaluation_large_50um.csv')
    results=[]
    for branch,kb,path,prediction,records in frozen:
        observations=[]
        # These are conditional error assumptions, not estimated FE error bars.
        epsilon_um=[0.,.05,.2]
        excluded={str(e):0 for e in epsilon_um};ambiguous={str(e):0 for e in epsilon_um}
        for record in records:
            result=record['result'];row=record['input'];label=labels[row['record_id']]
            item={'record_id':row['record_id'],'observed_label':label,'input':row,'observer':record['slice_observer'],
                  'central_apical_displacement_um':height*result['central_apical_displacement'],
                  'central_basal_displacement_um':height*result['central_basal_displacement'],
                  'maximum_column_shortening_um':height*result['maximum_activated_column_shortening']}
            for e in epsilon_um:
                if result['single_valued_surface_graphs']:
                    top=np.asarray(result['top_profile']);bottom=np.asarray(result['bottom_profile'])
                    try:possible=observe(top[:,0],top[:,1],bottom[:,0],bottom[:,1],uniform_surface_error=e/height)['unknown_slice_label_set']
                    except ValueError:possible=[0,1]
                else:possible=[0,1]
                excluded[str(e)]+=label not in possible
                ambiguous[str(e)]+=len(possible)>1
            observations.append(item)
        pairs=[]
        for observed in scalar:
            matched=[r for r in records if abs(r['input']['activation_magnitude']-float(observed['activation_magnitude']))<1e-12 and abs(r['input']['reported_diameter_um']-float(observed['reported_diameter_um']))<1e-8]
            if len(matched)!=1:raise RuntimeError('No unique source-input match for scalar evaluation')
            predicted=height*matched[0]['result']['maximum_activated_column_shortening']
            pairs.append({'record_id':observed['record_id'],'prediction_um':predicted,'observed_um':float(observed['maximum_shortening_um']),'prediction_source':matched[0]['input']['record_id']})
        residuals=[r['prediction_um']-r['observed_um'] for r in pairs]
        stats={'RMSE_um':math.sqrt(math.fsum(r*r for r in residuals)/len(residuals)),'MAE_um':math.fsum(abs(r) for r in residuals)/len(residuals),'signed_bias_um':math.fsum(residuals)/len(residuals)}
        results.append({'branch':branch,'basal_modulus':kb,'gain':prediction['gain'],'frozen_prediction_sha256':sha(path),
                        'scalar_evaluation_n':len(pairs),'scalar_scores':stats,'scalar_pairs':pairs,
                        'label_n':len(records),'conditional_nominal_label_exclusions_by_assumed_surface_error_um':excluded,
                        'ambiguous_label_sets_by_assumed_surface_error_um':ambiguous,
                        'maximum_force_residual':max(r['result']['scaled_free_force_linf'] for r in records),
                        'minimum_local_J':min(r['result']['minimum_J'] for r in records),
                        'maximum_local_J':max(r['result']['maximum_J'] for r in records),
                        'maximum_absolute_global_volume_change':max(abs(r['result']['relative_volume_change']) for r in records),
                        'records':observations})
    receipt={'status':'CONDITIONAL_DISCRETE_DIAGNOSTIC_SCORED','branches':results,
             'evidence':{source:sha(source) for source in ['SPATIAL-PROTOCOL.json','scripts/score_spatial.py','data/derived/evaluation_large_50um.csv','data/derived/evaluation_larger_bending.csv']},
             'error_assumption_warning':'0/.05/.2 um are illustrative assumed uniform surface-error bounds, NOT measured or certified errors. Label exclusions remain conditional on slice/body/time assumptions.',
             'missing_validation':['specimen initial geometry and material bounds','complete calibration identified set','source transient/slice observation operator','rigorous continuum accuracy','ellipse/non-axisymmetric response'],
             'novelty':'None established; this is model reproduction/stress testing and conditional observation mathematics'}
    (ROOT/'evidence/spatial_transfer_scores.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    plot(results)
    print(json.dumps({'status':receipt['status'],'branches':[{k:r[k] for k in ['branch','gain','scalar_scores','conditional_nominal_label_exclusions_by_assumed_surface_error_um','minimum_local_J','maximum_force_residual']} for r in results]}))


def plot(results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('dark_background')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    for r,color in zip(results,['#8ad4ef','#e8bd80','#bfa4e9']):
        x=[p['observed_um'] for p in r['scalar_pairs']];y=[p['prediction_um'] for p in r['scalar_pairs']]
        axes[0].scatter(x,y,color=color,label=f'Basal modulus {r["basal_modulus"]:g}')
        axes[1].plot([0,.05,.2],[r['conditional_nominal_label_exclusions_by_assumed_surface_error_um'][str(e)] for e in [0.,.05,.2]],'o-',color=color,label=f'Basal modulus {r["basal_modulus"]:g}')
    axes[0].plot([0,30],[0,30],'--',color='#aaaaaa')
    axes[0].set(xlabel='Observed maximum shortening (µm)',ylabel='Predicted material-column maximum (µm)',title='Larger-region scalar diagnostic; proxy mismatch explicit')
    axes[0].legend(fontsize=8)
    axes[1].set(xlabel='ASSUMED uniform surface error (µm)',ylabel='Observed labels excluded by possible-label set',title='Conditional observation sensitivity, not error certification')
    axes[1].legend(fontsize=8)
    fig.suptitle('Finite-height discrete model: conditional stress test, no novelty or continuum-validation claim',fontsize=12)
    fig.savefig(ROOT/'docs/assets/spatial_transfer.png',dpi=170)
    plt.close(fig)


if __name__=='__main__':
    main()
