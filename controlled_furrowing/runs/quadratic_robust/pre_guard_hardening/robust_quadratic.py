"""Verified-state adaptive continuation for the frozen P2 pointwise energy.

Only numerical policy changes: larger iteration budgets and rejection/bisection
of failed activation increments. Energy, gradient and geometry are inherited.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import hashlib,json,time
import numpy as np
import torch
from axisymmetric_quadratic import Parameters,AxisymmetricModel as FrozenModel


@dataclass(frozen=True)
class SolvePolicy:
    maxiter:int=1200
    newton_iterations:int=25
    minimum_step:float=1e-3
    maximum_attempts:int=128
    force_tolerance:float=1e-7

    def __post_init__(self):
        if self.maxiter<1 or self.newton_iterations<1 or self.maximum_attempts<1:raise ValueError('positive iteration/attempt budgets required')
        if not np.isfinite(self.minimum_step) or self.minimum_step<=0:raise ValueError('positive finite minimum_step required')
        if self.force_tolerance!=1e-7:raise ValueError('the force threshold is locked at1e-7')


def _vector_hash(vector):
    return hashlib.sha256(np.asarray(vector,dtype='<f8').tobytes()).hexdigest()


class AxisymmetricModel(FrozenModel):
    def __init__(self,parameters=Parameters(),policy=SolvePolicy()):
        self.policy=policy
        super().__init__(parameters)

    def _newton_polish(self,vector,strength,iterations=None):
        return super()._newton_polish(vector,strength,iterations=self.policy.newton_iterations)

    def _attempt(self,strength,initial):
        # Direct inherited fixed-strength solve; no adaptive recursion here.
        return super().solve(strength,continuation_steps=0,initial=initial,
                             continuation_start_strength=strength,maxiter=self.policy.maxiter)

    def _verify_state(self,result):
        if not result.get('converged',False):raise ValueError('previous result is not an accepted verified state')
        if result.get('parameters')!=asdict(self.parameters):raise ValueError('previous state has different model parameters')
        vector=np.asarray(result['free_displacement'],float)
        strength=float(result['activation_strength'])
        if vector.shape!=(len(self.free),) or not np.all(np.isfinite(vector)):raise ValueError('invalid previous free vector')
        _,g=self.objective(vector,strength)
        cert=self.orientation_certificate(self.coordinates(torch.tensor(vector,dtype=torch.float64)))
        if np.max(abs(g))>=self.policy.force_tolerance or not cert['certified_positive']:raise ValueError('previous state fails replayed force/geometry verification')
        return vector,strength

    def advance(self,target_strength,previous=None,output_directory=None):
        """Reach target by bisecting failed increments from the last verified state.

        Returns the last verified state with converged=False if target cannot be
        reached within declared limits. Failed state vectors are never reused.
        """
        if not np.isfinite(target_strength) or target_strength<0:raise ValueError('nonnegative finite target required')
        output=Path(output_directory) if output_directory is not None else None
        if output is not None:output.mkdir(parents=True,exist_ok=True)
        if previous is None:
            current=self._attempt(0.,np.zeros(len(self.free)))
            if not current['converged']:raise RuntimeError('zero-activation reference did not verify')
        else:current=previous
        vector,strength=self._verify_state(current)
        if target_strength<strength:raise ValueError('this policy currently supports nondecreasing activation only')
        initial_energy=self.objective(vector,strength,False)*self.energy_scale
        attempts=[];accepted=[];pending=[] if target_strength==strength else [float(target_strength)]
        stop_reason='target reached';started=time.monotonic()
        total_work=0.;total_drop=0.
        while pending:
            if len(attempts)>=self.policy.maximum_attempts:
                stop_reason='maximum_attempts reached';break
            candidate=pending.pop(0)
            before_vector=vector.copy();before_strength=strength
            before_hash=_vector_hash(before_vector)
            record={'attempt':len(attempts),'from_strength':strength,'to_strength':candidate,'initial_verified_vector_sha256':before_hash}
            tick=time.monotonic()
            try:
                result=self._attempt(candidate,before_vector)
                # The fixed solve already enforces force+geometry. Replay avoids
                # relying on stale flags in any future alternate solver backend.
                passed=bool(result['converged'])
                if passed:self._verify_state(result)
                record.update({'accepted':passed,'final_force_linf':result['scaled_free_force_linf'],'geometry_certified':result['orientation_certificate']['certified_positive'],'final_vector_sha256':_vector_hash(result['free_displacement'])})
                if output is not None:
                    path=output/f"attempt_{record['attempt']:03d}_gain_{candidate:.12g}.json"
                    path.write_text(json.dumps(result,indent=2)+'\n');record['result_path']=str(path)
            except (RuntimeError,ValueError,np.linalg.LinAlgError) as error:
                result=None;passed=False;record.update({'accepted':False,'error_type':type(error).__name__,'error':str(error)})
            record['seconds']=time.monotonic()-tick
            attempts.append(record)
            if passed:
                before_old=self.objective(before_vector,before_strength,False)*self.energy_scale
                before_new=self.objective(before_vector,candidate,False)*self.energy_scale
                work=before_new-before_old;drop=before_new-result['energy']['total']
                accepted.append({'from_strength':before_strength,'to_strength':candidate,'external_activation_potential_work':work,'relaxation_potential_drop':drop,'result_vector_sha256':record['final_vector_sha256']})
                total_work+=work;total_drop+=drop
                current=result;vector=np.asarray(result['free_displacement'],float).copy();strength=candidate
            else:
                # Neither result nor its optimizer vector enters the next solve.
                half_step=(candidate-before_strength)/2
                if half_step<self.policy.minimum_step:
                    stop_reason='failed step cannot be bisected above minimum_step';break
                pending=[before_strength+half_step,candidate]+pending
            if output is not None:
                (output/'attempt_ledger.json').write_text(json.dumps({'policy':asdict(self.policy),'requested_strength':target_strength,'last_verified_strength':strength,'attempts':attempts,'accepted':accepted},indent=2)+'\n')
        answer=dict(current);reached=bool(not pending and strength==target_strength)
        # A break after popping the last failed candidate can leave pending empty;
        # equality with target is the decisive completion test.
        reached=bool(strength==target_strength)
        answer['last_verified_state_converged']=bool(current['converged'])
        answer['converged']=bool(reached and current['converged'])
        answer['requested_activation_strength']=float(target_strength)
        answer['reached_target']=reached
        answer['robust_solver_policy']=asdict(self.policy)
        answer['robust_continuation']={'attempts':attempts,'accepted_increments':accepted,'stop_reason':stop_reason,'last_verified_strength':strength,'seconds':time.monotonic()-started,'failed_state_reuse':False}
        end_energy=answer['energy']['total']
        answer['accepted_path_work_ledger']={'initial_energy':initial_energy,'final_energy':end_energy,'external_activation_potential_work':total_work,'relaxation_potential_drop':total_drop,'balance_residual':end_energy-initial_energy-total_work+total_drop,'interpretation':'Chosen mechanical potential, accepted states only; failed numerical trials are not a physical trajectory or ATP/optical energy expenditure.'}
        if output is not None:(output/'final.json').write_text(json.dumps(answer,indent=2)+'\n')
        return answer


def build_model(parameters=None,policy=SolvePolicy(),**kwargs):
    return AxisymmetricModel(Parameters(**kwargs) if parameters is None else parameters,policy)


def _run():
    root=Path(__file__).resolve().parents[1];out=root/'runs/quadratic_robust';out.mkdir(parents=True,exist_ok=True)
    base={'nr':12,'nz':6,'height':1.,'outer_radius':3.,'patch_radius':9.35512037/33.5,'activation_z_min':8.5/33.5,'activation_z_max':28.5/33.5,'mu':1.,'basal_area_modulus':2.,'apical_area_modulus':.2,'basal_surface_tension':2.,'apical_surface_tension':0.}
    protocol={'scope':'Target-free small-geometry robust numerical continuation, no outcomes fitted','geometry_height_um':33.5,'parameters_except_lambda':base,'lambda':[500.,5000.],'targets':[1.,5.,10.],'policy':asdict(SolvePolicy()),'verification':['Energy and gradient bit-for-bit against frozen fullP2 on random states','Every accepted state replayed force<1e-7 and sufficient positiveJ','Every failed attempt retained; next attempt initial hash must be last accepted hash','Bounds are not attained compression; no continuum-error claim']}
    (root/'evidence/quadratic_robust_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    rng=np.random.default_rng(19111989);model=build_model(Parameters(bulk=500.,**base));frozen=FrozenModel(Parameters(bulk=500.,**base));equality=[]
    for gain in [0.,1.,5.]:
        vector=rng.normal(size=len(model.free))*1e-5;e,g=model.objective(vector,gain);ef,gf=frozen.objective(vector,gain)
        equality.append({'gain':gain,'energy_exactly_equal':e==ef,'gradient_exactly_equal':bool(np.array_equal(g,gf))})
    rows=[]
    for lam in protocol['lambda']:
        model=build_model(Parameters(bulk=lam,**base));previous=None
        for target in protocol['targets']:
            result=model.advance(target,previous,out/f'lambda_{lam}_target_{target}')
            row={'lambda':lam,'target':target,'reached':result['reached_target'],'verified_strength':result['activation_strength'],'force':result['scaled_free_force_linf'],'sampled_J_min':result['sampled_minimum_J'],'sampled_J_max':result['sampled_maximum_J'],'shortening':result['maximum_activated_column_shortening'],'attempts':len(result['robust_continuation']['attempts']),'failed_attempts':sum(not q['accepted'] for q in result['robust_continuation']['attempts']),'seconds':result['robust_continuation']['seconds'],'stop_reason':result['robust_continuation']['stop_reason']}
            rows.append(row);print(json.dumps(row),flush=True)
            receipt={'status':'COMPLETE' if len(rows)==6 and all(q['reached'] for q in rows) else 'IN_PROGRESS_OR_STOPPED','energy_gradient_equality':equality,'all_energy_gradients_bit_identical':all(q['energy_exactly_equal'] and q['gradient_exactly_equal'] for q in equality),'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parents':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['axisymmetric_mechanics.py','axisymmetric_prestress.py','axisymmetric_mixed.py','axisymmetric_quadratic.py']},'scope':protocol['scope']}
            (root/'evidence/quadratic_robust_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
            if not result['reached_target']:break
            previous=result


if __name__=='__main__':
    torch.set_num_threads(1)
    _run()
