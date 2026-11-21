"""Cold audit of the four-target exact symmetric-preconditioner caller."""

from base64 import b64decode
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch
from zipfile import ZipFile

import numpy as np


HERE = Path(__file__).resolve().parent
PREP = HERE.parent
SPD = PREP.parent
EXACT = SPD.parent
ROOT = EXACT.parents[3]
SOURCE = ROOT/'scripts/transfer_ambiguity_preconditioner_spd.py'
MAIN = '67ef2be40fc423d721390fd360f6244ecfc529923aa5c41a04981edec4cd20d7'
PROTOCOL = '1b6369802375d358abf14d8db75ae1b0e3d418be8705c2040395ce2cc434c905'
PREFLIGHT = '489cce78693ac12b6a0851c48b39895a0a65022611de356367b0b5275a2846cc'
RANK = 611


def require(test,message):
    if not test:
        raise ValueError(message)


def read(path, digest):
    raw = path.read_bytes()
    require(sha256(raw).hexdigest()==digest,f'failed exact-byte pin {path}')
    return json.loads(raw)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def checked_path(name):
    path=(ROOT/name).resolve()
    require(path.is_relative_to(ROOT) and path.relative_to(ROOT).as_posix()==name,'noncanonical authority path')
    return path


def npy_identity(packet):
    block=packet['frozen_binary_cache']
    payload=b64decode(block['data'],validate=True)
    require(sha256(payload).hexdigest()==block['sha256'],'changed embedded cache bytes')
    with ZipFile(BytesIO(payload)) as z, z.open('preconditioner.npy') as f:
        require(np.lib.format.read_magic(f)==(1,0),'unexpected frozen NPY format')
        shape,fortran,dtype=np.lib.format.read_array_header_1_0(f)
        raw=f.read()
    require(shape==(RANK,RANK) and dtype.str=='<f8' and len(raw)==8*RANK*RANK,'changed frozen R format')
    # Canonical logical C-order SHA, even when the archive stores column-major entries.
    if fortran:
        digest=sha256()
        for i in range(RANK):
            digest.update(b''.join(raw[offset:offset+8] for offset in range(8*i,len(raw),8*RANK)))
        return digest.hexdigest(), 'Fortran'
    return sha256(raw).hexdigest(), 'C'


def verify_inputs(protocol):
    links=protocol['source_and_input_sha256']
    for name,h in links.items():
        require(sha(checked_path(name))==h,'changed input closure')
    catalog=read(EXACT/'preparation/target_manifest.json',protocol['target_manifest_sha256'])
    rows=protocol['ordered_targets']
    require(tuple(r['pilot_id'] for r in rows)==tuple(catalog['pilot_order']) and len(rows)==4
            and len({r['pilot_id'] for r in rows})==4,'exactly four distinct admitted targets')
    checked=[]
    for row,entry in zip(rows,catalog['targets'],strict=True):
        ident=row['identity']
        contract=read(checked_path(entry['model_contract']),entry['model_contract_sha256'])
        context_row=read(checked_path(entry['row']),entry['row_sha256'])
        packet=read(checked_path(contract['point_cache']),ident['cache_packet_sha256'])
        prior=read(checked_path(contract['cached_root_independent_receipt']),row['point_proof_sha256'])
        h,storage=npy_identity(packet)
        cell=contract['parameter_cell']
        cell_hash=sha256(json.dumps(cell,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        p=contract['exact_requested_drive']
        scale=Q(contract['energy_scale_binary_exact'])
        req=Q(p['activation_strength_exact'])
        dose=Q(p['optical_dose_exact'])
        gain=Q(p['common_gain_exact'])
        require(0<req==dose*gain==Q(ident['drive_exact'])==Q(packet['activation_strength_exact'])==Q(prior['exact_strength'])
                and req!=Q(packet['activation_strength_binary_exact']) and 0<scale==Q(packet['normalization_binary_exact'])
                ==Q(ident['energy_scale_binary_exact']) and (dose,gain)==(Q(ident['dose_exact']),Q(ident['gain_exact']))
                and h==contract['binary_preconditioner_sha256']==ident['R_binary_sha256']
                and packet['frozen_binary_cache']['sha256']==ident['embedded_cache_sha256']==prior['cache_array_contents_sha256']
                and prior['cache_packet_sha256']==ident['cache_packet_sha256']==context_row['point_cache_sha256']
                and contract['source_field_sha256']==ident['field_sha256']==packet['input_state_sha256']==prior['field_sha256']
                and ident['pilot_id']==entry['pilot_id']==contract['pilot_id']==prior['pilot_id']
                and ident['record_id']==contract['record_id']==context_row['record_id']
                and ident['parameter_cell_sha256']==cell_hash==context_row['parameter_cell_sha256']
                and all(Q(lo)==Q(hi) for lo,hi in cell.values()) and tuple(map(Q,cell['common_gain']))==(gain,gain)
                and ident['model_contract_sha256']==entry['model_contract_sha256']
                and row['larger_radius_box_sha256']==entry['box_1_sha256']
                and 0<Q(row['cached_reference_R_norm_upper_exact'])==Q(prior['preconditioner_infinity_norm_exact']),
                'mixed model, exact drive, matrix, energy scale, proof, or cell')
        checked.append({'pilot_id':ident['pilot_id'],'record_id':ident['record_id'],'binary_R_sha256':h,
                        'disk_array_storage_order':storage,'verified_canonical_order':'C','cache_sha256':ident['cache_packet_sha256'],
                        'dose_exact':str(dose),'gain_exact':str(gain),'drive_exact':str(req),'energy_scale_exact':str(scale),
                        'cell_sha256':cell_hash})
    return checked


def exact_oracle(a,k):
    """Small exact, direct six-factor expansion of K sym(R) K^T."""
    n=len(a)
    r=[[Q.from_float(float(x)) for x in row] for row in a]
    t=[[Q.from_float(float(x)) for x in row] for row in k]
    B=[[sum(t[i][l]*(r[l][m]+r[m][l])*t[j][m]/2 for l in range(n) for m in range(n))
        for j in range(n)] for i in range(n)]
    q=max(sum(abs(Q(i==j)-B[i][j]) for j in range(n)) for i in range(n))
    return q,max(sum(map(abs,row)) for row in r)


def rejected(f):
    try:
        f()
    except ValueError:
        return True
    return False


def numeric_fixtures(c):
    a=np.array([[2.,1.],[2.**-53,2.]],dtype='<f8')
    k=np.eye(2,dtype='<f8')/2
    exact,expected_norm=exact_oracle(a,k)
    require(exact==Q(5,8)+Q(1,2**56),'analytic binary oracle inaccurate')
    prod,_=c.producer.certify_sympart(a,k,640)
    q,rows=c.replay.congruence_bounds(a,k,640)
    qp=Q(prod['q_infinity_upper_exact'])
    require(Q(5,8)<exact<=min(qp,q)<=max(qp,q)<1 and max(rows)==q
            and prod['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD'
            and c.replay.norm_R_integer(a)==expected_norm==Q(3),
            'a rounded sym(R) could falsely match the outward oracle')
    improper=np.array([[1.,4.],[0.,1.]],dtype='<f8')
    I=np.eye(2,dtype='<f8')
    invalid,_=c.producer.certify_sympart(improper,I,640)
    bad,_=c.replay.congruence_bounds(improper,I,640)
    singular=np.diag([1.,0.]).astype('<f8')
    boundary,_=c.producer.certify_sympart(I,singular,640)
    qb,_=c.replay.congruence_bounds(I,singular,640)
    require(exact_oracle(improper,I)[0]==bad==Q(invalid['q_infinity_upper_exact'])==2
            and invalid['status']!='CERTIFIED_EXACT_SYMMETRIC_PART_SPD'
            and qb==Q(boundary['q_infinity_upper_exact'])==1
            and boundary['status']!='CERTIFIED_EXACT_SYMMETRIC_PART_SPD','strict condition or negative oracle failed')
    huge=np.array([[np.finfo(np.float64).max]],dtype='<f8')
    failed,proposal=c.producer.certify_sympart(huge,precision_bits=640)
    require(failed['status']=='UNRESOLVED_FLOATING_PROPOSAL' and proposal is None
            and failed['positive_or_negative_eigenvalue_inference'] is False,
            'a floating-point proposal failure must remain inconclusive')
    return {'rounded_proxy_q_exact':'5/8','exact_dyadic_oracle_q':str(exact),'producer_q':str(qp),'alternate_q':str(q),
            'indefinite_q':str(bad),'singular_K_border_q':str(qb),'exact_norm_R':str(expected_norm)}


def factory_fixture(c):
    pilot=c.PILOTS[0]
    a=np.array([[2.,.25],[-.125,3.]],dtype='<f8')
    a.flags.writeable=False
    signature=sha256(a.tobytes(order='C')).hexdigest()
    row={'pilot_id':pilot,'identity':{'pilot_id':pilot,'record_id':'fixture','field_sha256':'field',
          'cache_packet_sha256':'cache','embedded_cache_sha256':'embedded','R_binary_sha256':signature,
          'model_contract_sha256':'contract','parameter_cell_sha256':'cell','dose_exact':'2/7',
          'gain_exact':'3/5','drive_exact':'6/35','energy_scale_binary_exact':'3/2'}}
    target=SimpleNamespace(packet={'input_state_sha256':'field'},energy_scale=Q(3,2))
    fake=SimpleNamespace(R=a,u=np.zeros(2),free=np.arange(2),target=target,pilot_id=pilot,
        record_id='fixture',row={'point_cache_sha256':'cache','frozen_binary_cache_sha256':'embedded',
                                 'parameter_cell_sha256':'cell'},dose=Q(2,7),base={'common_gain':Q(3,5)},drive=Q(6,35))
    events=[]
    def boot():
        events.append('bootstrap')
        c.EX.bootstrap()
    def open_context(*_):
        require(events==['bootstrap'] and c.arb_context.prec==640 and c.EX.T.core.A.cache_info().currsize==0,
                'fresh shared Arb context must precede the factory')
        events.append('factory')
        return fake,None,{'model_contract_sha256':'contract'}
    orig=c.find_matrix
    c.arb_context.prec=128
    c.EX.T.core.A(Q(1,7))
    require(c.EX.T.core.A.cache_info().currsize>0,'synthetic cache was not populated')
    with patch.object(c,'NFREE',2),patch.object(c,'boot',side_effect=boot), \
         patch.object(c.EX,'open_context',side_effect=open_context), \
         patch.object(c,'find_matrix',side_effect=lambda m:orig(m,2)):
        matrix,data=c.actual_matrix(pilot,(None,[row],None,None))
    require(matrix is a and data==row and events==['bootstrap','factory'],'fixture matrix routing changed')
    return True


def main():
    require(sha(SOURCE)==sha(PREP/'transfer_ambiguity_preconditioner_spd_frozen.py')==MAIN,
            'current and frozen SPD caller differ')
    protocol=read(PREP/'protocol.json',PROTOCOL)
    pf=read(PREP/'preflight.json',PREFLIGHT)
    require(pf['status']=='PASS' and all(v is True for v in pf['checks'].values())
            and type(pf['actual_production_matrix_jobs']) is int and pf['actual_production_matrix_jobs']==0
            and pf['source_sha256']==MAIN and pf['protocol_sha256']==PROTOCOL,
            'changed cold preflight')
    identities=verify_inputs(protocol)
    sys.path.insert(0,str(ROOT/'scripts'))
    import transfer_ambiguity_preconditioner_spd as c
    require(Path(c.__file__).resolve()==SOURCE,'explicit reviewed module required')
    oracle=numeric_fixtures(c)
    factory_fixture(c)
    stop=AssertionError('attempted production matrix or AD')
    with patch.object(c.EX,'open_context',side_effect=stop), \
         patch.object(c.producer,'certify_sympart',side_effect=stop), \
         patch.object(c.replay,'congruence_bounds',side_effect=stop), \
         patch.object(c.EX.T.core,'assemble_point',side_effect=stop), \
         patch.object(c.EX.T.core,'local_terms',side_effect=stop):
        snapshot=c.authority()
        require(c.protocol(snapshot)==protocol,'live cold contract differs')
        bindings={'source_sha256':MAIN,'protocol_sha256':PROTOCOL,'preflight_sha256':PREFLIGHT,
                  'minimum_design_sha256':protocol['minimum_design_sha256'],
                  'cached_root_review_sha256':protocol['cached_root_review_sha256'],
                  'target_manifest_sha256':protocol['target_manifest_sha256']}
        fakepath=HERE/'unused_synthetic_review.json'
        real_pinned=c.pinned
        synthetic={'status':c.REVIEW_STATUS,'checks':{'fixture':True},**bindings}
        def substitute(path,digest):
            return synthetic if path==fakepath else real_pinned(path,digest)
        with patch.object(c,'pinned',side_effect=substitute),patch.object(c,'authority',return_value=snapshot):
            _,admitted=c.gate(fakepath,'7'*64)
            require(admitted=={**bindings,'prelaunch_review_sha256':'7'*64},'review gate status or binding changed')
            synthetic['target_manifest_sha256']='0'*64
            require(rejected(lambda:c.gate(fakepath,'7'*64)),'cross-manifest approval accepted')
    identity={'pilot_id':c.PILOTS[0],'pid':7,'supervisor_pid':3,'bindings':bindings}
    observed={**identity,'observed_exit_code':0}
    terminal={**identity,'status':'FINISHED','exit_code':0,'result_sha256':'a'*64}
    c.check_child(observed,terminal,identity,'a'*64)
    require(all(rejected(f) for f in (
        lambda:c.check_child({**observed,'observed_exit_code':False},terminal,identity,'a'*64),
        lambda:c.check_child({**observed,'pid':7.0},terminal,identity,'a'*64),
        lambda:c.check_child(observed,{**terminal,'exit_code':False},identity,'a'*64),
        lambda:c.check_child(observed,{**terminal,'supervisor_pid':4},identity,'a'*64),
        lambda:c.check_child(observed,terminal,identity,'b'*64))),
        'OS identity, exit or exact result pin was not enforced')
    terminal_pid_alias = not rejected(lambda:c.check_child(observed,{**terminal,'pid':7.0},identity,'a'*64))
    require(not c.ACTUAL.exists(),'production attempt already exists during cold review')
    result={'status':c.REVIEW_STATUS,**bindings,
            'checks':{'all_current_inputs_frozen_and_four_binary_R_identities_independently_verified':True,
                      'per_target_energy_exact_drive_cache_and_degenerate_cell_verified':True,
                      'cold_authority_matches_preparation_and_refuses_actual_work':True,
                      'independent_exact_dyadic_oracle_two_congruences_and_R_norm':True,
                      'reject_indefinite_and_singular_K_equality_and_float_proposal_inconclusive':True,
                      'explicit_matrix_factory_bootstrap_and_cache_clear_verified_by_synthetic_route':True,
                      'separate_frozen_source_protocol_preflight_and_manifest_gate_tested':True,
                      'supervisor_records_real_child_wait_and_exclusive_single_run_and_worker_outputs':True,
                      'reject_nonliteral_observed_exit_or_worker_mismatch_and_bind_result_hash':True},
            'review_source_sha256':sha(Path(__file__)),'reviewed_current_input_count':len(protocol['source_and_input_sha256']),
            'independently_matched_targets':identities,'synthetic_matrix_readback':oracle,'actual_production_matrix_jobs':0,
            'actual_force_Hessian_AD_or_solve_jobs':0,'matrix_scope':'1- and 2-dimensional synthetic congruences.',
            'nonblocking_schema_observation':{'accepts_terminal_numeric_pid_alias_under_equal_typed_OS_PID':terminal_pid_alias,
                'detail':'The frozen terminal field is produced from os.getpid() and verified against a typed real observed PID. Later actual independent review should also require literal integer terminal and started PIDs.'},
            'scope':'Cold prelaunch approval for four exact symmetric-part sufficient tests. The source uses two different product associations and the same frozen binary R and K. This approval establishes no actual SPD, root, minimum, stability or continuum result.'}
    path=HERE/'receipt.json'
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({'status':result['status'],'sha256':sha(path),'targets':len(identities),'checked_inputs':result['reviewed_current_input_count']}))


if __name__=='__main__':
    main()
