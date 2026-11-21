"""Explicit-record exact-integration pilot, with separate reviewed phases.

The first attempt uses the original predictor cube and its audited cached-Q
bounds. A failed sufficient test remains unresolved. No parameter fitting,
new equilibrium solve, global identified-set or continuum claim is made.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,datetime,hashlib,json,os,time,traceback
from flint import arb,arb_mat,ctx as precision
import exact_cohort_context as C
import exact_cohort_integrator as E
import root_body_geometry as G
S=E.SUR;core=E.core
ROOT=Path(__file__).resolve().parents[1];SOURCE=Path(__file__)
BASE=ROOT/'runs/exact_integration_cohort_v1/records'
PRECISION=640;MARGIN_DIVISOR=10000

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def serial(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):serial(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [serial(x)for x in v]
 return v
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():raise ValueError('preserve existing output '+str(p))
 p.write_text(json.dumps(serial(d),indent=2,allow_nan=False)+'\n',encoding='utf8')
def alltrue(v):return isinstance(v,dict)and bool(v)and all(x is True for x in v.values())
def outdir(record):return BASE/record/'original_cube_v1'
def assert_inputs(p):
 for name,h in p['source_and_input_hashes'].items():
  if sha(ROOT/name)!=h:raise ValueError('changed frozen input '+name)
def bind(d,path,expected=None):
 path=Path(path).resolve();key=path.relative_to(ROOT.resolve()).as_posix();h=sha(path)
 if(expected is not None and h!=expected)or(key in d and d[key]!=h):raise ValueError('conflicting authoritative input '+key)
 d[key]=h
def gate(out,phase):
 r=read(out/(phase+'_prelaunch_independent_review.json'));p=out/'protocol.json';pre=out/'preflight.json'
 if r.get('status')!='PASS'or not alltrue(r.get('checks'))or r.get('source_sha256')!=sha(SOURCE)or r.get('protocol_sha256')!=sha(p)or r.get('preflight_sha256')!=sha(pre):raise ValueError('bound independent '+phase+' review required')
 if phase=='integration'and r.get('domains_result_sha256')!=sha(out/'domains.json'):raise ValueError('current reviewed geometry/selectors required')

def incidence(context):
 count=[0]*611;rowcount=[0]*611
 groups=[list(map(int,t))for t in context.target.a['triangles']]
 groups += [list(map(int,t))for k in(0,1)for t in context.target.a[f'surface_{k}_edges']]
 for nodes in groups:
  ids=[context.target.index[2*n+k]for n in nodes for k in(0,1)if 2*n+k in context.target.index]
  for i in ids:count[i]+=1;rowcount[i]+=len(ids)
 return count,rowcount

def prepare(record):
 context=C.load_context(record);out=outdir(record)
 if(out/'protocol.json').exists():raise ValueError('preserve prior preparation')
 inputs=dict(context.input_hashes)
 expected={ROOT/'evidence/exact_cohort_context_independent_review.json':'898fb6828c16fc7f0db410c9d677373cf1696db5d0d418c08ef0004a30e215c5',ROOT/'evidence/exact_cohort_integrator_independent_review.json':'1d233c61e5b7283ddfeb856e206a064d28cb2731246f13bff716e9a4063dbfa8'}
 for path,h in expected.items():
  bind(inputs,path,h);v=read(path)
  if not v['status'].startswith('PASS')or not alltrue(v['checks']):raise ValueError('nonempty generic review required')
 deps=read(ROOT/'runs/exact_integration_cohort_v1/integrator_preparation/protocol.json')['source_hashes_and_inputs']
 for path,h in deps.items():bind(inputs,ROOT/path,h)
 for path in [SOURCE,Path(C.__file__),Path(E.__file__),Path(G.__file__),Path(S.__file__)]:bind(inputs,path)
 capsule=read(ROOT/context.row['independent_root_capsule_file']);audit=next(v for v in capsule['record_audits']if v['record_id']==record)
 if not alltrue(audit['audit']['checks']):raise ValueError('audited cached root required')
 bind(inputs,ROOT/audit['review_file'],audit['review_sha256']);old=read(ROOT/context.row['root_certificate_file'])
 eta=max(F(audit['audit']['bounds']['eta']['upper_exact']),F(old['uniform_cell_result']['eta_uniform']['upper_exact']))
 z=F(audit['audit']['bounds']['Z']['upper_exact'])
 if not 0<=z<1:raise ValueError('audited old contraction required')
 Rnorm=max(sum((abs(F.from_float(float(v)))for v in row),F(0))for row in context.R)
 counts,rows=incidence(context);fc,hc=max(counts),max(rows)
 if min(fc,hc,Rnorm)<=0:raise ValueError('nonempty finite positive incidence/norm required')
 fbudget=(1-z)*context.old_radius/(MARGIN_DIVISOR*Rnorm*fc)
 hbudget=(1-z)/(MARGIN_DIVISOR*Rnorm*hc)
 box=serial(C.exact_node_box(context,context.old_radius));box['domain_id']='exactFE:'+record+':original_cube_v1'
 write(out/'coordinate_box.json',box);bind(inputs,out/'coordinate_box.json')
 protocol={'utc':now(),'source_and_input_hashes':inputs,'record_id':record,'row':context.row,'precision_bits':PRECISION,'normal_radius_exact':str(context.old_radius),'parameter_cell_sha256':context.row['parameter_cell_sha256'],
 'coordinate_box_file':(out/'coordinate_box.json').relative_to(ROOT).as_posix(),'coordinate_box_sha256':sha(out/'coordinate_box.json'),
 'cached_eta_upper_exact':str(eta),'cached_Z_upper_exact':str(z),'R_infinity_norm_exact':str(Rnorm),
 'force_incidence':counts,'H_row_incidence':rows,'max_force_incidence':fc,'max_H_row_incidence':hc,
 'force_tail_budget_exact':str(fbudget),'H_entry_tail_budget_exact':str(hbudget),'margin_divisor':MARGIN_DIVISOR,
 'global_force_tail_budget_exact':str(Rnorm*fc*fbudget),'global_derivative_tail_budget_exact':str(Rnorm*hc*hbudget),
 'tail_policy':'Each local force tail <=(1-oldZ)*oldradius/(10000*||R||*max_F_incidence); each local H-entry tail <=(1-oldZ)/(10000*||R||*max_H_row_incidence). Signed assembly then norm; actual tests decide sufficiency.',
 'max_bulk_order':32,'max_surface_order':64,'surface_domain_depth':8,
 'target_scope':'Entire explicit independent nodal box times SAME C; candidate root cube remains original p0(theta)+B(oldradius). All known cached bounds reused only on that exact original cube.',
 'no_new_solve_or_material_fit':True,'root_asserted':False}
 write(out/'protocol.json',protocol)
 checks={'same640_precision_policy':PRECISION==E.PRECISION_BITS,'all650_coordinates':len(box['current_centers_exact'])==325,'original_cube_no_shift':all(F(x)==0 for x in box['shift_exact']),
 'incidence_tail_accounting_exact':Rnorm*fc*fbudget==(1-z)*context.old_radius/MARGIN_DIVISOR and Rnorm*hc*hbudget==(1-z)/MARGIN_DIVISOR,
 'root_not_assumed':box['root_membership_assumed']is False}
 if not alltrue(checks):raise AssertionError(checks)
 write(out/'preflight.json',{'status':'PASS_PREPARATION_ONLY','checks':checks,'source_sha256':sha(SOURCE),'protocol_sha256':sha(out/'protocol.json'),'actual_body_or_force_job_run':False})
 print(json.dumps({'phase':'prepared','record':record,'source_sha256':sha(SOURCE),'protocol_sha256':sha(out/'protocol.json'),'preflight_sha256':sha(out/'preflight.json'),'force_tail':float(fbudget),'H_tail':float(hbudget)}),flush=True)

def surface_selector(context,fields,domain,rule,face,edge,p):
 if domain['status']!='CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN':raise ValueError('positive closed-edge domain unresolved')
 jet=E.JET.surface_squared_jacobian_jet(fields['q_physical_R'],fields['B_physical_R']);_,_,K,T=S.material_boxes(context.row,face)
 dmax=max(abs(t-k)for t in T for k in K);q0=F(domain['q0_exact']);rho=F(domain['rho_exact'])
 sqrtlo=F(str(core.A(q0).sqrt().lower().fmpq()))
 if sqrtlo<=0 or sqrtlo*sqrtlo>q0:raise ValueError('outward positive square-root lower bound required')
 mass=F(str(core.upper(context.target.normalization()*core.A((fields['Rhi']**2-fields['Rlo']**2)/2)+core.A(sum(rule['surface_weights'][face][edge],F(0)))).fmpq()))
 g=max(S.poly_norm_upper(v,fields['Rlo'],fields['Rhi'])for v in jet['gradient']);h=max(S.poly_norm_upper(v,fields['Rlo'],fields['Rhi'])for v in jet['hessian_upper'].values())
 attempts=[]
 for n in range(2,p['max_surface_order']+1):
  a=F(dmax,2)*E.SUR.S.binomial_tail(F(-1,2),n,rho)/sqrtlo;b=F(dmax,4)*E.SUR.S.binomial_tail(F(-3,2),n,rho)/(q0*sqrtlo)
  ft,ht=mass*a*g,mass*(a*h+b*g*g);attempts.append({'degree':n,'force_tail_upper_exact':str(ft),'H_entry_tail_upper_exact':str(ht)})
  if ft<=F(p['force_tail_budget_exact'])and ht<=F(p['H_entry_tail_budget_exact']):return n,attempts
 raise ValueError('surface order policy unresolved')

def domains(record):
 out=outdir(record);gate(out,'domains');p=read(out/'protocol.json');assert_inputs(p)
 if(out/'domains_started.json').exists():raise ValueError('single domain run only')
 write(out/'domains_started.json',{'pid':os.getpid(),'utc':now()});started=time.monotonic();precision.prec=PRECISION;core.A.cache_clear()
 context=C.load_context(record);box=read(ROOT/p['coordinate_box_file']);Y=S.common_nodes(context.target,box);rule=context.target.auxiliary_rule()
 body=G.certify_box(box['reference_nodes_exact'],box['current_centers_exact'],box['triangles'],box['coordinate_radii_exact'])
 write(out/'body_geometry.json',body)
 if body['status']!='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING':raise ValueError('whole-body geometry unresolved')
 selectors=[E.prepare_bulk_selector(context,Y,rule,body['triangles'][i],i,F(p['force_tail_budget_exact']),F(p['H_entry_tail_budget_exact']),max_order=p['max_bulk_order'])for i in range(144)]
 cases=[]
 for face in(0,1):
  for edge in range(12):
   fields=context.target.surface(face,edge,nodes=Y);domain=S.uniform_G_bounds(fields,max_depth=p['surface_domain_depth']);n,attempts=surface_selector(context,fields,domain,rule,face,edge,p)
   sites=[S.local_site_domain(context.target,box,face,edge,iq)for iq in range(12)]
   cases.append({'record_id':record,'face':face,'edge':edge,'global_indices':fields['global_indices'],'domain':domain,'order':n,'order_attempts':attempts,'cache_site_domains':sites})
 assert_inputs(p)
 result={'status':'PASS_FRESH_DOMAINS_AND_SELECTORS_UNREVIEWED','record_id':record,'protocol_sha256':sha(out/'protocol.json'),'coordinate_box_sha256':sha(ROOT/p['coordinate_box_file']),'body_sha256':sha(out/'body_geometry.json'),'bulk_selectors':selectors,'surface_cases':cases,'actual_force_H_integrals_run':False,'seconds':time.monotonic()-started}
 write(out/'domains.json',result);write(out/'domains_terminal.json',{'status':'FINISHED','exit_code':0,'pid':os.getpid(),'utc':now(),'results_sha256':sha(out/'domains.json')})
 print(json.dumps({'phase':'domains_complete','record':record,'body_J':body['J_bound'],'bulk_orders':sorted(set(v['degree']for v in selectors)),'surface_orders':sorted(set(v['order']for v in cases)),'seconds':result['seconds'],'result_sha256':sha(out/'domains.json')}),flush=True)

def integrals(record,part):
 out=outdir(record);gate(out,'integration');p=read(out/'protocol.json');assert_inputs(p);d=read(out/'domains.json')
 if sha(out/'body_geometry.json')!=d['body_sha256']or sha(ROOT/p['coordinate_box_file'])!=d['coordinate_box_sha256']:raise ValueError('changed reviewed body/box')
 if(out/(part+'_started.json')).exists():raise ValueError('single part job only')
 write(out/(part+'_started.json'),{'pid':os.getpid(),'utc':now(),'domains_sha256':sha(out/'domains.json')});precision.prec=PRECISION;core.A.cache_clear();started=time.monotonic()
 context=C.load_context(record);Y=S.common_nodes(context.target,read(ROOT/p['coordinate_box_file']));rule=context.target.auxiliary_rule();body=read(out/'body_geometry.json');files=[]
 log=(out/(part+'_stdout.log')).open('x',encoding='utf8')
 try:
  for index,case in enumerate(d['bulk_selectors']if part=='bulk'else d['surface_cases']):
   if part=='bulk':
    b=body['triangles'][case['element']]['J_bound']
    if case['J_bounds_exact']!=[b['lower_exact'],b['upper_exact']]:raise ValueError('selector not bound to its certified body element')
    value=E.bulk_component(context,Y,rule,case);name=f'element_{case["element"]:03d}.json'
   else:value=E.surface_component(context,Y,rule,case);name=f'face{case["face"]}_edge{case["edge"]:02d}.json'
   path=out/part/name;write(path,value);files.append({'file':path.relative_to(ROOT).as_posix(),'sha256':sha(path)})
   line=json.dumps({'phase':part,'complete':index+1,'total':144 if part=='bulk' else 24,'seconds':time.monotonic()-started});log.write(line+'\n');log.flush();print(line,flush=True)
  assert_inputs(p);result={'status':'COMPONENT_PART_COMPLETE_UNREVIEWED','record_id':record,'part':part,'domains_sha256':sha(out/'domains.json'),'files':files,'seconds':time.monotonic()-started}
  write(out/(part+'_results.json'),result);write(out/(part+'_terminal.json'),{'status':'FINISHED','exit_code':0,'pid':os.getpid(),'utc':now(),'results_sha256':sha(out/(part+'_results.json'))})
 except BaseException as error:
  write(out/(part+'_failure.json'),{'error':str(error),'traceback':traceback.format_exc()});write(out/(part+'_terminal.json'),{'status':'FAILED','exit_code':1,'pid':os.getpid(),'utc':now()});raise
 finally:log.close()

def assemble(record):
 out=outdir(record);gate(out,'integration');p=read(out/'protocol.json');assert_inputs(p);d=read(out/'domains.json');context=C.load_context(record);precision.prec=PRECISION;core.A.cache_clear()
 g=[arb(0)for _ in range(611)];H={};packets=[]
 for part,total in [('bulk',144),('surface',24)]:
  result=read(out/(part+'_results.json'));terminal=read(out/(part+'_terminal.json'))
  if terminal['exit_code']!=0 or terminal['results_sha256']!=sha(out/(part+'_results.json'))or len(result['files'])!=total or result['domains_sha256']!=sha(out/'domains.json'):raise ValueError('complete current terminal parts required')
  for item in result['files']:
   if sha(ROOT/item['file'])!=item['sha256']:raise ValueError('changed local packet')
   v=read(ROOT/item['file']);ids=v['global_indices'];n=len(ids)
   if v['record_id']!=record or len(v['force_integral_minus_cache'])!=n or len(v['H_integral_minus_cache_upper'])!=n*(n+1)//2:raise ValueError('same record/complete local matrix required')
   ball=lambda b:core.A(F(b['lower_exact'])).union(core.A(F(b['upper_exact'])))
   for i,x in zip(ids,v['force_integral_minus_cache']):g[i]+=ball(x)
   for x in v['H_integral_minus_cache_upper']:
    key=tuple(sorted((ids[x['i']],ids[x['j']])));H[key]=H.get(key,arb(0))+ball(x['enclosure'])
   packets.append(item)
 R=arb_mat(611,611,[core.A(float(v))for v in context.R.ravel()]);RF=R*arb_mat(611,1,g);rows=[arb(0)for _ in g]
 for(i,j),v in H.items():
  rows[i]+=core.upper(v)
  if i!=j:rows[j]+=core.upper(v)
 AR=arb_mat(611,611,[core.upper(v)for v in R.entries()]);BH=AR*arb_mat(611,1,rows)
 RFbound=F(str(max(core.upper(RF[i,0])for i in range(611)).fmpq()));b=F(str(max(core.upper(BH[i,0])for i in range(611)).fmpq()))
 z=F(p['cached_Z_upper_exact'])+b;eta=F(p['cached_eta_upper_exact'])+RFbound;r=F(p['normal_radius_exact']);ok=z<1 and eta+z*r<r
 result={'status':'CANDIDATE_UNIFORM_EXACT_FE_ROOT_UNREVIEWED'if ok else'UNRESOLVED_SUFFICIENT_ROOT_TEST','record_id':record,'protocol_sha256':sha(out/'protocol.json'),'domains_sha256':sha(out/'domains.json'),'component_packets':packets,
 'force_integral_minus_cache':[core.pack(v)for v in g],'H_integral_minus_cache_upper':[{'i':i,'j':j,'enclosure':core.pack(v)}for(i,j),v in sorted(H.items())],
 'preconditioned_force_defect':[core.pack(RF[i,0])for i in range(611)],'H_absolute_row_bounds':[core.pack(v)for v in rows],
 'RF_defect_uniform_upper_exact':str(RFbound),'b_upper_exact':str(b),'z_upper_exact':str(z),'eta_upper_exact':str(eta),'selfmap_lhs_exact':str(eta+z*r),'r_exact':str(r),'sufficient_conditions_pass':bool(ok),
 'root_error_upper_exact':str(eta/(1-z))if ok else None,'root_assertion_independently_reviewed':False,
 'scope':'One original predictor cube and unchanged C for this source dose only. No integration root assertion before independent domain/kernel/reduction/composition review; no full physical identified set, PDE, larger prediction or biological validation.'}
 write(out/'assembled.json',result);print(json.dumps({'phase':'assembled','status':result['status'],'RF':float(RFbound),'b':float(b),'z':float(z),'lhs':float(eta+z*r),'r':float(r),'sha256':sha(out/'assembled.json')}),flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('command',choices=['prepare','domains','bulk','surface','assemble']);a.add_argument('record');v=a.parse_args()
 if v.command=='prepare':prepare(v.record)
 elif v.command=='domains':domains(v.record)
 elif v.command in('bulk','surface'):integrals(v.record,v.command)
 else:assemble(v.record)
