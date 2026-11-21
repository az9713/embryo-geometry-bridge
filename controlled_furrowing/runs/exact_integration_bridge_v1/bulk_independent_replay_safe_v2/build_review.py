"""Version independent local point implementations without changing old files."""
from pathlib import Path
import ast,hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    old=ROOT/'runs/exact_integration_bridge_v1/bulk_pilot_independent_replay/replay.py';s=old.read_text();s=s[:s.index('def prepare():')]
    s=s.replace('Independent four-element POINT replay. Never an assembled/uniform bridge.','SAFE V2 independent local point library; no automatic execution or old-output premises.')
    s=s.replace('bridge_design_probe','safe_kernel_v2').replace('import exact_fe_target as T','import exact_fe_target_safe as T').replace('import polynomial_jet as DUAL','import polynomial_jet_safe as DUAL')
    lines=[x for x in s.splitlines()if not x.startswith(('ROOTOUT=','EXPECTED_DRIVER=','EXPECTED_PROTOCOL=','EXPECTED_RESULTS='))]
    dest=HERE/'kernel.py';assert not dest.exists();dest.write_text('\n'.join(lines)+'\n',encoding='utf8')
    src=ROOT/'runs/exact_integration_bridge_v1/bulk_assembly_B_independent_review/audit.py';text=src.read_text();node=next(x for x in ast.parse(text).body if isinstance(x,ast.FunctionDef)and x.name=='sample_element');sample=ast.get_source_segment(text,node)
    needle="a,b=I.series(I.T.sub(jet['value'],I.T.ONE),18);vol=I.DUAL.contract_volumetric_defect(jet,a,b,fun.defect);drive="
    replace="a,b=I.series(I.T.sub(jet['value'],I.T.ONE),18);vol=I.DUAL.contract_volumetric_defect(jet,a,b,fun.defect);direct_g=fun.apply(I.T.mul(a,jet['gradient'][0]));direct_h=fun.apply(I.T.add(I.T.mul(a,jet['hessian_upper'][0,0]),I.T.mul(b,I.T.mul(jet['gradient'][0],jet['gradient'][0]))));direct_checks={'gradient0_direct_moment_product_overlaps_dual':direct_g.overlaps(vol['gradient'][0]),'H00_direct_moment_product_overlaps_dual':direct_h.overlaps(vol['hessian_upper'][0,0])};assert all(direct_checks.values());drive="
    assert needle in sample;sample=sample.replace(needle,replace)
    sample=sample.replace("'element':index,'order':18", "'element':index,'checks':direct_checks,'order':18")
    f=HERE/'sample.py';assert not f.exists();f.write_text('from fractions import Fraction as F\nfrom flint import arb\nimport time\n'+sample+'\n',encoding='utf8')
    p=HERE/'code_lineage.json';assert not p.exists();p.write_text(json.dumps({'old_helper_sha256':sha(old),'old_sample_driver_sha256':sha(src),'safe_kernel_sha256':sha(dest),'safe_sample_sha256':sha(f),'changes':['Safe target/helper imports; old execution contexts removed','Direct gradient0/H00 monomial contractions added to each new sample','No historical point output is a validated premise']},indent=2)+'\n')
    print(json.dumps({'kernel_sha256':sha(dest),'sample_sha256':sha(f)}))
if __name__=='__main__':main()
