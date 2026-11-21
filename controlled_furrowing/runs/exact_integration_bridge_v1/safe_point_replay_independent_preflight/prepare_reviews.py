"""Version own previous preparation-only checks for safe point libraries."""
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    plans=[('bulk','bulk_assembly_B_independent_preflight','bulk_assembly_B','bulk_assembly_B_safe_v2','exact_bulk_bridge_assembly.py','exact_bulk_bridge_assembly_safe.py','41fc1735566f72ec42e249aa00fbc61f9e26b9ef55ea96f3e790bb82adac83ff','84a36fe274e3606d77adfe668fcdfc975e84996319d7bfc3cf4883d31de94950','8b14c671ca877d96256225014dcc42fef97df6ea6f7cabc3af860b857398cb26','92488086a01c04d8e0e307f606484a05a701328d192465b7ea9757a97c64d015'),('surface','surface_assembly_B_independent_preflight','surface_assembly_B','surface_assembly_B_safe_v2','exact_surface_bridge_assembly.py','exact_surface_bridge_assembly_safe.py','3d6c73ed55b710fe3ffa5b34a6834e0c3d3601c737c502348aaf75c248ae95ab','1012f6f679d48661890310daa0a4d05064a380930a55843e22a88995b0680fb3','2bba6d5e11ead452ce195a035effb612ff726d38039998cfa957a368b53068c3','55e0df292e44081d58863dc8c17b1b568efd6e1ece7177f38da203de6f0d24e7')]
    for kind,folder,oldout,newout,oldsource,newsource,oldhash,newhash,oldprot,newprot in plans:
        original=ROOT/'runs/exact_integration_bridge_v1'/folder/'review.py';src=original.read_text()
        src=src.replace('C=H.parents[2]',"C=next(p for p in H.parents if (p/'scripts/exact_fe_target_safe.py').exists())")
        src=src.replace('/'+oldout+"'",'/'+newout+"'").replace(oldsource,newsource).replace(oldhash,newhash).replace(oldprot,newprot)
        src=src.replace("pre['status']=='PASS_PREPARATION_ONLY'","pre['status']=='PASS_SAFE_V2_PREPARATION_ONLY'")
        src=src.replace("'ONE saved center' in p['scope']","'original saved B center' in p['scope']")
        src=src.replace("'ONE saved B-center point' in p['scope']","'original saved B center' in p['scope']")
        src=src.replace('import exact_surface_bridge_pilot as P','import exact_surface_bridge_safe as P')
        if kind=='bulk':
            src=src.replace('ctx.prec=512;I.C.A.cache_clear();',"import exact_fe_target_safe as SAFE_T\n    I.T=SAFE_T  # Every independently formed jet now uses exact-zero-only arithmetic.\n    ctx.prec=512;I.C.A.cache_clear();")
        else:src=src.replace('1787b0f3409aa40b7e6525269775c6c78bbab55de57acbd33df258abe15f22f3','cbdca7015f7bdbbd5be092509ea8824b68190279d5ee967b22d17e6f778fa07f')
        # No launch approval can be written unless common semantic checks pass.
        src=src.replace('def main():\n',"def main():\n    common=read(H.parent/'common_review.json')\n    assert common['status']=='PASS' and all(v is True for v in common['checks'].values())\n",1)
        src=src.replace("'checks':checks,'source_sha256'","'checks':{**checks,**common['checks']},'common_review_sha256':sha(H.parent/'common_review.json'),'source_sha256'")
        src=src.replace("(H/'review.json').write_text","assert not (H/'review.json').exists() and not (OUT/'prelaunch_independent_review.json').exists()\n    (H/'review.json').write_text")
        dest=HERE/kind/'review.py';dest.parent.mkdir(parents=True,exist_ok=True);assert not dest.exists();dest.write_text(src,encoding='utf8')
    print(json.dumps({'bulk_review_sha256':sha(HERE/'bulk/review.py'),'surface_review_sha256':sha(HERE/'surface/review.py')}))
if __name__=='__main__':main()
