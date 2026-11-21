"""Read-only transitive local-import/hash coverage for the prepared direct job.

Static import discovery is conservative across all AST branches. It does not
claim to discover dynamically constructed imports or pin external libraries.
"""
from pathlib import Path
import ast,hashlib,json

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=ROOT/'scripts'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def closure(entry):
    todo=[entry];graph={}
    while todo:
        name=todo.pop()
        if name in graph:continue
        imports=[]
        for node in ast.walk(ast.parse((SCRIPTS/name).read_text(encoding='utf-8'))):
            modules=([node.module] if isinstance(node,ast.ImportFrom) else
                     [alias.name for alias in node.names] if isinstance(node,ast.Import) else [])
            for module in modules:
                if module and (SCRIPTS/(module+'.py')).is_file():imports.append(module+'.py')
        graph[name]=sorted(set(imports));todo+=sorted(set(imports)-graph.keys())
    return graph


def main():
    path=ROOT/'runs/material_cohort_direct_v1/protocol.json'
    p=json.loads(path.read_text());graph=closure('material_cohort_direct.py');locked=p['source_hashes']
    missing=sorted(set(graph)-locked.keys())
    changed=[name for name,value in locked.items() if digest(SCRIPTS/name)!=value]
    checks={'all_transitive_static_local_imports_guarded':not missing,
            'all_declared_source_hashes_current':not changed,
            'entry_point_hash_bound':locked.get('material_cohort_direct.py')==digest(SCRIPTS/'material_cohort_direct.py')}
    receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
             'protocol_sha256':digest(path),'verifier_sha256':digest(__file__),
             'local_import_graph':graph,'uncovered_local_modules':missing,'changed_locked_sources':changed,
             'all_locked_source_hashes':locked,
             'scope':'Static local-import dependency coverage and actual source hashes before direct replay; no numerical, dynamic-import, environment or five-gate completion claim.'}
    out=ROOT/'evidence/material_direct_dependency_closure.json'
    out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'local_modules':len(graph),'locked_modules':len(locked),'missing':missing,'changed':changed}))
    if receipt['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
