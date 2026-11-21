"""Retain relevant package versions without machine/account identifiers."""
from importlib.metadata import version
import json
from pathlib import Path
import platform

ROOT=Path(__file__).resolve().parents[1]
names=['numpy','scipy','torch','sympy','mpmath','matplotlib','openpyxl','markdown','playwright','Pillow']
versions={name:version(name) for name in names}
result={'python':platform.python_version(),'platform':platform.system(),'packages':versions,'scope':'Actual host environment used; not a fresh-environment replay or a dependency-support recommendation'}
(ROOT/'evidence/controlled_environment.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
index='--extra-index-url https://download.pytorch.org/whl/cu128\n' if '+cu128' in versions['torch'] else ''
(ROOT/'requirements.txt').write_text(index+'\n'.join(f'{name}=={value}' for name,value in versions.items())+'\n',encoding='utf-8')
print(json.dumps(result))
