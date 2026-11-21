"""Recognized text-pattern and raster-metadata privacy checks; no values emitted."""
from pathlib import Path
from hashlib import sha256
from collections import Counter
import argparse
import io
import json
import re
import tokenize

ROOT=Path(__file__).resolve().parents[1]
EXT={'.py','.md','.html','.json','.txt','.csv','.js','.css','.svg','.xml','.toml','.yml','.yaml'}
EMAIL=re.compile(rb'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}')
PATTERNS={
 'absolute_user_path':re.compile(rb'(?i)(?:[a-z]:[\\/]{1,4}(?:users|documents and settings)[\\/]{1,4}|/(?:home|Users)/[^/\s]+/)'),
 'private_workspace_url':re.compile(rb'https?://chatgpt\.com/space/'),
 'recognized_credential':re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|sk-[A-Za-z0-9_-]{24,})')}
def email(data):
    # Bound the search around the required @, avoiding quadratic scans of long
    # exact rational numbers and embedded binary archives.
    return any(EMAIL.search(data[max(0,m.start()-254):m.start()+255]) for m in re.finditer(b'@',data))
def main():
    cli=argparse.ArgumentParser();cli.add_argument('--output',type=Path,required=True);args=cli.parse_args()
    bad=[];count=Counter()
    for p in ROOT.rglob('*'):
        if not p.is_file() or '.git' in p.relative_to(ROOT).parts or p.suffix.lower() not in EXT:continue
        data=p.read_bytes();count['text_files']+=1;count['text_bytes']+=len(data)
        for label,pattern in PATTERNS.items():
            if pattern.search(data):bad.append({'path':p.relative_to(ROOT).as_posix(),'category':label})
        if email(data):
            match=True
            if p.suffix.lower()=='.py':
                match=any(email(t.string.encode()) for t in tokenize.tokenize(io.BytesIO(data).readline)
                          if t.type in (tokenize.STRING,tokenize.COMMENT))
            if match:bad.append({'path':p.relative_to(ROOT).as_posix(),'category':'email_like_text'})
    from PIL import Image
    for p in ROOT.rglob('*'):
        if p.suffix.lower() not in ('.png','.jpg','.jpeg') or '.git' in p.relative_to(ROOT).parts:continue
        with Image.open(p) as im:
            count['raster_files']+=1
            meta={**im.info,**{str(k):str(v) for k,v in im.getexif().items()}}
            for k,v in meta.items():
                if k.lower() in ('author','artist','gps','gpsinfo','ownername','cameraownername'):
                    bad.append({'path':p.relative_to(ROOT).as_posix(),'category':'identifying_raster_metadata'})
                if isinstance(v,str):
                    raw=v.encode()
                    if email(raw) or any(x.search(raw) for x in PATTERNS.values()):
                        bad.append({'path':p.relative_to(ROOT).as_posix(),'category':'recognized_raster_metadata_pattern'})
    record={'status':'PASS_RECOGNIZED_PRIVACY_SCREEN' if not bad else 'FAIL_PRIVACY_SCREEN',
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'counts':dict(count),
            'flagged_paths_and_categories':bad,'matching_values_printed':False,
            'scope':'Recognized user paths, private workspace URLs, emails and credential patterns; inspected common raster metadata. Not a universal secrecy guarantee.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':record['status'],'counts':record['counts'],'flagged_files':len(bad)}))
    if bad:raise SystemExit(1)
if __name__=='__main__':main()
