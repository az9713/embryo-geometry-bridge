"""Official supplemental movies: immutable acquisition and full-frame decoding.

Source images/contact sheets are local inspection assets, not public adaptations.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import datetime, hashlib, json, re, subprocess, urllib.request, xml.etree.ElementTree as ET
import fitz
from PIL import Image,ImageDraw

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
RAW=ROOT/'data/raw'
def digest(p,algo='sha256'):
    h=hashlib.new(algo)
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
def run(args):
    r=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',errors='replace')
    if r.returncode:raise RuntimeError({'args':args,'exit':r.returncode,'stderr':r.stderr})
    return r
meta=json.loads((RAW/'PMC12361397.1.json').read_text(encoding='utf-8'))
listing=ET.parse(RAW/'oa_object_listing.xml').getroot()
objects={x.find('{*}Key').text.split('/')[-1]:{'size':int(x.find('{*}Size').text),'etag':x.find('{*}ETag').text}for x in listing.findall('{*}Contents')}
md5s={x.split('/')[-1].split('?')[0]:x.split('?md5=')[-1] for x in meta['media_urls']if '?md5=' in x}
script_sha=digest(Path(__file__))
protocol={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':script_sha,'official_base':'https://pmc-oa-opendata.s3.amazonaws.com/PMC12361397.1/','filenames':[f'41467_2025_62483_MOESM{i}_ESM.avi'for i in range(8,15)],'metadata_sha256':digest(RAW/'PMC12361397.1.json'),'object_listing_sha256':digest(RAW/'oa_object_listing.xml'),'primary_xml_sha256':digest(RAW/'primary.xml'),'description_pdf_sha256':digest(RAW/'41467_2025_62483_MOESM3_ESM.pdf'),'policy':'Download exact official bytes; require size and official MD5. ffprobe all streams/tags/count_frames; ffmpeg -xerror map0 full decode to framemd5. Nine-frame local-only contacts include first/middle/last. No contours, source fit, model or publication changes.'}
if (OUT/'protocol.json').exists():raise RuntimeError('Do not overwrite protocol')
save(OUT/'protocol.json',protocol)
pdf=fitz.open(RAW/'41467_2025_62483_MOESM3_ESM.pdf')
(OUT/'description_pdf_text.txt').write_text('\n\n'.join(f'PAGE {i+1}\n'+p.get_text()for i,p in enumerate(pdf)),encoding='utf-8')
def audit(filename):
    index=int(re.search(r'MOESM(\d+)',filename).group(1));movie=index-4
    folder=OUT/f'movie_{movie:02d}';folder.mkdir(exist_ok=True)
    p=RAW/filename;url=protocol['official_base']+filename;http={}
    if not p.exists():
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 source-contract-audit'})
        with urllib.request.urlopen(req,timeout=120)as f:
            b=f.read();http={'status':f.status,'final_url':f.url,'headers':dict(f.headers)}
        if len(b)!=objects[filename]['size']:raise RuntimeError('download size mismatch '+filename)
        if hashlib.md5(b).hexdigest()!=md5s[filename]:raise RuntimeError('official MD5 mismatch '+filename)
        p.write_bytes(b)
    assert p.stat().st_size==objects[filename]['size']
    assert digest(p,'md5')==md5s[filename]
    probe=run(['ffprobe','-v','error','-count_frames','-show_streams','-show_format','-of','json',str(p)])
    (folder/'ffprobe.json').write_text(probe.stdout,encoding='utf-8');data=json.loads(probe.stdout)
    result=run(['ffmpeg','-v','error','-xerror','-i',str(p),'-map','0','-f','framemd5','-'])
    (folder/'decoded_framemd5.txt').write_text(result.stdout,encoding='utf-8');(folder/'decode_stderr.txt').write_text(result.stderr,encoding='utf-8')
    lines=[x for x in result.stdout.splitlines()if x and not x.startswith('#')];counts={}
    for line in lines:idx=line.split(',')[0].strip();counts[idx]=counts.get(idx,0)+1
    videos=[x for x in data['streams']if x['codec_type']=='video'];assert len(videos)==1
    v=videos[0];n=int(v['nb_read_frames']);assert counts[str(v['index'])]==n
    indices=sorted(set(round(i*(n-1)/8)for i in range(9)))
    for j in indices:
        target=folder/f'frame_{j:04d}.png'
        r=run(['ffmpeg','-v','error','-xerror','-i',str(p),'-vf',f'select=eq(n\\,{j})','-frames:v','1','-vsync','0',str(target)])
    images=[Image.open(folder/f'frame_{j:04d}.png').convert('RGB')for j in indices]
    w,h=images[0].size;thumbw=min(w,900);thumbh=round(h*thumbw/w);sheet=Image.new('RGB',(3*thumbw,3*(thumbh+30)),(15,19,25));draw=ImageDraw.Draw(sheet)
    for k,(j,img)in enumerate(zip(indices,images)):
        x=(k%3)*thumbw;y=(k//3)*(thumbh+30);sheet.paste(img.resize((thumbw,thumbh)),(x,y+30));draw.text((x+8,y+8),f'Movie {movie}: decoded frame {j}/{n-1}',fill=(230,230,230))
    sheet.save(folder/'contact_sheet.png')
    rec={'movie':movie,'file':str(p.relative_to(ROOT)).replace('\\','/'),'url':url,'http':http,'expected_size':objects[filename]['size'],'actual_size':p.stat().st_size,'official_expected_md5':md5s[filename],'actual_md5':digest(p,'md5'),'actual_sha256':digest(p),'s3_etag_not_assumed_md5':objects[filename]['etag'],'ffprobe_streams':data['streams'],'ffprobe_format':data['format'],'full_decode_exit_zero':True,'decoded_frame_counts_by_stream':counts,'frames_extracted':indices,'first_middle_last_in_contact':all(x in indices for x in [0,round((n-1)/2),n-1]),'assets':[{'file':str(q.relative_to(ROOT)).replace('\\','/'),'sha256':digest(q)}for q in sorted(folder.iterdir())if q.is_file()]}
    save(folder/'receipt.json',rec);print(json.dumps({'movie':movie,'frames':n,'sha256':rec['actual_sha256']},ensure_ascii=True),flush=True)
    return rec
with ThreadPoolExecutor(max_workers=3)as ex:rows=list(ex.map(audit,protocol['filenames']))
save(OUT/'acquisition_decode_results.json',{'status':'PASS_ACQUISITION_AND_FULL_DECODE','protocol_sha256':digest(OUT/'protocol.json'),'records':rows,'visual_inspection_pending':True,'scope':'Download, container metadata and complete decoding only; biological interpretation/visual inspection separate.'})
