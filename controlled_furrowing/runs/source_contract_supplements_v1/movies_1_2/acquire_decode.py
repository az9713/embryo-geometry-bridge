from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import datetime,hashlib,json,os,subprocess,traceback
import requests,av,cv2
import numpy as np
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);data=json.dumps(d,indent=2,allow_nan=False)+'\n'
    if p.exists():raise ValueError('retain existing artifact: '+str(p))
    p.write_text(data,encoding='utf8')
def acquire(row,protocol):
    name=row['filename'];folder=OUT/('movie_'+str(row['movie']));folder.mkdir()
    url=protocol['url_prefix']+name;target=ROOT/'data/raw'/name;temporary=folder/(name+'.part')
    if target.exists():raise ValueError('source already exists; do not overwrite')
    started=now()
    with requests.get(url,stream=True,timeout=(20,90)) as response:
        headers=dict(response.headers);response.raise_for_status();status=response.status_code;final=response.url
        with temporary.open('xb') as f:
            for chunk in response.iter_content(1<<18):
                if chunk:f.write(chunk)
    raw=temporary.read_bytes();md5=hashlib.md5(raw).hexdigest()
    http={'url':url,'final_url':final,'status_code':status,'response_headers':headers,'requested_utc':started,'finished_utc':now(),
        'received_bytes':len(raw),'expected_official_bytes':row['expected_bytes'],'md5':md5,'expected_official_md5':row['expected_official_md5'],
        'SHA256':hashlib.sha256(raw).hexdigest()}
    write(folder/'http.json',http)
    if len(raw)!=row['expected_bytes'] or md5!=row['expected_official_md5'] or headers.get('ETag','').strip('"')!=md5:raise ValueError('officialsize/MD5/ETag mismatch')
    temporary.rename(target)
    command=['ffprobe','-v','error','-count_frames','-show_format','-show_streams','-of','json',str(target)]
    completed=subprocess.run(command,capture_output=True,text=True,encoding='utf8')
    (folder/'ffprobe.stderr.txt').write_text(completed.stderr,encoding='utf8')
    if completed.returncode:raise ValueError('ffprobe failed')
    probe=json.loads(completed.stdout);write(folder/'ffprobe.json',probe)
    frames=[]
    with av.open(str(target)) as container:
        stream=container.streams.video[0]
        for i,frame in enumerate(container.decode(stream)):
            rgb=frame.to_ndarray(format='rgb24')
            frames.append({'index':i,'pts':frame.pts,'time_base':str(frame.time_base),'media_time_seconds':frame.time,
                'shape_rgb':list(rgb.shape),'canonical_RGB24_SHA256':hashlib.sha256(np.ascontiguousarray(rgb).tobytes()).hexdigest()})
    n=len(frames)
    if not n:raise ValueError('no decoded frames')
    cv=cv2.VideoCapture(str(target),cv2.CAP_FFMPEG);cvhashes=[]
    while True:
        success,bgr=cv.read()
        if not success:break
        rgb=cv2.cvtColor(bgr,cv2.COLOR_BGR2RGB);cvhashes.append(hashlib.sha256(np.ascontiguousarray(rgb).tobytes()).hexdigest())
    cv.release()
    picks={'first':0,'middle':n//2,'last':n-1}
    contact=sorted(set(round(j*(n-1)/8) for j in range(9))|{n//2});wanted=set(picks.values())|set(contact);images={}
    with av.open(str(target)) as container:
        for i,frame in enumerate(container.decode(container.streams.video[0])):
            if i in wanted:images[i]=Image.fromarray(frame.to_ndarray(format='rgb24'))
    selected={}
    for label,index in picks.items():
        file=folder/(label+'.png');images[index].save(file)
        selected[label]={**frames[index],'PNG_file':file.relative_to(ROOT).as_posix(),'PNG_SHA256':sha(file)}
    width=720;header=26;ratio=images[0].height/images[0].width;height=round(width*ratio)
    sheet=Image.new('RGB',(width*3,(height+header)*((len(contact)+2)//3)),(15,20,29));draw=ImageDraw.Draw(sheet)
    for k,index in enumerate(contact):
        x=(k%3)*width;y=(k//3)*(height+header);draw.text((x+5,y+5),'frame '+str(index)+'; media '+str(frames[index]['media_time_seconds'])+'s',fill=(235,241,252))
        sheet.paste(images[index].resize((width,height)),(x,y+header))
    file=folder/'contact-sheet.png';sheet.save(file)
    write(folder/'all_frames.json',{'PyAV_frames':frames,'OpenCV_RGB24_hashes':cvhashes,'PyAV_version':av.__version__,'OpenCV_version':cv2.__version__})
    stream=next(s for s in probe['streams'] if s['codec_type']=='video')
    checks={'official_bytes_md5_etag_match':True,'PyAV_and_OpenCV_counts_agree':n==len(cvhashes),
        'ffprobe_actual_frame_count_agrees':int(stream.get('nb_read_frames','-1'))==n,
        'all_RGB24_frame_hashes_agree_between_bindings':[f['canonical_RGB24_SHA256'] for f in frames]==cvhashes,
        'first_mid_last_saved':len(selected)==3,'ffprobe_exit_zero':completed.returncode==0}
    result={'movie':row['movie'],'filename':name,'source_file':target.relative_to(ROOT).as_posix(),'source_SHA256':sha(target),'source_bytes':target.stat().st_size,
        'checks':checks,'decoded_frame_count':n,'selected_frames':selected,'contact_sheet_file':file.relative_to(ROOT).as_posix(),'contact_sheet_SHA256':sha(file),
        'contact_indices':contact,'ffprobe_file':(folder/'ffprobe.json').relative_to(ROOT).as_posix(),'ffprobe_SHA256':sha(folder/'ffprobe.json'),
        'HTTP_file':(folder/'http.json').relative_to(ROOT).as_posix(),'HTTP_SHA256':sha(folder/'http.json'),
        'all_frame_evidence_file':(folder/'all_frames.json').relative_to(ROOT).as_posix(),'all_frame_evidence_SHA256':sha(folder/'all_frames.json'),
        'decoder_scope':'Two bindings both useFFmpeg; this is complete-stream consistency, not two independent codec implementations.',
        'timing_scope':'Container presentation seconds are not assumed to equal the biological timestamps printed in the movie.',
        'local_preview_only_not_public_adaptation':True}
    write(folder/'decode_receipt.json',result)
    print(json.dumps({'movie':row['movie'],'frames':n,'checks':checks}),flush=True)
    return result
if __name__=='__main__':
    protocol=json.loads((OUT/'protocol.json').read_text())
    for file,digest in protocol['inputs'].items():
        if sha(ROOT/file)!=digest:raise ValueError('frozen source metadata changed')
    write(OUT/'job.json',{'status':'RUNNING','pid':os.getpid(),'started_utc':now(),'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json')})
    try:
        with ThreadPoolExecutor(max_workers=2) as executor:results=list(executor.map(lambda row:acquire(row,protocol),protocol['movies']))
        write(OUT/'acquisition_decode.json',{'status':'COMPLETE' if all(all(m['checks'].values()) for m in results) else 'CHECK_DIFFERENCES',
            'movies':results,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),
            'no_model_or_biological_annotation_inference':True,'visual_source_interpretation_pending':True})
        write(OUT/'terminal.json',{'status':'FINISHED','pid':os.getpid(),'finished_utc':now(),'results_sha256':sha(OUT/'acquisition_decode.json')})
    except BaseException:
        write(OUT/'terminal.json',{'status':'FAILED','pid':os.getpid(),'failed_utc':now(),'traceback':traceback.format_exc()});raise
