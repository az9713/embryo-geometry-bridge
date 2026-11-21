from pathlib import Path
import hashlib,json
import av,cv2,numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
file=ROOT/'data/raw/41467_2025_62483_MOESM5_ESM.avi';cv=cv2.VideoCapture(str(file),cv2.CAP_FFMPEG)
count=0;pixels=0;changed_pixels=0;channel_changed=np.zeros(3,dtype=np.int64);sums=np.zeros(3,dtype=np.int64);maxima=np.zeros(3,dtype=np.int64);selected={}
with av.open(str(file)) as container:
    for i,frame in enumerate(container.decode(container.streams.video[0])):
        ok,bgr=cv.read()
        if not ok:raise ValueError('comparison binding ended early')
        a=frame.to_ndarray(format='rgb24');b=cv2.cvtColor(bgr,cv2.COLOR_BGR2RGB)
        if a.shape!=b.shape:raise ValueError('pixel dimensions differ')
        diff=np.abs(a.astype(np.int16)-b.astype(np.int16))
        maxima=np.maximum(maxima,diff.max(axis=(0,1)));sums+=diff.sum(axis=(0,1),dtype=np.int64);channel_changed+=(diff!=0).sum(axis=(0,1),dtype=np.int64)
        changed_pixels+=int(np.any(diff!=0,axis=2).sum());pixels+=a.shape[0]*a.shape[1];count+=1
        if i in (0,79,158):selected[str(i)]={'max_abs_difference_per_RGB_channel':diff.max(axis=(0,1)).tolist(),
            'mean_abs_difference_per_RGB_channel':diff.mean(axis=(0,1)).tolist(),'fraction_pixels_any_channel_changed':float(np.any(diff!=0,axis=2).mean())}
    extra,_=cv.read()
cv.release()
result={'status':'DECODER_DIFFERENCE_QUANTIFIED_NOT_FORCED_IDENTICAL','source_SHA256':hashlib.sha256(file.read_bytes()).hexdigest(),
'frame_count':count,'extra_OpenCV_frame':bool(extra),'compared_pixels':pixels,'max_abs_difference_per_RGB_channel':maxima.tolist(),
'mean_abs_difference_per_RGB_channel':(sums/pixels).tolist(),'fraction_changed_per_RGB_channel':(channel_changed/pixels).tolist(),
'fraction_pixels_any_RGB_channel_changed':changed_pixels/pixels,'selected_frames':selected,
'PyAV_version':av.__version__,'PyAV_library_versions':{k:list(v) for k,v in av.library_versions.items()},
'OpenCV_version':cv2.__version__,'OpenCV_FFmpeg_build_lines':[s.strip() for s in cv2.getBuildInformation().splitlines() if any(x in s for x in ['FFMPEG:','avcodec:','avformat:','avutil:','swscale:'])],
'scope':'Numerical decoder/color-conversion difference only. Neither RGB agreement nor disagreement is an anatomical calibration, noise bound or biological finding.'}
target=OUT/'movie_1/decoder_difference.json';assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['frame_count','max_abs_difference_per_RGB_channel','mean_abs_difference_per_RGB_channel','fraction_pixels_any_RGB_channel_changed']}))
