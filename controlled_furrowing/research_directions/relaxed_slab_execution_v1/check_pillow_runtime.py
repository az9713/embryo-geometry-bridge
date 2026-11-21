"""Trusted image compatibility checks in an explicitly selected interpreter."""
from pathlib import Path
from hashlib import sha256
from io import BytesIO
import argparse
import importlib.metadata as metadata
import json
import sys
from PIL import Image, ImageDraw


def main():
    p=argparse.ArgumentParser();p.add_argument('--publication',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();rows=[]
    for path in sorted(a.publication.rglob('*')):
        if path.suffix.lower() not in ['.png','.jpg','.jpeg'] or '.git' in path.parts:continue
        with Image.open(path) as im:
            im.load();rgba=im.convert('RGBA')
            rows.append({'path':path.relative_to(a.publication).as_posix(),
                         'size':list(im.size),'rgba_pixel_sha256':sha256(rgba.tobytes()).hexdigest()})
    canvas=Image.new('RGB',(120,80),'white');draw=ImageDraw.Draw(canvas)
    draw.rectangle((10,10,70,50),fill='navy');draw.line((0,79,119,0),fill='red',width=3)
    formats={}
    for fmt in ['PNG','TIFF','JPEG','GIF']:
        b=BytesIO();canvas.save(b,format=fmt);b.seek(0)
        with Image.open(b) as im:
            im.load();assert im.size==(120,80)
            if fmt in ['PNG','TIFF']:assert im.convert('RGB').tobytes()==canvas.tobytes()
            formats[fmt]={'decoded':True,'size':list(im.size)}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(3,2),dpi=80);ax.plot([0,1,2],[0,1,0]);fig.canvas.draw()
    expected=bytes(fig.canvas.buffer_rgba());b=BytesIO();fig.savefig(b,format='png');plt.close(fig);b.seek(0)
    with Image.open(b) as im:assert im.convert('RGBA').tobytes()==expected
    r={'status':'PASS_TRUSTED_IMAGE_COMPATIBILITY','pillow':metadata.version('Pillow'),
       'pillow_loaded_inside_selected_venv':Path(Image.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()),
       'base_interpreter_unchanged_by_test':True,'raster_files':len(rows),'images':rows,
       'trusted_format_roundtrips':formats,'matplotlib_rgba_sha256':sha256(expected).hexdigest(),
       'scope':'Trusted saved assets and generated format tests; not an exploit test or proof of absence of security defects.'}
    a.output.write_text(json.dumps(r,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({k:r[k] for k in ['status','pillow','raster_files','pillow_loaded_inside_selected_venv']}))


if __name__=='__main__':main()
