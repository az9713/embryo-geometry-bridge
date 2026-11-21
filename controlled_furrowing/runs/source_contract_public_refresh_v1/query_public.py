"""Read public repository search metadata only. No dataset file fetches."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import requests

HERE=Path(__file__).resolve().parent
TITLE="Endogenous OptoRhoGEFs reveal biophysical principles of epithelial tissue furrowing"
DOI="10.1038/s41467-025-62483-6"
GITHUB="https://api.github.com/search/repositories"
ZEN="https://zenodo.org/api/records"
FIG="https://api.figshare.com/v2/articles/search"
EMBL="https://www.ebi.ac.uk/biostudies/api/v1/search"
EMPIAR="https://www.ebi.ac.uk/ebisearch/ws/rest/empiar"
DATACITE="https://api.datacite.org/dois"
TASKS={
    "zenodo_doi":("GET",ZEN,{"q":'related.identifiers.identifier:"'+DOI+'"',"size":12}),
    "zenodo_title":("GET",ZEN,{"q":'"OptoRhoGEFs" OR "'+TITLE+'"',"size":12}),
    "figshare_title":("POST",FIG,{"search_for":":title: Endogenous OptoRhoGEFs","page_size":12}),
    "figshare_related_doi":("POST",FIG,{"resource_doi":DOI,"page_size":12}),
    "biostudies_title":("GET",EMBL,{"query":"OptoRhoGEFs","pageSize":12}),
    "biostudies_doi":("GET",EMBL,{"query":DOI,"pageSize":12}),
    "biostudies_pmc_candidate":("GET","https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC12361397",{}),
    "empiar_title":("GET",EMPIAR,{"query":"OptoRhoGEFs","format":"json","size":12}),
    "github_project":("GET",GITHUB,{"q":"OptoRhoGEF OR OptoRhoGEFs in:name,description,readme","per_page":12}),
    "github_doi":("GET",GITHUB,{"q":'"'+DOI+'" in:readme',"per_page":12}),
    "github_lab":("GET",GITHUB,{"q":"KaszaLab in:name,description","per_page":12}),
    "datacite_related_doi":("GET",DATACITE,{"query":'relatedIdentifiers.relatedIdentifier:"'+DOI+'"',"page[size]":12}),
    "datacite_title":("GET",DATACITE,{"query":'titles.title:"OptoRhoGEFs"',"page[size]":12}),
    "europepmc_datalinks":("GET","https://www.ebi.ac.uk/europepmc/webservices/rest/MED/40825942/datalinks",{"format":"json"}),
    "europepmc_full_text":("GET","https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12361397/fullTextXML",{}),
    "europepmc_article_metadata":("GET","https://www.ebi.ac.uk/europepmc/webservices/rest/search",{"query":"EXT_ID:40825942","format":"json"}),
}

def request(item):
    key,(method,url,params)=item
    row={"name":key,"method":method,"public_url":url,"metadata_query":params}
    try:
        args={"json" if method=="POST" else "params":params}
        r=requests.request(method,url,timeout=(7,18),headers={"User-Agent":"public-repository-readback/1.0"},**args)
        row.update(status=r.status_code,final_url=r.url,media_type=r.headers.get("Content-Type"),bytes=len(r.content),
                   response_sha256=sha256(r.content).hexdigest())
        if len(r.content)>2_000_000:
            row["retained_bytes"]=False
            return row
        suffix="xml" if "xml" in (row["media_type"] or "").lower() else "json" if "json" in (row["media_type"] or "").lower() else "txt"
        fp=HERE/"metadata"/(key+"."+suffix)
        with fp.open("xb") as f:
            f.write(r.content)
        row["retained_path"]=fp.relative_to(HERE).as_posix()
    except requests.RequestException as exc:
        row.update(status=None,transport_failure=f"{type(exc).__name__}: {exc}")
    return row


if __name__=="__main__":
    (HERE/"metadata").mkdir(exist_ok=False)
    with ThreadPoolExecutor(max_workers=8) as workers:
        result=list(workers.map(request,TASKS.items()))
    with (HERE/"metadata_query_receipt.json").open("x",encoding="utf-8",newline="\n") as f:
        json.dump({"utc":datetime.now(timezone.utc).isoformat(),"doi":DOI,"requests":result,
                   "scope":"Public metadata query and article text only; no artifacts, data files or contacting authors."},f,indent=2)
        f.write("\n")
    for row in result:
        print(row["name"],row["status"],row.get("media_type",""))

