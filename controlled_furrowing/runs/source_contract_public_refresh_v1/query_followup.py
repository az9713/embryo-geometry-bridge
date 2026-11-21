"""Follow up repository hits with narrow public queries and the lab listing."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import json
from query_public import HERE, DOI, FIG, EMPIAR, request

MORE={
    "figshare_distinctive_token":("POST",FIG,{"search_for":":title: OptoRhoGEFs","page_size":20}),
    "figshare_default_distinctive_token":("POST",FIG,{"search_for":"OptoRhoGEFs","page_size":20}),
    "empiar_doi":("GET",EMPIAR,{"query":DOI,"format":"json","size":12}),
    "github_lab_repository_list":("GET","https://api.github.com/orgs/kaszalab/repos",{"type":"public","per_page":100}),
    "crossref_work_relation":("GET","https://api.crossref.org/works/"+DOI,{}),
    "crossref_prior_cited_work":("GET","https://api.crossref.org/works/10.15252/embj.2018100170",{}),
}

if __name__=="__main__":
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows=list(pool.map(request,MORE.items()))
    with (HERE/"metadata_followup_receipt.json").open("x",encoding="utf-8",newline="\n") as f:
        json.dump({"utc":datetime.now(timezone.utc).isoformat(),"requests":rows,"scope":"Public index reads only."},f,indent=2)
        f.write("\n")
    for row in rows:
        print(row["name"],row["status"])

