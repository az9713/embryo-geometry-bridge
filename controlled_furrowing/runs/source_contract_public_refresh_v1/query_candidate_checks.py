"""Resolve an escaped index query and inspect only small public README texts."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import json
from query_public import HERE, DOI, EMPIAR, request

QUERIES={
    "empiar_quoted_doi":("GET",EMPIAR,{"query":f'"{DOI}"',"format":"json","size":12}),
    "lab_fork_readme":("GET","https://raw.githubusercontent.com/kaszalab/image_translation/main/README.md",{}),
    "lab_fork_setup":("GET","https://raw.githubusercontent.com/kaszalab/image_translation/main/kaszalab_setup.md",{}),
    "lab_fork_download_manifest":("GET","https://raw.githubusercontent.com/kaszalab/image_translation/main/download_data.sh",{}),
    "lab_xenopus_readme":("GET","https://raw.githubusercontent.com/kaszalab/xenopus-protein-localization/main/README.md",{}),
}

if __name__=="__main__":
    with ThreadPoolExecutor(max_workers=5) as pool:
        rows=list(pool.map(request,QUERIES.items()))
    with (HERE/"metadata_candidate_checks.json").open("x",encoding="utf-8",newline="\n") as f:
        json.dump({"utc":datetime.now(timezone.utc).isoformat(),"requests":rows,"scope":"Public index and small text file reads only. No downloader executed."},f,indent=2)
        f.write("\n")
    for row in rows:
        print(row["name"],row["status"])

