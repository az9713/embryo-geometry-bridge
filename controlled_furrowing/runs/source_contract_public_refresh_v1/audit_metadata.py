"""Audit saved public metadata and identify the known article supplements by size/name."""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PREVIOUS=ROOT/"evidence/public_supplement_asset_manifest.json"
CONTRACT=ROOT/"docs/JOINT-OBSERVATION-CONTRACT.md"
DOI="10.1038/s41467-025-62483-6"


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def check(p,h):
    assert sha(p)==h,p


def doc(p):
    return json.loads(p.read_bytes())


def load(key):
    return doc(HERE/"metadata"/(key+".json"))


def main():
    f1="b5c36764d5cecbdE6c1fcd2aa21a0d2c3cc6ecdfd1501d73b397ff098b7bf4af".lower()
    f2="d805ff0f9c757d35fdd61283ae1591be3873d73d7daadbaac521020454654d04"
    check(PREVIOUS,f1)
    check(CONTRACT,f2)
    receipt_files=[HERE/n for n in ("metadata_query_receipt.json","metadata_followup_receipt.json","metadata_candidate_checks.json")]
    attempts=[row for p in receipt_files for row in doc(p)["requests"]]
    assert len(attempts)==27 and len({row["name"] for row in attempts})==27
    assert dict(Counter(row["status"] for row in attempts))=={200:26,400:1}
    for row in attempts:
        check(HERE/row["retained_path"],row["response_sha256"])
    listing=load("biostudies_pmc_candidate")
    hits=load("biostudies_title")
    assert listing["accno"]=="S-EPMC12361397" and hits["totalHits"]==1
    assert hits["hits"][0]["accession"]==listing["accno"] and hits["hits"][0]["files"]==16
    original={r["filename"]:r["bytes"] for r in doc(PREVIOUS)["assets"]}
    new={f["path"]:f["size"] for f in listing["section"]["files"]}
    assert new==original and len(new)==16
    assert Counter(Path(k).suffix for k in new)=={".pdf":4,".xlsx":2,".avi":10}
    root=ET.parse(HERE/"metadata/europepmc_full_text.xml").getroot()
    xlink="{http://www.w3.org/1999/xlink}href"
    pub=[x for x in root.findall("./front/article-meta/article-id") if x.get("pub-id-type")=="doi" and x.text==DOI]
    assert pub
    materials=[e for e in root.findall(".//supplementary-material") if (e.get("id") or "").startswith("MOESM")]
    assert len(materials)==16
    assert {e.find("media").get(xlink) for e in materials}==set(original)
    links=[e.attrib.get(xlink,"") for e in root.findall(".//ext-link")]
    assert not any(any(t in u.lower() for t in ("zenodo","figshare","github","datadryad","empiar","osf.io","biostudies")) for u in links)
    d=load("europepmc_datalinks")
    categories={c["Name"]:c["CategoryLinkCount"] for c in d["dataLinkList"]["Category"]}
    assert d["hitCount"]==15 and categories=={"RRID":11,"Data Citations":1,"Altmetric":1,
            "BioStudies: supplemental material and supporting data":1,"FlyBase":1}
    github=load("github_lab_repository_list")
    assert {r["full_name"] for r in github}=={"kaszalab/xenopus-protein-localization","kaszalab/image_translation"}
    zen=load("zenodo_title")["hits"]
    assert zen["total"]==1 and zen["hits"][0]["doi"]=="10.5281/zenodo.7944191"
    zero={
        "Zenodo related DOI":load("zenodo_doi")["hits"]["total"],
        "Figshare resource DOI":len(load("figshare_related_doi")),
        "Figshare distinctive token, title":len(load("figshare_distinctive_token")),
        "Figshare distinctive token, metadata":len(load("figshare_default_distinctive_token")),
        "BioStudies unqualified DOI":load("biostudies_doi")["totalHits"],
        "EMPIAR distinctive title token":load("empiar_title")["hitCount"],
        "EMPIAR quoted DOI":load("empiar_quoted_doi")["hitCount"],
        "GitHub targeted project words":load("github_project")["total_count"],
        "GitHub exact DOI in README":load("github_doi")["total_count"],
        "DataCite related DOI":load("datacite_related_doi")["meta"]["total"],
        "DataCite token title":load("datacite_title")["meta"]["total"],
    }
    assert all(n==0 for n in zero.values())
    p=HERE/"receipt.json"
    new_data={"status":"COMPLETE_BOUNDED_PUBLIC_REFRESH_GATE2_REMAINS_OPEN","doi":DOI,
        "query_receipt_sha256":{f.name:sha(f) for f in receipt_files},
        "previous_audit_sha256":{str(PREVIOUS.relative_to(ROOT)):f1,str(CONTRACT.relative_to(ROOT)):f2},
        "input_and_source_sha256":{str(f.relative_to(HERE)):sha(f) for f in (*HERE.glob("query_*.py"),Path(__file__),HERE/"PUBLIC-REFRESH.md")},
        "biostudies_supplement_accession":"S-EPMC12361397","same_name_and_size_count":16,"files_by_name_bytes":new,
        "publication_xml_supplement_names_equal_previous_manifest":True,
        "primary_dataset_files_downloaded_in_refresh":0,"europe_pmc_data_link_categories":categories,
        "bounded_successful_zero_queries":zero,"successful_queries":26,"retained_400_empiar_unquoted_retry":True,
        "matched_profile_or_tracking_or_reference_deposit_identified":False,
        "specimen_join_registration_time_and_uncertainty_established":False,
        "outside_contact_or_messages_sent":False,"held_out_large_response_values_opened":False,"gate2_complete":False,
        "scope":"Bounded public metadata inspection, with no inference of universal nonexistence. Name/size cross-check, no second supplement download."}
    with p.open("x",encoding="utf-8",newline="\n") as out:
        json.dump(new_data,out,indent=2)
        out.write("\n")
    print(json.dumps({"receipt_sha256":sha(p),"same_16":len(new),"requests":len(attempts),"http_200":26}))


if __name__=="__main__":
    main()

