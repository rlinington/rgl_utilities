#!/usr/bin/env python3

"""Tools for inserting deposited data from NP Atlas site to NP Atlas curator"""


import json
from pathlib import Path
import requests
import time
import ast


def submission_to_curator(submission_dir):
    """Convert submitted jsons into format required for insertion to NP Atlas curator"""

    insertion_data = {"instructions": "User submissions to 09/2022",
                      "articles": []
                      }

    for submission in submission_dir.glob("depositions_*.json"):
        with open(submission, "rb") as f:
            input_data = json.load(f)

        for deposition in input_data:
            payload = {'doi': deposition["doi"]}
            r = requests.post('https://litapi.liningtonlab.org/article/', params=payload)
            if r.status_code != 200:
                payload = {'pii': deposition["doi"]}
                r = requests.post('https://litapi.liningtonlab.org/article/', params=payload)
                if r.status_code != 200:
                    print("Error for " + str(deposition["id"]))
                    continue
            if r.status_code == 200:
                print("Processing " + deposition["doi"])
            else:
                print("ERROR for " + deposition["doi"])

            citation = r.json()
            curator_article = {"pmid": citation["pmid"],
                               "doi": citation["doi"],
                               "journal": citation["journal"],
                               "year": citation["year"],
                               "volume": citation["volume"],
                               "issue": citation["issue"],
                               "pages": citation["pages"],
                               "title": citation["title"],
                               "authors": citation["authors"],
                               "abstract": citation["abstract"],
                               "notes": "",
                               "compounds": []
                               }
            for compound in json.loads(deposition["compounds"]):
                curator_article["compounds"].append({"name": compound["name"],
                                                     "smiles": compound["smiles"],
                                                     "source_organism": compound["genus"] + " " + compound["species"]})
            insertion_data["articles"].append(curator_article)
            time.sleep(1)

    export_file = submission_dir / "full_depositions_curator_input.json"
    with open(export_file, "w") as g:
        json.dump(insertion_data, g)


if __name__ == "__main__":

    depositions_dir = Path("/Users/roger/Documents/Chemistry/NP_Atlas/Depositions")
    submission_to_curator(depositions_dir)