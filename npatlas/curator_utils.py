#!/usr/bin/env python3

"""Tools for inserting deposited data from NP Atlas site to NP Atlas curator"""


import json
from pathlib import Path
import requests
import time
import json


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


def find_bad_character(name_string):
    """Evaluate string for non-utf-8 characters"""
    cleaned = ''.join(c for c in name_string if c.isprintable())
    if len(name_string) != len(cleaned):
        print(name_string)
        print(cleaned)


if __name__ == "__main__":

    # Create depositions from jsons
    # depositions_dir = Path("/Users/roger/Documents/Chemistry/NP_Atlas/Depositions")
    # submission_to_curator(depositions_dir)

    # Find bad names strings

    with open("/Users/roger/Git/rgl_utilities/data/NP_Atlas/response_1666116851678.json", "rb") as f:
        name_list = json.load(f)
    for name_string in name_list:
        find_bad_character(name_string)

    # test_strings = ["(18R)-18-O--d-glucopyranosyl-(1-3)-B[-d-glucopyranosyl-(1-2)]-B-d-glucopyranoside of allo-murolic acid",
    #                 "(18R)-18-O-B-d-glucopyranosyl-(1-2)-B-d-glucopyranoside of allo-murolic acid",
    #                 "(18R)-O-B\f-d-glucopyranosyl-(1-3)-B\f-d-glucopyranoside of murolic acid",
    #                 "����",
    #                 "🿾🿿𯿾𯿿𿿾𿿿񏿾񏿿񟿾񟿿񯿾񯿿񿿾񿿿򏿾򏿿"]
    # for test_string in test_strings:
    #     find_bad_character(test_string)
