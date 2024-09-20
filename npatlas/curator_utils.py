#!/usr/bin/env python3

"""Tools for inserting deposited data from NP Atlas site to NP Atlas curator"""


import json
from pathlib import Path
import requests
import time
import json
import sys
import csv
from collections import defaultdict
from rdkit import Chem


CURATORS = {'Roger': 19, 'Ella': 832}


def submission_to_curator(submission_dir, title, scan_entries=False):
    """Convert submitted jsons into format required for insertion to NP Atlas curator
    NOTE: It is best to run scan_entries=True first to check that the entries pass basic literature API validation
    before curating all compounds. If you do this you must run it again with scan_entries=False to get the full curator
    output
    Args:
        submission_dir: Directory containing json dump of all new depositions from Atlas
        title (str): Title of dataset to be displayed in curator app
        scan_entries: Bool to allow initial Atlas query of all depositions to make sure there are no SMILES format
        errors. Returns an error in terminal if Atlas API does not return 200 status code
    Returns:
        json of entries for NP Atlas curator app"""

    insertion_data = {"instructions": title,
                      "articles": []
                      }

    # NOTE: This glob is very promiscuous. Currently, output file is also a json, so running code twice will create
    # problems because it will try to process the output file as an input.
    for submission in submission_dir.glob("*.json"):
        with open(submission, "rb") as f:
            input_data = json.load(f)

        for deposition in input_data:
            payload = {'doi': deposition["doi"]}
            retry = True
            while retry:
                r = requests.post('https://litapi.liningtonlab.org/article/', params=payload)
                if r.status_code == 200:
                    print("Processing " + deposition["doi"])
                    retry = False
                elif r.status_code != 200:
                    payload = {'pii': deposition["doi"]}
                    r = requests.post('https://litapi.liningtonlab.org/article/', params=payload)
                    if r.status_code != 200:
                        print(f"{r.status_code} Error for " + str(deposition["id"]))
                        retry_response = input("Retry article (y/n)?")
                        if retry_response == "n":
                            retry = False
                else:
                    print(f"{r.status_code} ERROR for " + deposition["doi"])
                    retry_response = input("Retry article (y/n)?")
                    if retry_response == "n":
                        retry = False
            if r.status_code == 200:
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
                    print(f"Name: {compound['name']}")
                    curate_compound = True
                    if compound["new_compound"] is False:
                        old_curate_response = input(f"Compound {compound['name']} is not new. Curate (y/n)?")
                        if old_curate_response == "n":
                            continue

                    retry_compound = True
                    while retry_compound:
                        atlas_r = requests.post(f'https://www.npatlas.org/api/v1/compounds/basicSearch?smiles='
                                                f'{compound["smiles"]}&method=full&threshold=1&origin_type=all&rank='
                                                f'all&orderby=npaid&ascending=true&skip=0&limit=10')
                        if atlas_r.status_code == 200:
                            retry_compound = False
                            if scan_entries:
                                pass
                            else:
                                if len(atlas_r.json()) > 0 and atlas_r.json()[0]["npaid"] != "NPA000001":
                                    atlas_json = atlas_r.json()[0]
                                    if atlas_json["original_name"] == compound["name"]:
                                        print(f"MATCH: for compound {compound['name']}")
                                        curate_compound = False
                                    else:
                                        include_response = input(f"Structure match but not name match to "
                                                                 f"{atlas_json['npaid']} for "
                                                                 f"{atlas_json['original_name']} vs {compound['name']}."
                                                                 f" Curate compound (y/n)?")
                                        if include_response in ["n", "N", "no", "No", "NO"]:
                                            curate_compound = False
                                        elif include_response in ["y", "Y", "yes", "Yes", "YES"]:
                                            pass
                                        else:
                                            print("Invalid response. Must be 'y' or 'n'")
                                            sys.exit()
                                else:
                                    print("No structure returned from NP Atlas")
                        elif atlas_r.status_code != 200:
                            print(f"Atlas error {atlas_r.status_code} for {deposition['id']} with compound "
                                  f"{compound['smiles']}")
                            retry_response = input("Retry compound (y/n)?")
                            if retry_response == "n":
                                retry_compound = False
                                curate_response = input("Curate compound (y/n)?")
                                if curate_response == "n":
                                    curate_compound = False
                        else:
                            print(f"Unknown atlas error {atlas_r.status_code} for {deposition['id']} "
                                  f"with compound {compound['smiles']}")
                            retry_response = input("Retry compound (y/n)?")
                            if retry_response == "n":
                                retry_compound = False
                                curate_response = input("Curate compound (y/n)?")
                                if curate_response == "n":
                                    curate_compound = False

                    if curate_compound:
                        if compound["genus"] and compound["species"]:
                            organism = compound["genus"] + " " + compound["species"]
                        elif compound["genus"] and not compound["species"]:
                            organism = compound["genus"]
                        else:
                            organism = None
                        curator_article["compounds"].append({"name": compound["name"],
                                                             "smiles": compound["smiles"],
                                                             "source_organism": organism})
                if len(curator_article["compounds"]) > 0:
                    insertion_data["articles"].append(curator_article)
            time.sleep(0.2)

    export_file = submission_dir / "full_depositions_curator_input.json"
    with open(export_file, "w") as g:
        json.dump(insertion_data, g)


def find_bad_character(name_string):
    """Evaluate string for non-utf-8 characters"""
    cleaned = ''.join(c for c in name_string if c.isprintable())
    if len(name_string) != len(cleaned):
        print(name_string)
        print(cleaned)


def submission_csv_to_json(csv_path: Path, curator: str):
    """Convert list of compounds from csv into standard deposition json format. Used in submission_to_curator to make
    new curator datasets
    Args:
        csv_path (Path): the path to the csv file containing the new compounds. Must contain only columns
        ['compound', 'doi', 'smiles', 'taxonomy', 'origin_type']
        curator (str): The name of the person who curated the data.
    Returns:
        Returns a json file containing all compounds from csv file"""
    try:
        user_id = CURATORS[curator]
    except KeyError:
        print(f"{curator} is not in the list of known curator IDs.")
        sys.exit()

    input_data = defaultdict(list)
    row_counter = 2

    with open(csv_path) as f:
        csv_f = csv.reader(f)
        headers = next(csv_f)
        if headers != ['compound', 'doi', 'smiles', 'taxonomy', 'origin_type']:
            print("WARNING: csv headers not correct. Must be ['compound', 'doi', 'smiles', 'taxonomy', 'origin_type']")
            sys.exit()
        for row in csv_f:
            # Perform QC on each row and insert if passes QC check
            if not isinstance(row[0], str) or len(row[0]) == 0:
                print(f"ERROR for compound name in row {row_counter}")
                sys.exit()
            try:
                mol = Chem.MolFromSmiles(row[2])
            except:
                print(f"ERROR for SMILES in row {row_counter}")
                sys.exit()
            if not isinstance(row[3], str) or len(row[3]) == 0:
                print(f"ERROR for taxonomy in row {row_counter}")
                sys.exit()
            if row[4] not in ['Bacterium', 'Fungi']:
                print(f"ERROR for origin_type in row {row_counter}. Options are Bacterium or Fungi")
                sys.exit()
            if row[3] not in ['Unknown-bacterium', 'Unknown-fungus']:
                split_taxonomy = row[3].split(" ", 1)
                genus = split_taxonomy[0]
                species = split_taxonomy[1]
            else:
                genus = row[3]
                species = None
            input_data[row[1]].append({'name': row[0],
                                       'smiles': row[2],
                                       'genus': genus,
                                       'species': species,
                                       'type': row[4],
                                       'new_compound': True})
            row_counter += 1

    article_data = []
    article_id = 1
    for doi, compound_list in input_data.items():
        article_data.append({'id': article_id,
                             'handled': True,
                             'doi': doi,
                             'compounds': json.dumps(compound_list),
                             'user': user_id})
        article_id += 1

    output_path = Path(csv_path.parent, csv_path.stem + "_reformatted.json")
    with open(output_path, "w") as g:
        json.dump(article_data, g)


if __name__ == "__main__":

    # Create depositions from jsons
    # depositions_dir = Path("/Users/roger/Documents/Chemistry/NP_Atlas/Depositions/20240425")
    # submission_to_curator(depositions_dir, "Depositions up to 20240531", scan_entries=False)

    # Convert csv to deposition json
    csv_path = Path("/Users/roger/Documents/Chemistry/NP_Atlas/Reassignments/20240828/Atlas_to_be_added_2_doi_forced.csv")
    submission_csv_to_json(csv_path, "Ella")
    json_dir = Path(csv_path.parent)
    submission_to_curator(json_dir, "reassignment additions, 20240828", scan_entries=False)

    # Find bad names strings

    # with open("/Users/roger/Git/rgl_utilities/data/NP_Atlas/response_1666116851678.json", "rb") as f:
    #     name_list = json.load(f)
    # for name_string in name_list:
    #     find_bad_character(name_string)

    # test_strings = ["(18R)-18-O--d-glucopyranosyl-(1-3)-B[-d-glucopyranosyl-(1-2)]-B-d-glucopyranoside of allo-murolic acid",
    #                 "(18R)-18-O-B-d-glucopyranosyl-(1-2)-B-d-glucopyranoside of allo-murolic acid",
    #                 "(18R)-O-B\f-d-glucopyranosyl-(1-3)-B\f-d-glucopyranoside of murolic acid",
    #                 "����",
    #                 "🿾🿿𯿾𯿿𿿾𿿿񏿾񏿿񟿾񟿿񯿾񯿿񿿾񿿿򏿾򏿿"]
    # for test_string in test_strings:
    #     find_bad_character(test_string)
