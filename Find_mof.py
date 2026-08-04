import json

with open(r"F:\Qmof\qmof_database\qmof.json") as f:
    qmof = json.load(f)

for entry in qmof:
    formula = entry["info"]["formula_reduced"]

    if formula in [
        "Zn4C24H12O13",
        "C24H12O13Zn4",
        "Zn4O(C8H4O4)3"
    ]:
        print(entry["qmof_id"])
        print("Formula:", entry["info"]["formula"])
        print("Reduced:", formula)
        print("Name:", entry["name"])
        print("DOI:", entry["info"]["doi"])
        print("-"*60)