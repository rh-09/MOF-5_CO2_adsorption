import json
import os
from pymatgen.core import Structure

# ------Settings------#
struct_json_path = "qmof_structure_data.json"  # path to structure json
cif_folder_path = "MOF5"  # path to folder where CIFs will be stored
write_site_props = True  # if site properties should be written to CIF
only_ddec_charge = True  # set to True if you only want _atom_site_charge flags
# ------Settings------#

# Make new folder to store CIFs
os.makedirs(cif_folder_path, exist_ok=True)

# Read in structure data
with open(struct_json_path, "r") as f:
    qmof_struct_data = json.load(f)

# Loop over structures and write each one out to a CIF
TARGET_ID = "qmof-a2d95c3"
for entry in qmof_struct_data:
    if entry["qmof_id"] != TARGET_ID:
        continue

    qmof_id = entry["qmof_id"]  # name for CIF
    print(f"Writing {qmof_id}")

    struct = Structure.from_dict(entry["structure"])  # Pymatgen structure
    cif_path = os.path.join(cif_folder_path, f"{qmof_id}.cif")  # path to write CIF
    struct.to(filename=cif_path)  # write initial CIF

    if write_site_props:
        properties = dict(sorted(struct.site_properties.items()))  # fetch site properties
        new_cif_lines = []
        i = 0
        prop_lines = False

        with open(cif_path, "r") as f:
            for line in f:
                if "_atom_site_occupancy" in line:
                    new_cif_lines.append(line)
                    if only_ddec_charge:
                        new_cif_lines.append(" _atom_site_charge\n")
                    else:
                        for key in properties.keys():
                            new_cif_lines.append(f" _atom_site_{key}\n")
                    prop_lines = True
                    continue

                if i == len(struct):
                    prop_lines = False

                if prop_lines:
                    # Modify the atomic coordinate lines to include the property values
                    new_line = line.strip()
                    if only_ddec_charge:
                        new_line += f"  {properties['pbe_ddec_charge'][i]}"
                    else:
                        for value_sets in properties.values():
                            new_line += f"  {value_sets[i]}"
                    new_cif_lines.append(new_line + "\n")
                    i += 1
                else:
                    new_cif_lines.append(line)

        # Overwrite CIF with the modified lines
        with open(cif_path, "w") as f:
            f.writelines(new_cif_lines)