from pathlib import Path
import sys


def yield_chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def prepare(file_dict: dict):
    # read input idx file
    with open(file_dict["idx_file"], "r") as f:
        idx_file = f.readlines()

    # make a dict out of input idx file
    idx_dict = {}
    key = "header"
    for line in idx_file:
        if line.strip().startswith("["):
            key = line.strip(" []\t\n")
            idx_dict[key] = []
        else:
            idx_dict[key].extend(line.split())

    # create new idx groups and keep their names
    idx_append = []
    rotation_groups = []
    for key in ["ACE_&_CH3", "NME_&_CH3"]:
        assert idx_dict.get(
            key
        ), f"[ {key} ] section not in index file. Please add this group to the index file."
        name = key[:3]
        chunks = list(yield_chunks(idx_dict[key], 3))
        for i, chunk in enumerate(chunks):
            rotation_groups.append(f"{name}_{i}")
            idx_append.append(f"[ {name}_{i} ]")
            idx_append.append(" ".join(chunk))

    # write new idx groups to file which gets appended to index.ndx
    with open("append.ndx", "w") as f:
        f.write("\n".join(idx_append))

    # create enforced_rotation mdp section
    enforced_rotation = {
        "rotation": "Yes",
        "rot-nstrout": "1000",
        "rot-ngroups": len(rotation_groups),
        "rot-type": "rm-pf",
        "rot-vec": "0 0 1",
        "rot-rate": "0.0",
        "rot-k": "200",
        "rot-fit-method": "norm",
    }
    mdp_append = []
    mdp_append.append("\n\n; Enforced rotation general settings")
    mdp_append.append(f"rotation                 = {enforced_rotation['rotation']}")
    mdp_append.append(
        f"rot-nstrout              = {enforced_rotation['rot-nstrout']}           ; Output frequency"
    )
    mdp_append.append(
        f"rot-ngroups              = {enforced_rotation['rot-ngroups']}           ; Number of rotation groups"
    )
    mdp_append.append(f"\n; Per group settings")
    mdp_append.append(f"; rot-group: Rotation group name")
    mdp_append.append(
        f"; rot-type: Rotation potential. Can be iso, iso-pf, pm, pm-pf, rm, rm-pf, rm2, rm2-pf, flex, flex-t, flex2, flex2-t"
    )
    mdp_append.append(f"; rot-vec: Rotation vector, will get normalized")
    mdp_append.append(f"; rot-rate: Rotation rate [degree/ps]")
    mdp_append.append(f"; rot-k: Rotation force constant [kJ/(mol*nm^2)]")
    mdp_append.append(
        f"; rot-fit-method: Fitting method to determine angle of rotation group (rmsd, norm, or potential)"
    )
    for i, group in enumerate(rotation_groups):
        mdp_append.append(f"rot-group{i}              = {group}")
        mdp_append.append(
            f"rot-type{i}               = {enforced_rotation['rot-type']}"
        )
        mdp_append.append(f"rot-vec{i}                = {enforced_rotation['rot-vec']}")
        mdp_append.append(
            f"rot-rate{i}               = {enforced_rotation['rot-rate']}"
        )
        mdp_append.append(f"rot-k{i}                  = {enforced_rotation['rot-k']}")
        mdp_append.append(
            f"rot-fit-method{i}         = {enforced_rotation['rot-fit-method']}"
        )

    # create enforced_rotation mdp section
    with open("append.mdp", "w") as f:
        f.write("\n".join(mdp_append))


if __name__ == "__main__":
    try:
        file_dict = {"idx_file": Path(sys.argv[1])}
    except:
        file_dict = {"idx_file": Path("index_new.ndx")}
    prepare(file_dict)
