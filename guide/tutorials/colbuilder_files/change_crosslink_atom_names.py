def main(run_info):
    # Imports:
    from functions_v4 import get_data_from_file, store_linelist_to_file 
    
    # Extracting data
    top_all, top_array = get_data_from_file(run_info["path_top_file"])
    gro_all, gro_array = get_data_from_file(run_info["path_gro_file"])

    # Dictionaries of corsslinks
    LYX = {"CB":"CBYX", "CG":"CGYX", "CD":"CDYX", "CE":"CEYX"}
    LY3 = {"CB":"CBY3", "CG":"CGY3"}
    LY2 = {"CB":"CBY2"}
    PYD = {"LYX":LYX, "LY3":LY3, "LY2":LY2}
    
    L5Y = {"CB":"CB5Y", "CG":"CG5Y", "CD":"CD5Y", "CE":"CE5Y"}
    L4Y = {"CB":"CB4Y", "CG":"CG4Y", "CD":"CD4Y", "CE":"CE4Y"}
    HLKNL = {"L5Y":L5Y, "L4Y":L4Y}
    
    crosslinks = {"PYD":PYD, "HLKNL":HLKNL}
    
    # Renaming atoms #Could be multithreaded?
    for file in [top_all, gro_all]:
        for _, crosslink in crosslinks.items():
            for i in range(len(file)):    
                for residue_name, residue in crosslink.items():
                    for atom_name, new_atom_name in residue.items():
                            if (residue_name in file[i] and atom_name in file[i]):
                                file[i] = file[i].replace(("  " + atom_name), new_atom_name) #adding two new letter
        
    
    #Saving modified gro and top files
    store_linelist_to_file(top_all, "top_with_renamed_atoms.top")
    store_linelist_to_file(gro_all, "gro_with_renamed_atoms.gro")
    
if __name__ == "__main__":
    import sys
    run_info = {'path_top_file':sys.argv[1],'path_gro_file':sys.argv[2]}
    main(run_info)
