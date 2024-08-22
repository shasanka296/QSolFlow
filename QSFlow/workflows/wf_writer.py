from fireworks import Workflow
from QSFlow.workflows.FW_classes import *
from QSFlow.workflows.envwf import meta_dir


# Copyright 2021, University of Kentucky

def Gen_name_dic(number_of_systems, **kwargs):
    names = {}
    for i in range(number_of_systems):
        names[f"names{i + 1}"] = kwargs.get(f"WF_name{i + 1}")
    return names


def gen_ligpargen(number_of_systems, name_dic, **kwargs):
    ligpargen_fws = []
    name_dic = name_dic
    for typ, name, smiles, charge in zip(kwargs.get("type_list"), kwargs.get("name_list"),
                                         kwargs.get("smiles_list"),
                                         kwargs.get("charge_list")):
        ligpargen_fws.append(
            Ligpargen_FW(name=name, smiles=smiles, con=int(charge), Type=typ, di=meta_dir,
                         **kwargs))
        for i in range(number_of_systems):
            if str((i + 1)) in name:
                name_dic[f"names{i + 1}"] = name_dic[f"names{i + 1}"] + f"_{name}"
    return ligpargen_fws, name_dic


def edit_name_dic(name_dic, number_of_systems, **kwargs):
    name_dic = name_dic
    for name in kwargs.get("name_list"):
        for i in range(number_of_systems):
            if str((i + 1)) in name:
                name_dic[f"names{i + 1}"] = name_dic[f"names{i + 1}"] + f"_{name}"
    return name_dic


def md_wf(**kwargs):
    key_dic = kwargs.get("key_dic")
    print(key_dic)
    key_mat = []
    for i, j in key_dic.items():
        key_mat.append(j)

    number_of_systems = int(kwargs.get("num_systems"))
    name_dic = Gen_name_dic(number_of_systems, **kwargs)
    titration = kwargs.get("is_titration")
    fire_workdir = {}
    ligpargen_fws = []
    matrix_of_titration = []

    if not titration:

        if not kwargs.get("inital_sys"):
            ligpargen_fws, name_dic = gen_ligpargen(number_of_systems, name_dic, **kwargs)

        if kwargs.get("inital_sys"):
            name_dic = edit_name_dic(name_dic, number_of_systems, **kwargs)

        for i in range(number_of_systems):
            fw_pack_key = f"fw_pack{i + 1}"
            fire_workdir[fw_pack_key] = Pack_FW(
                name=name_dic[f"names{i + 1}"] + 'pack',
                parents=ligpargen_fws,
                solute_name=kwargs.get(f"solute_name{i + 1}"),
                solvent_name=kwargs.get(f"solvent_name{i + 1}"),
                solute_smiles=kwargs.get(f"solute_smiles{i + 1}"),
                solvent_smiles=kwargs.get(f"solvent_smiles{i + 1}"),
                x=kwargs.get(f"x{i + 1}"),
                y=kwargs.get(f"y{i + 1}"),
                z=kwargs.get(f"z{i + 1}"),
                di=meta_dir,
                conmatrix=kwargs.get(f"conmatrix{i + 1}"),
                den=kwargs.get(f"den{i + 1}"), key=key_mat[i], intial=kwargs.get("inital_sys"),
                own_path=kwargs.get("own_path"),
                multi=kwargs.get(f'multiplicity{i + 1}'), charge=kwargs.get(f'charge{i + 1}')
            )

            fw_em_key = f"fw_em{i + 1}"
            fire_workdir[fw_em_key] = EM_FW(
                name=name_dic[f"names{i + 1}"] + "EM",
                parents=fire_workdir[fw_pack_key], key=key_mat[i],
                **kwargs
            )

            fw_nvt_key = f"fw_nvt{i + 1}"
            fire_workdir[fw_nvt_key] = NVT_FW(
                name=name_dic[f"names{i + 1}"] + "NVT",
                parents=fire_workdir[fw_em_key], key=key_mat[i],
                **kwargs
            )

            fw_npt_key = f"fw_npt{i + 1}"
            fire_workdir[fw_npt_key] = NPT_FW(
                name=name_dic[f"names{i + 1}"] + "NPT",
                parents=fire_workdir[fw_nvt_key], key=key_mat[i],
                **kwargs
            )

            fw_den_key = f"fw_den{i + 1}"
            fire_workdir[fw_den_key] = Density_FW(
                name=name_dic[f"names{i + 1}"] + "DEN",
                parents=fire_workdir[fw_npt_key], key=key_mat[i],
                **kwargs
            )

            fw_prod_key = f"fw_prod{i + 1}"
            fire_workdir[fw_prod_key] = Check_FW(
                name=name_dic[f"names{i + 1}"] + "prod",
                parents=fire_workdir[fw_den_key], key=key_mat[i], den=kwargs.get(f"den{i + 1}"),
                mm=kwargs.get(f"MM{i + 1}"),
                **kwargs
            )

            fw_trj_key = f"fw_trj{i + 1}"
            fire_workdir[fw_trj_key] = TR_FW(
                name=name_dic[f"names{i + 1}"] + "TRJ",
                parents=fire_workdir[fw_prod_key], key=key_mat[i],
                **kwargs
            )

            fw_index_key = f"fw_index{i + 1}"
            fire_workdir[fw_index_key] = Index_FW(
                name=name_dic[f"names{i + 1}"] + "Index",
                parents=fire_workdir[fw_trj_key], key=key_mat[i],
                **kwargs
            )

            fw_resi_key = f"fw_resi{i + 1}"
            fire_workdir[fw_resi_key] = RES_FW(
                name=name_dic[f"names{i + 1}"] + "RESI",
                parents=fire_workdir[fw_index_key], key=key_mat[i],
                **kwargs
            )

            fw_rdf_key = f"fw_rdf{i + 1}"
            fire_workdir[fw_rdf_key] = RDF_FW(
                name=name_dic[f"names{i + 1}"] + "RDF",
                parents=fire_workdir[fw_resi_key], key=key_mat[i],
                **kwargs
            )

            fw_cord_key = f"fw_cord{i + 1}"
            fire_workdir[fw_cord_key] = CORD_FW(
                name=name_dic[f"names{i + 1}"] + "CORD",
                parents=fire_workdir[fw_rdf_key], key=key_mat[i],
                **kwargs
            )

    if titration:
        global number_of_titrations
        global outer_system
        number_of_titrations = len(kwargs.get("titartion_list"))
        outer_system = int(number_of_systems / number_of_titrations)

        if not kwargs.get("inital_sys"):
            ligpargen_fws, name_dic = gen_ligpargen(number_of_systems, name_dic, **kwargs)

        if kwargs.get("inital_sys"):
            name_dic = edit_name_dic(name_dic, number_of_systems, **kwargs)

        print(f"the name dic is {name_dic}")
        print(
            f'type list:{kwargs.get("type_list")} name_lis:{kwargs.get("name_list")} smiles_list:{kwargs.get("smiles_list")} ')

        for i in range(outer_system):
            print(key_mat)
            Index_key_to_pull = i * number_of_titrations
            fw_pack_key = f"fw_pack{i + 1}"
            fire_workdir[fw_pack_key] = Pack_FW(
                name=name_dic[f"names{Index_key_to_pull + 1}"] + 'pack',
                parents=ligpargen_fws,
                solute_name=kwargs.get(f"solute_name{i + 1}"),
                solvent_name=kwargs.get(f"solvent_name{i + 1}"),
                solute_smiles=kwargs.get(f"solute_smiles{i + 1}"),
                solvent_smiles=kwargs.get(f"solvent_smiles{i + 1}"),
                x=kwargs.get(f"x{i + 1}"),
                y=kwargs.get(f"y{i + 1}"),
                z=kwargs.get(f"z{i + 1}"),
                di=meta_dir,
                conmatrix=kwargs.get(f"conmatrix{i + 1}"),
                den=kwargs.get(f"den{i + 1}"), key=key_mat[Index_key_to_pull], intial=kwargs.get("inital_sys"),
                own_path=kwargs.get("own_path"),
                multi=kwargs.get(f'multiplicity{i + 1}'), charge=kwargs.get(f'charge{i + 1}')
            )
            matrix_of_titration.append(Titrate(
                name=name_dic[f"names{Index_key_to_pull + 1}"] + 'Titrate',
                parents=fire_workdir[fw_pack_key],
                solute_name=kwargs.get(f"solute_name{i + 1}"),
                solvent_name=kwargs.get(f"solvent_name{i + 1}"),
                solute_smiles=kwargs.get(f"solute_smiles{i + 1}"),
                solvent_smiles=kwargs.get(f"solvent_smiles{i + 1}"),
                x=kwargs.get(f"x{i + 1}"),
                y=kwargs.get(f"y{i + 1}"),
                z=kwargs.get(f"z{i + 1}"),
                di=meta_dir,
                conmatrix=kwargs.get(f"conmatrix{i + 1}"),
                den=kwargs.get(f"den{i + 1}"), key=key_mat[Index_key_to_pull],
                titration_list=kwargs.get("titartion_list"),
                intial=kwargs.get("inital_sys"),
                own_path=kwargs.get("own_path")
            ))
            print(f'in the workflow this is what is being passed for titration list{kwargs.get("titartion_list")}')
        densities = [kwargs.get(f"den{i + 1}") for i in range(outer_system) for _ in range(number_of_titrations)]
        print(f"what is passed in as  density {densities[0]}")
        print(f"waht is passed for deinstiy matrix{densities}")
        molarmasses = [kwargs.get(f"MM{i + 1}") for i in range(outer_system) for _ in range(number_of_titrations)]  #

        for i in range(number_of_systems):
            print("this is key")
            print(key_mat)

            fw_em_key = f"fw_em{i + 1}"
            fire_workdir[fw_em_key] = EM_FW(
                name=name_dic[f"names{i + 1}"] + "EM",
                parents=matrix_of_titration, key=key_mat[i],
                **kwargs
            )

            fw_nvt_key = f"fw_nvt{i + 1}"
            fire_workdir[fw_nvt_key] = NVT_FW(
                name=name_dic[f"names{i + 1}"] + "NVT",
                parents=fire_workdir[fw_em_key], key=key_mat[i],
                **kwargs
            )

            fw_npt_key = f"fw_npt{i + 1}"
            fire_workdir[fw_npt_key] = NPT_FW(
                name=name_dic[f"names{i + 1}"] + "NPT",
                parents=fire_workdir[fw_nvt_key], key=key_mat[i],
                **kwargs
            )

            fw_den_key = f"fw_den{i + 1}"
            fire_workdir[fw_den_key] = Density_FW(
                name=name_dic[f"names{i + 1}"] + "DEN",
                parents=fire_workdir[fw_npt_key], key=key_mat[i],
                **kwargs
            )

            fw_prod_key = f"fw_prod{i + 1}"
            fire_workdir[fw_prod_key] = Check_FW(
                name=name_dic[f"names{i + 1}"] + "prod",
                parents=fire_workdir[fw_den_key], key=key_mat[i], den=densities[i],
                mm=molarmasses[i],
                **kwargs
            )

            fw_trj_key = f"fw_trj{i + 1}"
            fire_workdir[fw_trj_key] = TR_FW(
                name=name_dic[f"names{i + 1}"] + "TRJ",
                parents=fire_workdir[fw_prod_key], key=key_mat[i],
                **kwargs
            )

            fw_index_key = f"fw_index{i + 1}"
            fire_workdir[fw_index_key] = Index_FW(
                name=name_dic[f"names{i + 1}"] + "Index",
                parents=fire_workdir[fw_trj_key], key=key_mat[i],
                **kwargs
            )

            fw_resi_key = f"fw_resi{i + 1}"
            fire_workdir[fw_resi_key] = RES_FW(
                name=name_dic[f"names{i + 1}"] + "RESI",
                parents=fire_workdir[fw_index_key], key=key_mat[i],
                **kwargs
            )

            fw_rdf_key = f"fw_rdf{i + 1}"
            fire_workdir[fw_rdf_key] = RDF_FW(
                name=name_dic[f"names{i + 1}"] + "RDF",
                parents=fire_workdir[fw_resi_key], key=key_mat[i],
                **kwargs
            )

            fw_cord_key = f"fw_cord{i + 1}"
            fire_workdir[fw_cord_key] = CORD_FW(
                name=name_dic[f"names{i + 1}"] + "CORD",
                parents=fire_workdir[fw_rdf_key], key=key_mat[i],
                **kwargs
            )
            print(f"{i + 1} iteration")
    print("made regualr")
    regula = []
    plotter = []

    for fw, values in fire_workdir.items():
        regula.append(values)
    print(f"the lig dict {len(ligpargen_fws)} done")
    print(f"lig fw: {ligpargen_fws}")
    print(f"matrix_titration: {matrix_of_titration}")
    if titration:

        for i in range(outer_system):
            Index_key_to_pull = i * number_of_titrations
            fw_graph_key = f"Graph{i + 1}"
            plotter.append(Plotter(
                name=fw_graph_key,
                parents=regula,
                key=key_mat[Index_key_to_pull], path_to_folder=os.path.join(meta_dir, "InputGrofiles"), **kwargs
            ))

    key_fw = key_GEN(**kwargs, parents=plotter if titration else regula)
    fws = [key_fw] + ligpargen_fws + regula
    fws.extend(matrix_of_titration)
    fws.extend(plotter)
    print(fws)

    wf = Workflow(fws, name=kwargs.get("populate_name"))

    return wf
