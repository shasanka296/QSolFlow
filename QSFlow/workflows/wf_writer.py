from fireworks import Workflow
from QSFlow.workflows.FW_classes import *
from QSFlow.workflows.envwf import meta_dir


# Copyright 2021, University of Kentucky


regula=[]

def Gen_name_dic(number_of_systems, **kwargs):
    names = {}
    for i in range(number_of_systems):
        names[f"names{i + 1}"] = kwargs.get(f"WF_name{i + 1}")
    return names


def gen_ligpargen(number_of_systems, name_dic, **kwargs):
    ligpargen_fws = []
    name_dic = name_dic
    for typ, name, smiles, charge in zip(
        kwargs.get("type_list"),
        kwargs.get("name_list"),
        kwargs.get("smiles_list"),
        kwargs.get("charge_list"),
    ):
        ligpargen_fws.append(
            Ligpargen_FW(
                name=name,
                smiles=smiles,
                con=int(charge),
                Type=typ,
                di=meta_dir,
                **kwargs,
            )
        )
        print(f"appending {typ, name, smiles, charge}")
        print(
            f'the total shit{kwargs.get("type_list"), kwargs.get("name_list"),kwargs.get("smiles_list"), kwargs.get("charge_list")}'
        )
        for i in range(number_of_systems):
            if str((i + 1)) in name:
                name_dic[f"names{i + 1}"] = name_dic[f"names{i + 1}"] + f"_{name}"
    print(f"in the fuction lig_fws{ligpargen_fws}")
    return ligpargen_fws, name_dic



def edit_name_dic(name_dic, number_of_systems, **kwargs):
    name_dic = name_dic
    for name in kwargs.get("name_list"):
        for i in range(number_of_systems):
            if str((i + 1)) in name:
                name_dic[f"names{i + 1}"] = name_dic[f"names{i + 1}"] + f"_{name}"
    return name_dic


def make_pac(
    number_of_systems, name_dic, firework_dic, ligpargen_fws, key_mat, **kwargs
):
    fw_pack_keys = []
    for i in range(number_of_systems):
        fw_pack_key = f"fw_pack{i + 1}"
        fw_pack_keys.append(fw_pack_key)
        firework_dic[fw_pack_key] = Pack_FW(
            name=name_dic[f"names{i + 1}"] + "pack",
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
            den=kwargs.get(f"den{i + 1}"),
            key=key_mat[i],
            intial=kwargs.get("inital_sys"),
            own_path=kwargs.get("own_path"),
            multi=kwargs.get(f"multiplicity{i + 1}"),
            charge=kwargs.get(f"charge{i + 1}"),
        )

    return fw_pack_keys


def make_em(number_of_systems, name_dic, firework_dic, pack_keys_or_titration, key_mat,IST=False, **kwargs):
    fw_em_keys = []
    for i in range(number_of_systems):
        fw_pack_key = pack_keys_or_titration if IST else pack_keys_or_titration[i]
        fw_em_key = f"fw_em{i + 1}"
        fw_em_keys.append(fw_em_key)
        firework_dic[fw_em_key] = EM_FW(
            name=name_dic[f"names{i + 1}"] + "EM",
            parents=fw_pack_key if IST else firework_dic[fw_pack_key] ,
            key=key_mat[i],
            **kwargs,
        )
    return fw_em_keys


def make_nvt(number_of_systems, name_dic, firework_dic, em_key, key_mat, **kwargs):
    fw_nvt_keys = []
    for i in range(number_of_systems):
        fw_em_key = em_key[i]
        fw_nvt_key = f"fw_nvt{i + 1}"
        fw_nvt_keys.append(fw_nvt_key)
        firework_dic[fw_nvt_key] = NVT_FW(
            name=name_dic[f"names{i + 1}"] + "NVT",
            parents=firework_dic[fw_em_key],
            key=key_mat[i],
            **kwargs,
        )
    return fw_nvt_keys


def make_npt(number_of_systems, name_dic, firework_dic, nvt_key, key_mat, **kwargs):
    fw_npt_keys = []
    for i in range(number_of_systems):
        fw_nvt_key = nvt_key[i]
        fw_npt_key = f"fw_npt{i + 1}"
        fw_npt_keys.append(fw_npt_key)
        firework_dic[fw_npt_key] = NPT_FW(
            name=name_dic[f"names{i + 1}"] + "NPT",
            parents=firework_dic[fw_nvt_key],
            key=key_mat[i],
            **kwargs,
        )
    return fw_npt_keys


def make_den(number_of_systems, name_dic, firework_dic, npt_key, key_mat, **kwargs):
    fw_den_keys = []
    for i in range(number_of_systems):
        fw_den_key = f"fw_den{i + 1}"
        fw_den_keys.append(fw_den_key)
        fw_npt_key = npt_key[i]
        firework_dic[fw_den_key] = Density_FW(
            name=name_dic[f"names{i + 1}"] + "DEN",
            parents=firework_dic[fw_npt_key],
            key=key_mat[i],
            **kwargs,
        )
    return fw_den_keys


def make_prod(
    number_of_systems,
    name_dic,
    firework_dic,
    den_key,
    key_mat,
    molar_masses,
    densities,
    **kwargs,
):
    fw_prod_keys = []
    for i in range(number_of_systems):
        fw_prod_key = f"fw_prod{i + 1}"
        fw_prod_keys.append(fw_prod_key)
        fw_den_key = den_key[i]
        firework_dic[fw_prod_key] = Check_FW(
            name=name_dic[f"names{i + 1}"] + "prod",
            parents=firework_dic[fw_den_key],
            key=key_mat[i],
            den=densities[i],
            mm=molar_masses[i],
            **kwargs,
        )
    return fw_prod_keys


def make_trj(number_of_systems, name_dic, firework_dic, prod_key, key_mat, **kwargs):
    fw_trj_keys = []
    for i in range(number_of_systems):
        fw_trj_key = f"fw_trj{i + 1}"
        fw_trj_keys.append(fw_trj_key)
        fw_prod_key = prod_key[i]
        firework_dic[fw_trj_key] = TR_FW(
            name=name_dic[f"names{i + 1}"] + "TRJ",
            parents=firework_dic[fw_prod_key],
            key=key_mat[i],
            **kwargs,
        )
    return fw_trj_keys


def make_index(number_of_systems, name_dic, firework_dic, trj_key, key_mat, **kwargs):
    fw_index_keys = []
    for i in range(number_of_systems):
        fw_index_key = f"fw_index{i + 1}"
        fw_index_keys.append(fw_index_key)
        fw_trj_key = trj_key[i]
        firework_dic[fw_index_key] = Index_FW(
            name=name_dic[f"names{i + 1}"] + "Index",
            parents=firework_dic[fw_trj_key],
            key=key_mat[i],
            **kwargs,
        )
    return fw_index_keys


def make_resi(number_of_systems, name_dic, firework_dic, index_key, key_mat, **kwargs):
    fw_resi_keys = []
    for i in range(number_of_systems):
        fw_resi_key = f"fw_resi{i + 1}"
        fw_resi_keys.append(fw_resi_key)
        fw_index_key = index_key[i]
        firework_dic[fw_resi_key] = RES_FW(
            name=name_dic[f"names{i + 1}"] + "RESI",
            parents=firework_dic[fw_index_key],
            key=key_mat[i],
            **kwargs,
        )
    return fw_resi_keys


def make_rdf(number_of_systems, name_dic, firework_dic, resi_key, key_mat, **kwargs):
    fw_rdf_keys = []
    for i in range(number_of_systems):
        fw_rdf_key = f"fw_rdf{i + 1}"
        fw_rdf_keys.append(fw_rdf_key)
        fw_resi_key = resi_key[i]
        firework_dic[fw_rdf_key] = RDF_FW(
            name=name_dic[f"names{i + 1}"] + "RDF",
            parents=firework_dic[fw_resi_key],
            key=key_mat[i],
            **kwargs,
        )
    return fw_rdf_keys


def make_cord(number_of_systems, name_dic, firework_dic, rdf_key, key_mat, **kwargs):
    fw_cord_keys = []
    for i in range(number_of_systems):
        fw_cord_key = f"fw_rdf_{i + 1}"
        fw_cord_keys.append(fw_cord_key)
        fw_rdf_key = rdf_key[i]
        firework_dic[fw_cord_key] = CORD_FW(
            name=name_dic[f"names{i + 1}"] + "CORD",
            parents=firework_dic[fw_rdf_key],
            key=key_mat[i],
            **kwargs,
        )
    return fw_cord_keys


def matrix_of_titration_maker(
    titration_matrix,
    name_dic,
    outer_sys,
    key_mat,
    number_of_titrations,
    fire_workdir,
    packkey,
    **kwargs,
):
    for i in range(outer_sys):
        Index_key_to_pull = i * number_of_titrations
        titration_matrix.append(
            Titrate(
                name=name_dic[f"names{Index_key_to_pull + 1}"] + "Titrate",
                parents=fire_workdir[packkey[i]],
                solute_name=kwargs.get(f"solute_name{i + 1}"),
                solvent_name=kwargs.get(f"solvent_name{i + 1}"),
                solute_smiles=kwargs.get(f"solute_smiles{i + 1}"),
                solvent_smiles=kwargs.get(f"solvent_smiles{i + 1}"),
                x=kwargs.get(f"x{i + 1}"),
                y=kwargs.get(f"y{i + 1}"),
                z=kwargs.get(f"z{i + 1}"),
                di=meta_dir,
                conmatrix=kwargs.get(f"conmatrix{i + 1}"),
                den=kwargs.get(f"den{i + 1}"),
                key=key_mat[Index_key_to_pull],
                titration_list=kwargs.get("titartion_list"),
                intial=kwargs.get("inital_sys"),
                own_path=kwargs.get("own_path"),
            )
        )




def make_the_simulation(
    number_of_systems,
    name_dic,
    fire_workdir,
    key_mat,
    molarmasses,
    densities,
    ligpargen_fws,
    ist=False,
    outer_system=None,
    number_of_titrations=None,
    titration_mat=None,
    **kwargs,
):
    fw_pack_keys = make_pac(
        outer_system if ist else number_of_systems, name_dic, fire_workdir, ligpargen_fws, key_mat, **kwargs
    )
    if ist:
        matrix_of_titration_maker(
            titration_mat,
            name_dic,
            outer_system,
            key_mat,
            number_of_titrations,
            fire_workdir,
            fw_pack_keys,
            **kwargs,
        )
    em_keys = make_em(
        number_of_systems, name_dic, fire_workdir, titration_mat if ist else fw_pack_keys , key_mat,IST=ist, **kwargs
    )
    nvt_key = make_nvt(
        number_of_systems, name_dic, fire_workdir, em_keys, key_mat, **kwargs
    )
    npt_key = make_npt(
        number_of_systems, name_dic, fire_workdir, nvt_key, key_mat, **kwargs
    )
    den_key = make_den(
        number_of_systems, name_dic, fire_workdir, npt_key, key_mat, **kwargs
    )
    prod_key = make_prod(
        number_of_systems,
        name_dic,
        fire_workdir,
        den_key,
        key_mat,
        molarmasses,
        densities,
        **kwargs,
    )
    trj_key = make_trj(
        number_of_systems, name_dic, fire_workdir, prod_key, key_mat, **kwargs
    )
    index_key = make_index(
        number_of_systems, name_dic, fire_workdir, trj_key, key_mat, **kwargs
    )
    resi_key = make_resi(
        number_of_systems, name_dic, fire_workdir, index_key, key_mat, **kwargs
    )
    rdf_key = make_rdf(
        number_of_systems, name_dic, fire_workdir, resi_key, key_mat, **kwargs
    )
    cord_key = make_cord(number_of_systems, name_dic, fire_workdir, rdf_key, key_mat, **kwargs)
    regula.extend(fire_workdir.values())



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
        molarmasses = [kwargs.get(f"MM{i + 1}") for i in range(number_of_systems)]
        densities = [kwargs.get(f"den{i + 1}") for i in range(number_of_systems)]
        if not kwargs.get("inital_sys"):
            ligpargen_fws, name_dic = gen_ligpargen(
                number_of_systems, name_dic, **kwargs
            )
            print(f"HERE IS THE LIG_FWS{ligpargen_fws}")
        if kwargs.get("inital_sys"):
            name_dic = edit_name_dic(name_dic, number_of_systems, **kwargs)
        make_the_simulation(
            number_of_systems,
            name_dic,
            fire_workdir,
            key_mat,
            molarmasses,
            densities,
            ligpargen_fws,
            **kwargs,
        )

    fw_for_plotting = []
    if titration:
        number_of_titrations = len(kwargs.get("titartion_list"))
        outer_system = int(number_of_systems / number_of_titrations)
        if not kwargs.get("inital_sys"):
            ligpargen_fws, name_dic = gen_ligpargen(
                number_of_systems, name_dic, **kwargs
            )
        if kwargs.get("inital_sys"):
            name_dic = edit_name_dic(name_dic, number_of_systems, **kwargs)


        densities = [
            kwargs.get(f"den{i + 1}")
            for i in range(outer_system)
            for _ in range(number_of_titrations)
        ]
        molarmasses = [
            kwargs.get(f"MM{i + 1}")
            for i in range(outer_system)
            for _ in range(number_of_titrations)
        ]
        make_the_simulation(
            number_of_systems,
            name_dic,
            fire_workdir,
            key_mat,
            molarmasses,
            densities,
            ligpargen_fws,
            ist=True,
            outer_system=outer_system,
            number_of_titrations=number_of_titrations,
            titration_mat=matrix_of_titration, 
            **kwargs,
        )
        for i in range(outer_system):
            Index_key_to_pull = i * number_of_titrations
            fw_graph_key = f"Graph{i + 1}"
            fw_for_plotting.append(
                Plotter(
                    name=fw_graph_key,
                    parents=regula,
                    key=key_mat[Index_key_to_pull],
                    path_to_folder=os.path.join(meta_dir, "InputGrofiles"),
                    **kwargs,
                )
            )

    key_fw = key_GEN(**kwargs, parents=fw_for_plotting if titration else regula)
    fws = [key_fw] + ligpargen_fws + regula
    fws.extend(matrix_of_titration)
    fws.extend(fw_for_plotting)
    print(fws)
    print(f"this is regular {regula}")
    wf = Workflow(fws, name=kwargs.get("populate_name"))

    return wf
