from pathlib import Path
import numpy as np
from fireworks import Workflow
from d3tales_fw.workflows.D3TaLES_FW import *
from d3tales_fw.workflows.envwf import meta_dir
# from d3tales_api.Workflows.D3TaLES_FW import *
# from d3tales_api.Workflows.ParamSet import GausParamSet


# Copyright 2021, University of Kentucky


def d3tales_wf(paramset, identifier=None, smiles=None, wtune=True, solvent='acetonitrile',
               hf_mol_opt=False, email=None, username=None, wf_tag="", **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)

    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState
    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    opt_parents = f30 if wtune else f10
    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=opt_parents, **kwargs)
    f35 = Optimization(paramset=paramset.opt_anion1, species="anion1", parents=opt_parents, **kwargs)
    f36 = Optimization(paramset=paramset.opt_anion2, species="anion2", parents=opt_parents, **kwargs)
    f37 = Optimization(paramset=paramset.opt_cation1, species="cation1", parents=opt_parents, **kwargs)
    f38 = Optimization(paramset=paramset.opt_cation2, species="cation2", parents=opt_parents, **kwargs)

    f45 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="opt_anion1", parents=[f35],
                 **kwargs)
    f46 = Energy(paramset=paramset.energy_anion1, species="anion1", geometry="opt_groundState", parents=[f31], **kwargs)
    f47 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="opt_cation1", parents=[f37],
                 **kwargs)
    f48 = Energy(paramset=paramset.energy_cation1, species="cation1", geometry="opt_groundState", parents=[f31],
                 **kwargs)

    f51 = Energy(paramset=paramset.energy_groundState, parents=[f31], species='groundState', geometry="opt_groundState",
                 solvent=solvent, **kwargs)
    f55 = Energy(paramset=paramset.energy_anion1, parents=[f35], species='anion1', geometry="opt_anion1",
                 solvent=solvent, **kwargs)
    f56 = Energy(paramset=paramset.energy_anion2, parents=[f36], species='anion2', geometry="opt_anion2",
                 solvent=solvent, **kwargs)
    f57 = Energy(paramset=paramset.energy_cation1, parents=[f37], species='cation1', geometry="opt_cation1",
                 solvent=solvent, **kwargs)
    f58 = Energy(paramset=paramset.energy_cation2, parents=[f38], species='cation2', geometry="opt_cation2",
                 solvent=solvent, **kwargs)

    f61 = TDDFT(paramset=paramset.tddft_groundState, parents=[f31], **kwargs)
    f65 = TDDFT(paramset=paramset.tddft_anion1, species='anion1', parents=[f35], **kwargs)
    f66 = TDDFT(paramset=paramset.tddft_anion2, species='anion2', parents=[f36], **kwargs)
    f67 = TDDFT(paramset=paramset.tddft_cation1, species='cation1', parents=[f37], **kwargs)
    f68 = TDDFT(paramset=paramset.tddft_cation2, species='cation2', parents=[f38], **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f31, f35, f36, f37, f38, f45, f46, f47, f48, f51, f55, f56, f57, f58, f61, f65, f66, f67, f68, ]
    if wtune:
        fws.extend([f26, f30])

    if email:
        fws.append(EmailStart(identifier=identifier, email=email, username=username, parents=[f10]))
        fws.append(
            EmailEnd(identifier=identifier, email=email, username=username, parents=[f51, f55, f57, f61, f65, f67]))

    wf = Workflow(fws, name="{}gaus_{}".format(wf_tag, identifier or smiles))

    return wf


def d3tales_md_wf(param_file=None, **kwargs):

    # Establish calculation parameters from parm_file json file
    # param_file = param_file or os.path.join(Path(__file__).resolve().parent.parent, "parameters", 'md_gaus_parameter_file.json')
    # paramset = GausParamSet().from_json(param_file)

    key_dic = kwargs.get("key_dic")
    print(key_dic)
    key_mat = []
    for i, j in key_dic.items():
        key_mat.append(j)

    number_of_systems = int(kwargs.get("num_systems"))
    name_dic = {}
    titration=kwargs.get("is_titration")
    fire_workdic = {}
    # dft_fw = []
    folder=[]
    ligpargen_fws = []
    matrix_of_titration = []

    if not titration:
        for i in range(number_of_systems):
            name_dic[f"names{i + 1}"] = kwargs.get(f"WF_name{i + 1}")
        # if kwargs.get("own") == False and kwargs.get("inital_sys") ==False:
        #
        #     for smiles, (index,names) in zip(kwargs.get("smiles_list"), enumerate(kwargs.get("name_list"))):
        #         d = Optimization(paramset=paramset.opt_groundState,
        #                                                      species="groundState",
        #                                                      parents=folder[index-1] if index !=0 else None, s=f"{smiles}", submit=False, smiles=smiles)
        #         dft_fw.append(d)
        #         folder.append(fo(name=names,parents=d,**kwargs))

        if kwargs.get("inital_sys") == False:
            for typ, name, smiles, charge in zip(kwargs.get("type_list"), kwargs.get("name_list"),
                                                kwargs.get("smiles_list"),
                                                kwargs.get("charge_list")):
                ligpargen_fws.append(
                    Ligpargen_FW(name=name, smiles=smiles, con=int(charge), Type=typ, di=meta_dir, parents=folder,
                                 **kwargs))
                for i in range(number_of_systems):
                    if str((i + 1)) in name:
                        name_dic[f"names{i + 1}"] = name_dic[f"names{i + 1}"] + f"_{name}"
        if kwargs.get("inital_sys") == True:
            for  name in  kwargs.get("name_list"):
                for i in range(number_of_systems):
                    if str((i + 1)) in name:
                        name_dic[f"names{i + 1}"] = name_dic[f"names{i + 1}"] + f"_{name}"


        for i in range(number_of_systems):

            fw_pack_key = f"fw_pack{i + 1}"
            fire_workdic[fw_pack_key] = Pack_FW (
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
                den=kwargs.get(f"den{i + 1}"), key=key_mat[i], intial= kwargs.get("inital_sys"), own_path=kwargs.get("own_path"),
                multi=kwargs.get(f'multiplicity{i + 1}'), charge=kwargs.get(f'charge{i+1}')
            )

            fw_em_key = f"fw_em{i + 1}"
            fire_workdic[fw_em_key] = EM_FW(
                name=name_dic[f"names{i + 1}"] + "EM",
                parents=fire_workdic[fw_pack_key], key=key_mat[i],
                **kwargs
            )

            fw_nvt_key = f"fw_nvt{i + 1}"
            fire_workdic[fw_nvt_key] = NVT_FW (
                name=name_dic[f"names{i + 1}"] + "NVT",
                parents=fire_workdic[fw_em_key], key=key_mat[i],
                **kwargs
            )

            fw_npt_key = f"fw_npt{i + 1}"
            fire_workdic[fw_npt_key] = NPT_FW(
                name=name_dic[f"names{i + 1}"] + "NPT",
                parents=fire_workdic[fw_nvt_key], key=key_mat[i],
                **kwargs
            )

            fw_den_key = f"fw_den{i + 1}"
            fire_workdic[fw_den_key] = Density_FW(
                name=name_dic[f"names{i + 1}"] + "DEN",
                parents=fire_workdic[fw_npt_key], key=key_mat[i],
                **kwargs
            )

            fw_prod_key = f"fw_prod{i + 1}"
            fire_workdic[fw_prod_key] = Check_FW(
                name=name_dic[f"names{i + 1}"] + "prod",
                parents=fire_workdic[fw_den_key], key=key_mat[i], den=kwargs.get(f"den{i + 1}"),
                mm=kwargs.get(f"MM{i + 1}"),
                **kwargs
            )

            fw_trj_key = f"fw_trj{i + 1}"
            fire_workdic[fw_trj_key] = TR_FW(
                name=name_dic[f"names{i + 1}"] + "TRJ",
                parents=fire_workdic[fw_prod_key], key=key_mat[i],
                **kwargs
            )

            fw_index_key = f"fw_index{i + 1}"
            fire_workdic[fw_index_key] = Index_FW(
                name=name_dic[f"names{i + 1}"] + "Index",
                parents=fire_workdic[fw_trj_key], key=key_mat[i],
                **kwargs
            )

            fw_resi_key = f"fw_resi{i + 1}"
            fire_workdic[fw_resi_key] = RES_FW(
                name=name_dic[f"names{i + 1}"] + "RESI",
                parents=fire_workdic[fw_index_key], key=key_mat[i],
                **kwargs
            )

            fw_rdf_key = f"fw_rdf{i + 1}"
            fire_workdic[fw_rdf_key] = RDF_FW(
                name=name_dic[f"names{i + 1}"] + "RDF",
                parents=fire_workdic[fw_resi_key], key=key_mat[i],
                **kwargs
            )

            fw_cord_key = f"fw_cord{i + 1}"
            fire_workdic[fw_cord_key] = CORD_FW(
                name=name_dic[f"names{i + 1}"] + "CORD",
                parents=fire_workdic[fw_rdf_key], key=key_mat[i],
                **kwargs
            )


    if titration:
        global number_of_titrations
        global outer_system
        number_of_titrations = len(kwargs.get("titartion_list"))
        outer_system = int(number_of_systems / (number_of_titrations))
        for i in range(number_of_systems):
            name_dic[f"names{i + 1}"] = kwargs.get(f"WF_name{i + 1}")
        # if kwargs.get("own") == False:
        #
        #     for smiles, (index,names) in zip(kwargs.get("smiles_list"), enumerate(kwargs.get("name_list"))):
        #         d = Optimization(paramset=paramset.opt_groundState,
        #                                                      species="groundState",
        #                                                      parents=folder[index-1] if index !=0 else None, s=f"{smiles}", submit=False, smiles=smiles)
        #         dft_fw.append(d)
        #         folder.append(fo(name=names,parents=d,**kwargs))






        if kwargs.get("inital_sys") == False:
            for typ, name, smiles,charge in zip(kwargs.get("type_list"), kwargs.get("name_list"), kwargs.get("smiles_list"),
                     kwargs.get("charge_list") ):
                ligpargen_fws.append(Ligpargen_FW(name=name, smiles=smiles, con=int(charge), Type=typ, di=meta_dir,parents=folder, **kwargs))
                for i in range(number_of_systems):
                    if str((i + 1)) in name:
                        name_dic[f"names{i + 1}"] = name_dic[f"names{i + 1}"] + f"_{name}"
        if kwargs.get("inital_sys") == True:
            for  name in  kwargs.get("name_list"):
                for i in range(number_of_systems):
                    if str((i + 1)) in name:
                        name_dic[f"names{i + 1}"] = name_dic[f"names{i + 1}"] + f"_{name}"
        print(f"the name dic is {name_dic}")
        print(f'type list:{kwargs.get("type_list")} name_lis:{kwargs.get("name_list")} smiles_list:{ kwargs.get("smiles_list")} ')


        for i in range(outer_system):
            print(key_mat)
            Index_key_to_pull= i* number_of_titrations
            fw_pack_key = f"fw_pack{i + 1}"
            fire_workdic[fw_pack_key] = Pack_FW(
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
            matrix_of_titration.append( Titrate(
                name=name_dic[f"names{Index_key_to_pull + 1}"] + 'Titrate',
                parents=fire_workdic[fw_pack_key],
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
                titration_list=kwargs.get("titartion_list") ,
                intial= kwargs.get("inital_sys"),
                own_path=kwargs.get("own_path")
            )   )
            print(f'in the workflow this is what is being passed for titration list{kwargs.get("titartion_list")}')
        densities=[kwargs.get(f"den{i + 1}") for i in range(outer_system) for _ in range(number_of_titrations)]
        print(f"what is passed in as  density {densities[0]}")
        print(f"waht is passed for deinstiy matrix{densities}")
        molarmasses=[kwargs.get(f"MM{i + 1}") for i in range(outer_system) for _ in range(number_of_titrations)] #

        for i in range(number_of_systems):


            fw_em_key = f"fw_em{i + 1}"
            fire_workdic[fw_em_key] = EM_FW(
                name=name_dic[f"names{i + 1}"] + "EM",
                parents= matrix_of_titration, key=key_mat[i],
                **kwargs
            )

            fw_nvt_key = f"fw_nvt{i + 1}"
            fire_workdic[fw_nvt_key] = NVT_FW(
                name=name_dic[f"names{i + 1}"] + "NVT",
                parents=fire_workdic[fw_em_key], key=key_mat[i],
                **kwargs
            )

            fw_npt_key = f"fw_npt{i + 1}"
            fire_workdic[fw_npt_key] = NPT_FW(
                name=name_dic[f"names{i + 1}"] + "NPT",
                parents=fire_workdic[fw_nvt_key], key=key_mat[i],
                **kwargs
            )

            fw_den_key = f"fw_den{i + 1}"
            fire_workdic[fw_den_key] = Density_FW(
                name=name_dic[f"names{i + 1}"] + "DEN",
                parents=fire_workdic[fw_npt_key], key=key_mat[i],
                **kwargs
            )

            fw_prod_key = f"fw_prod{i + 1}"
            fire_workdic[fw_prod_key] = Check_FW(
                name=name_dic[f"names{i + 1}"] + "prod",
                parents=fire_workdic[fw_den_key], key=key_mat[i], den=densities[i],
                mm=molarmasses[i],
                **kwargs
            )

            fw_trj_key = f"fw_trj{i + 1}"
            fire_workdic[fw_trj_key] = TR_FW(
                name=name_dic[f"names{i + 1}"] + "TRJ",
                parents=fire_workdic[fw_prod_key], key=key_mat[i],
                **kwargs
            )

            fw_index_key = f"fw_index{i + 1}"
            fire_workdic[fw_index_key] = Index_FW(
                name=name_dic[f"names{i + 1}"] + "Index",
                parents=fire_workdic[fw_trj_key], key=key_mat[i],
                **kwargs
            )

            fw_resi_key = f"fw_resi{i + 1}"
            fire_workdic[fw_resi_key] = RES_FW(
                name=name_dic[f"names{i + 1}"] + "RESI",
                parents=fire_workdic[fw_index_key], key=key_mat[i],
                **kwargs
            )

            fw_rdf_key = f"fw_rdf{i + 1}"
            fire_workdic[fw_rdf_key] = RDF_FW(
                name=name_dic[f"names{i + 1}"] + "RDF",
                parents=fire_workdic[fw_resi_key], key=key_mat[i],
                **kwargs
            )

            fw_cord_key = f"fw_cord{i + 1}"
            fire_workdic[fw_cord_key] = CORD_FW(
                name=name_dic[f"names{i + 1}"] + "CORD",
                parents=fire_workdic[fw_rdf_key], key=key_mat[i],
                **kwargs
            )
            print(f"{i+1} iteration")
    print("made regualr")
    regula = []
    plotter=[]


    for fw, values in fire_workdic.items():
        regula.append(values)
    print(f"the lig dict {len(ligpargen_fws)} done")
    print(f"lig fw: {ligpargen_fws}")
    print(f"matrix_titration: {matrix_of_titration}")
    if titration:

        for i in range(outer_system):
            Index_key_to_pull = i * number_of_titrations
            fw_graph_key = f"Graph{i + 1}"
            plotter.append( Plotter(
                name=fw_graph_key,
                parents=regula,
                key=key_mat[Index_key_to_pull],path_to_folder= os.path.join(meta_dir,"InputGrofiles"), **kwargs
            ))

    key_fw = key_GEN(**kwargs, parents=plotter if titration else regula)
    fws = [key_fw] + ligpargen_fws + regula
    fws.extend(matrix_of_titration)
    fws.extend(plotter)
    fws.extend(folder)
    print(fws)

    wf = Workflow(fws, name=kwargs.get("populate_name"))


    return wf


def just_anion(paramset, identifier=None, smiles=None, wtune=False, solvent='acetonitrile',
               hf_mol_opt=False, email=None, username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)

    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState
    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)
    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=f30, **kwargs)

    opt_parents = f30 if wtune else f10

    f35 = Optimization(paramset=paramset.opt_anion1, species="anion1", parents=opt_parents, **kwargs)
    f36 = Optimization(paramset=paramset.opt_anion2, species="anion2", parents=opt_parents, **kwargs)

    f45 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="opt_anion1", parents=[f35],
                 **kwargs)
    f46 = Energy(paramset=paramset.energy_anion1, species="anion1", geometry="opt_groundState", parents=[opt_parents],
                 **kwargs)

    f55 = Energy(paramset=paramset.energy_anion1, parents=[f35], species='anion1', geometry="opt_anion1",
                 solvent=solvent, **kwargs)
    f56 = Energy(paramset=paramset.energy_anion2, parents=[f36], species='anion2', geometry="opt_anion2",
                 solvent=solvent, **kwargs)

    f65 = TDDFT(paramset=paramset.tddft_anion1, species='anion1', parents=[f35], **kwargs)
    f66 = TDDFT(paramset=paramset.tddft_anion2, species='anion2', parents=[f36], **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f35, f36, f45, f46, f55, f56, f65, f66]

    if wtune:
        fws.extend([f26, f30, f31])

    if email:
        fws.append(EmailStart(identifier=identifier, email=email, username=username, parents=[f10]))
        fws.append(EmailEnd(identifier=identifier, email=email, username=username, parents=[f55, f65, f66]))

    wf = Workflow(fws, name="anion_{}".format(identifier or smiles))

    return wf


def solv_wf(paramset, identifier=None, smiles=None, wtune=True, solvent='acetonitrile',
            hf_mol_opt=False, email=None, username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles,
                             **kwargs)

    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState
    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)

    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    opt_parents = f30 if wtune else f10
    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=opt_parents, solvent=solvent,
                       **kwargs)
    f37 = Optimization(paramset=paramset.opt_cation1, species="cation1", parents=opt_parents, solvent=solvent, **kwargs)
    f38 = Optimization(paramset=paramset.opt_cation2, species="cation2", parents=opt_parents, solvent=solvent, **kwargs)

    f47 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="solv_opt_cation1",
                 parents=[f37], solvent=solvent, **kwargs)
    f48 = Energy(paramset=paramset.energy_cation1, species="cation1", geometry="solv_opt_groundState", parents=[f31],
                 solvent=solvent, **kwargs)

    f61 = TDDFT(paramset=paramset.tddft_groundState, parents=[f31], solvent=solvent, solv_geom=True, **kwargs)
    f67 = TDDFT(paramset=paramset.tddft_cation1, species='cation1', parents=[f37], solvent=solvent, solv_geom=True,
                **kwargs)
    f68 = TDDFT(paramset=paramset.tddft_cation2, species='cation2', parents=[f38], solvent=solvent, solv_geom=True,
                **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f31, f37, f38, f47, f48, f61, f67, f68, ]
    if wtune:
        fws.extend([f26, f30])

    if email:
        fws.append(EmailStart(identifier=identifier, email=email, username=username, parents=[f10]))
        fws.append(EmailEnd(identifier=identifier, email=email, username=username, parents=[f61, f67]))

    wf = Workflow(fws, name="solv_{}".format(identifier or smiles))

    return wf


def hf_wf(paramset, identifier=None, smiles=None, hf_mol_opt=False, email=None, username=None, solvent=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState
    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    f36 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=f30, **kwargs)
    f37 = Optimization(paramset=paramset.hf_opt_groundState, species="groundState", parents=f30, name_tag='hf_',
                       **kwargs)
    f46 = Optimization(paramset=paramset.opt_cation1, species="cation1", parents=f30, **kwargs)
    f47 = Optimization(paramset=paramset.hf_opt_cation1, species="cation1", parents=f30, name_tag='hf_', **kwargs)
    f56 = Optimization(paramset=paramset.opt_cation2, species="cation2", parents=f30, **kwargs)
    f57 = Optimization(paramset=paramset.hf_opt_cation2, species="cation2", parents=f30, name_tag='hf_', **kwargs)

    f38 = Energy(paramset=paramset.energy_groundState, species="groundState", geometry="hf_opt_groundState",
                 parents=[f37], name_tag='hfgeom_', **kwargs)
    f39 = Energy(paramset=paramset.hfgeom_energy_groundState, species="groundState", geometry="opt_groundState",
                 parents=[f36], name_tag='hf_dftgeom_', **kwargs)
    f48 = Energy(paramset=paramset.energy_cation1, species="cation1", geometry="hf_opt_cation1", parents=[f47],
                 name_tag='hfgeom_', **kwargs)
    f49 = Energy(paramset=paramset.hfgeom_energy_cation1, species="cation1", geometry="opt_cation1", parents=[f46],
                 name_tag='hf_dftgeom_', **kwargs)
    f58 = Energy(paramset=paramset.energy_cation2, species="cation2", geometry="hf_opt_cation2", parents=[f57],
                 name_tag='hfgeom_', **kwargs)
    f59 = Energy(paramset=paramset.hfgeom_energy_cation2, species="cation2", geometry="opt_cation2", parents=[f56],
                 name_tag='hf_dftgeom_', **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f26, f30, f36, f37, f38, f39, f46, f47, f48, f49, f56, f57, f58, f59]
    wf = Workflow(fws, name="gaus_hf_{}".format(identifier or smiles))

    return wf


def huckaba_wf(paramset, identifier=None, smiles=None, solvent='DiMethylSulfoxide', wtune=True, hf_mol_opt=False,
               email=None, username=None, **kwargs):
    kwargs.update(dict(iop_str="0170700000"))
    # name = kwargs.get('mol_name')
    # species = "groundState" if "_p0" in name else "cation1" if "_p1" in name else "cation2"
    species = "groundState"

    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    f26 = MolOpt(paramset=paramset.opt_groundState, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)
    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=f30, **kwargs)

    opt_parents = f30 if wtune else f10

    f37 = Optimization(paramset=eval("paramset.opt_" + species), species=species, parents=opt_parents, solvent=solvent,
                       **kwargs)
    f63 = TDDFT(paramset=eval("paramset.tddft_" + species), parents=[f37], **kwargs)

    fws = [f10, f37, f63]
    if wtune:
        fws.extend([f26, f30, f31])

    wf = Workflow(fws, name="special_{}".format(identifier or smiles))
    return wf


def just_nmr(paramset, identifier=None, smiles=None, solvent='acetonitrile', hf_mol_opt=False, email=None,
             username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    f76 = Energy(paramset=paramset.nmr_groundState, parents=[f10], species='groundState', geometry="opt_groundState",
                 name_tag='nmr_', solvent=solvent, **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f76, ]
    wf = Workflow(fws, name="gaus_nmr_{}".format(identifier or smiles))

    return wf


def just_tddft(paramset, identifier=None, smiles=None, wtune=False, solvent=None,
               hf_mol_opt=False, email=None, username=None, check_if_already_run=True, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)

    mol_opt_params = paramset.hf_opt_groundState if hf_mol_opt else paramset.opt_groundState

    f26 = MolOpt(paramset=mol_opt_params, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)
    opt_parents = f30 if wtune else f10

    f31 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=opt_parents, solvent=solvent,
                       check_if_already_run=check_if_already_run, **kwargs)
    f37 = Optimization(paramset=paramset.opt_cation1, species="cation1", parents=opt_parents, solvent=solvent,
                       check_if_already_run=check_if_already_run, **kwargs)
    f38 = Optimization(paramset=paramset.opt_cation2, species="cation2", parents=opt_parents, solvent=solvent,
                       check_if_already_run=check_if_already_run, **kwargs)

    f61 = TDDFT(paramset=paramset.tddft_groundState, parents=[f31], solvent=solvent, solv_geom=True, **kwargs)
    f67 = TDDFT(paramset=paramset.tddft_cation1, species='cation1', parents=[f37], solvent=solvent, solv_geom=True,
                **kwargs)
    f68 = TDDFT(paramset=paramset.tddft_cation2, species='cation2', parents=[f38], solvent=solvent, solv_geom=True,
                **kwargs)

    # Establish fireworks in workflow
    fws = [f10, f31, f37, f38, f61, f67, f68, ]
    if wtune:
        fws.extend([f26, f30])

    if email:
        fws.append(EmailStart(identifier=identifier, email=email, username=username, parents=[f10]))
        fws.append(EmailEnd(identifier=identifier, email=email, username=username, parents=[f61, f67]))

    wf = Workflow(fws, name="tddft_{}".format(identifier or smiles))

    return wf


def just_initialize(paramset, identifier=None, smiles=None, solvent='acetonitrile', hf_mol_opt=False, email=None,
                    username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    fws = [f10, ]
    wf = Workflow(fws, name="init_{}".format(identifier or smiles))
    return wf


def dihed_rot(paramset, identifier=None, smiles=None, solvent='acetonitrile', hf_mol_opt=False, email=None,
              username=None, **kwargs):
    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    f12 = LowestEConformer(smiles=smiles, paramset=paramset.semi_empirical, parents=f10, **kwargs)
    f26 = Optimization(paramset=paramset.opt_groundState, geometry="conformer", species='groundState', parents=f12,
                       **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)
    f36 = Optimization(paramset=paramset.opt_groundState, species="groundState", parents=f30, **kwargs)

    # fws = [f10, f36, ]  # exclude omega tuning
    fws = [f10, f12, f26, f30, f36, ]  # include omega tuning

    for degree in np.linspace(start=0, stop=180, num=19):
        fw = DihedRot(paramset=paramset.opt_groundState, dihed_degree=degree, parents=f36, **kwargs)
        fws.append(fw)

    wf = Workflow(fws, name="dihedrot_{}".format(identifier or smiles))
    return wf


def specialized_job(paramset, identifier=None, smiles=None, solvent='acetonitrile', wtune=True, hf_mol_opt=False,
                    email=None, username=None, check_if_already_run=False, **kwargs):
    kwargs.update(dict(iop_str="0170700000"))
    # name = kwargs.get('mol_name')
    # species = "groundState" if "_p0" in name else "cation1" if "_p1" in name else "cation2"
    species = "groundState"

    f10 = InitializeMolecule(identifier=identifier, smiles=smiles, **kwargs)
    f26 = MolOpt(paramset=paramset.opt_groundState, parents=f10, **kwargs)
    f30 = WTuning(paramset=paramset.wtuning, parents=f26, **kwargs)

    opt_parents = f30 if wtune else f10

    f36 = Optimization(paramset=eval("paramset.opt_" + species), species=species, parents=opt_parents, **kwargs)
    f37 = Optimization(paramset=eval("paramset.opt_" + species), species=species, parents=opt_parents, solvent=solvent,
                       **kwargs)
    f63 = TDDFT(paramset=eval("paramset.tddft_" + species), parents=[f36], **kwargs)

    fws = [f10, f36, f37, f63]
    if wtune:
        fws.extend([f26, f30])

    wf = Workflow(fws, name="special_{}".format(identifier or smiles))
    return wf
