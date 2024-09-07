from fireworks import Firework
from QSFlow.workflows.Gromacs import *


# Copyright 2021-2022, University of Kentucky


class Ligpargen_FW(Firework):
    def __init__(
        self,
        name=None,
        priority=None,
        smiles=None,
        con=None,
        Type=None,
        di=None,
        parents=None,
        mult=None,
        **kwargs,
    ):
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "smiles": smiles,
                "name": name,
                "charge": con,
                "Type": Type,
            }
            if priority
            else {
                "_category": "gromacs",
                "smiles": smiles,
                "name": name,
                "charge": con,
                "Type": Type,
                "dir": di,
                **kwargs,
            }
        )  # # this is passed in as fw_spec when you do fw_spec.get()
        # this is retrived
        t = [Ligpargen(name=name, smile=smiles, charge=con, Type=Type, **kwargs)]
        super(Ligpargen_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class Pack_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        solute_name=None,
        solvent_name=None,
        solvent_smiles=None,
        solute_smiles=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        titration_constant=1,
        intial=False,
        own_path=None,
        multi=None,
        charge=None,
        **kwargs,
    ):
        if solute_smiles is None:
            solute_smiles = []
        if solvent_smiles is None:
            solvent_smiles = []
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [
            MDPrep(
                solute_name=solute_name,
                solvent_name=solvent_name,
                x=x,
                y=y,
                z=z,
                di=di,
                conmatrix=conmatrix,
                den=den,
                key=key,
                solvent_smiles=solvent_smiles,
                solute_smiles=solute_smiles,
                titration_constant=titration_constant,
                inital_sys=intial,
                own_path=own_path,
                multi=multi,
                charge=charge,
            )
        ]  ## this is passed for self, self.get() gets this
        super(Pack_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class Titrate(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        solute_name=None,
        solvent_name=None,
        solvent_smiles=None,
        solute_smiles=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        titration_list=None,
        intial=False,
        own_path=None,
        **kwargs,
    ):
        if solute_smiles is None:
            solute_smiles = []
        if solvent_smiles is None:
            solvent_smiles = []
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        print(f"in the details this is the titration list{titration_list}")
        t = [
            TitrationChargeScaler(
                solute_name=solute_name,
                solvent_name=solvent_name,
                x=x,
                y=y,
                z=z,
                di=di,
                conmatrix=conmatrix,
                den=den,
                key=key,
                solvent_smiles=solvent_smiles,
                solute_smiles=solute_smiles,
                titration_list=titration_list,
                inital_sys=intial,
                own_path=own_path,
            )
        ]  ## this is passed for self, self.get() gets this
        super(Titrate, self).__init__(t, parents=parents, spec=spec, name=name)


class EM_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [EnergyMinimization(key=key, **kwargs)]
        super(EM_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class NVT_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solute_name is None:
            solute_name = []
        if solvent_name is None:
            solvent_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [NVT(key=key, **kwargs)]
        super(NVT_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class NPT_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [NPT(key=key, **kwargs)]  ## this is passed for self, self.get() gets this
        super(NPT_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class Density_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [
            Density(key=key, **kwargs)
        ]  ## this is passed for self, self.get() gets this
        super(Density_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class TR_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solute_name is None:
            solute_name = []
        if solvent_name is None:
            solvent_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [
            trj_corrector(key=key, **kwargs)
        ]  ## this is passed for self, self.get() gets this
        super(TR_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class Index_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [Index(key=key, **kwargs)]  ## this is passed for self, self.get() gets this
        super(Index_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class RES_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [
            residue(key=key, **kwargs)
        ]  ## this is passed for self, self.get() gets this
        super(RES_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class RDF_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [rdf(key=key, **kwargs)]  ## this is passed for self, self.get() gets this
        super(RDF_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class CORD_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solvent_name is None:
            solvent_name = []
        if solute_name is None:
            solute_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [cord(key=key, **kwargs)]  ## this is passed for self, self.get() gets this
        super(CORD_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class Check_FW(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        mm=None,
        path_to_folder=None,
        **kwargs,
    ):
        if solute_name is None:
            solute_name = []
        if solvent_name is None:
            solvent_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                "MM": mm,
                "path_to_folder": path_to_folder,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [
            Den_checker(key=key, **kwargs)
        ]  ## this is passed for self, self.get() gets this
        print(f"in Details den {den}")
        super(Check_FW, self).__init__(t, parents=parents, spec=spec, name=name)


class key_GEN(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solute_name is None:
            solute_name = []
        if solvent_name is None:
            solvent_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [key_gen(**kwargs)]  ## this is passed for self, self.get() gets this
        super(key_GEN, self).__init__(t, parents=parents, spec=spec, name="key_gen")


class Plotter(Firework):
    def __init__(
        self,
        name=None,
        parents=None,
        priority=None,
        name_tag="",
        solute_name=None,
        solvent_name=None,
        x=None,
        y=None,
        z=None,
        di=None,
        conmatrix=None,
        den=None,
        key=None,
        **kwargs,
    ):
        if solute_name is None:
            solute_name = []
        if solvent_name is None:
            solvent_name = []
        spec = (
            {
                "_category": "processing",
                "_priority": priority,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "dir": di,
                **kwargs,
            }
            if priority
            else {
                "_category": "gromacs",
                "dir": di,
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "x": x,
                "y": y,
                "z": z,
                "conmatrix": conmatrix,
                "den": den,
                **kwargs,
            }
        )  ## this is passed in as fw_spec when you do fw_spec.get() this is retrived
        t = [
            Graph_plotter(This_key=key, **kwargs)
        ]  ## this is passed for self, self.get() gets this
        super(Plotter, self).__init__(
            t, parents=parents, spec=spec, name=f"plotter{key}"
        )
