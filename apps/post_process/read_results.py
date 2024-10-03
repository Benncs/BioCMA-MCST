import numpy as np
from typing import Dict, Optional, List
import h5py


class RawResults:
    def _read(self, _file):
        raise Exception("Not implemented yet")


class MainExportResult(RawResults):
    total_n_particle: int
    time: np.ndarray
    concentrations_liquid: np.ndarray
    volume_liquid: np.ndarray
    concentrations_gas: Optional[np.ndarray]
    volumes_gas: Optional[np.ndarray]
    file_list: List[str]
    t_per_flow_map: int
    weight: float
    n_compartment: int
    n_species: int
    n_export: int
    delta_time: float
    initial_biomass_concentration:float

    def _read(self, file):
        self.total_n_particle = int(file.get("/final_result/number_particles")[()])
        self.concentrations_liquid = np.array(file.get("/records/concentration_liquid"))
        self.volume_liquid = np.array(file.get("/records/volume_liquid"))
        self.initial_biomass_concentration = float(file.get("initial_parameters/initial_biomass_concentration")[()])
        self.concentrations_gas = file.get("/records/concentration_gas",None)
        self.volumes_gas = file.get("/records/gas_volume",None)
        files = file.get("files")
        self.file_list = [files[i].file.filename for i in files.keys()]
        self.t_per_flow_map = int(file.get("initial_parameters/t_per_flow_map")[()])
        self.weight = float(file.get("initial_parameters/initial_weight")[()])
        self.time = np.array(file.get("/records/time"))
        self.delta_time = float(file.get("initial_parameters/delta_time")[()])
        self.n_export = self.concentrations_liquid.shape[0]
        self.n_compartment = self.concentrations_liquid.shape[1]
        self.n_species = self.concentrations_liquid.shape[2]


class PartialResult(RawResults):
    probes: np.ndarray
    particle_repartition: np.ndarray
    extra_bioparam: Optional[List[Dict[str, np.ndarray]]] = None

    def _read(self, file):
        self.particle_repartition = file.get("/records/number_particle", None)
        self._read_bioparam(file)

    def _read_bioparam(self, file):
        user_export = []
        b_node = file.get("biological_model")


        def get_dict(name:str):
            node = file.get("biological_model")[name]
            d = {}
            for key in node:
                d[key] = np.array(node[key])
            user_export.append(d)

        get_dict("initial")

        for i_ds in range(1,len(b_node)-1):
            get_dict(f"{i_ds}")
        # for ds_name in file.get("biological_model"):
        #     if(ds_name!="initial" and ds_name!="final"):
        #         print(ds_name)
        #         get_dict(ds_name)

        get_dict("final")
        self.extra_bioparam = user_export


class Results:
    def __init__(
        self,
        main: MainExportResult,
        partials: List[PartialResult],
        repartition: np.ndarray,
    ) -> None:
        self._main = main
        self._partial = partials
        self.total_repartion = repartition  # FIXME

    @property
    def partial(self) -> List[PartialResult]:
        return self._partial

    @property
    def main(self) -> MainExportResult:
        return self._main

    @property
    def time(self)->np.ndarray:
        return self._main.time


def import_results(file_name: str):
    main = MainExportResult()

    with h5py.File(file_name, "r") as file:
        main._read(file)
    

    partials = [PartialResult() for _ in range(len(main.file_list))]
    total_particle_repatition = np.zeros((main.n_export, main.n_compartment))
    for i_file, file_name in enumerate(main.file_list):
        current = partials[i_file]
        with h5py.File(file_name, "r") as file:
            current._read(file)
            total_particle_repatition += current.particle_repartition

    return Results(main, partials, total_particle_repatition)
