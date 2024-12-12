from __future__ import annotations
import numpy as np
from typing import Dict, Optional, List
import h5py
import copy
import itertools


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
    initial_biomass_concentration: float

    def _read(self, file):
        self.total_n_particle = int(file.get("/final_result/number_particles")[()])
        self.concentrations_liquid = np.array(file.get("/records/concentration_liquid"))
        self.volume_liquid = np.array(file.get("/records/volume_liquid"))
        self.initial_biomass_concentration = float(
            file.get("initial_parameters/initial_biomass_concentration")[()]
        )
        if file.get("/records/concentration_gas", None) is None:
            self.concentrations_gas = None
        else:
            self.concentrations_gas = np.array(file.get("/records/concentration_gas"))

        if file.get("/records/gas_volume", None) is None:
            self.volumes_gas = None
        else:
            self.volumes_gas = file.get("/records/gas_volume", None)

        files = file.get("files")
        self.file_list = [files[i].file.filename for i in files.keys()]
        self.t_per_flow_map = int(file.get("initial_parameters/t_per_flow_map")[()])
        self.weight = float(file.get("initial_parameters/initial_weight")[()])
        self.time = np.array(file.get("/records/time"))
        self.delta_time = float(file.get("initial_parameters/delta_time")[()])
        self.n_export = self.concentrations_liquid.shape[0]
        self.n_compartment = self.concentrations_liquid.shape[1]
        self.n_species = self.concentrations_liquid.shape[2]

    @staticmethod
    def merge(*files: MainExportResult) -> MainExportResult:
        first = files[0]
        last = files[-1]
        merged_result = copy.deepcopy(first)
        merged_result.total_n_particle = last.total_n_particle

        n_export_total = 0
        for i in files:
            n_export_total += i.n_export

        merged_result.delta_time = 0  # FIXME

        merged_result.time = np.concatenate([arg.time for arg in files], axis=0)
        merged_result.volume_liquid = np.concatenate(
            [arg.volume_liquid for arg in files], axis=0
        )
        merged_result.concentrations_liquid = np.concatenate(
            [arg.concentrations_liquid for arg in files], axis=0
        )

        if merged_result.volumes_gas is not None:
            merged_result.volumes_gas = np.concatenate(
                [arg.volumes_gas for arg in files], axis=0
            )
            merged_result.concentrations_gas = np.concatenate(
                [arg.concentrations_gas for arg in files], axis=0
            )
            pass

        merged_result.n_export = n_export_total

        return merged_result


class PartialResult(RawResults):
    probes: Optional[np.ndarray]
    particle_repartition: np.ndarray
    extra_bioparam: Optional[List[Dict[str, np.ndarray]]] = None

    def _read(self, file):
        self.particle_repartition = np.array(file.get("/records/number_particle"))
        self._read_bioparam(file)

    @staticmethod
    def _concatenate_bioparam(
        extra_bioparam: List[Optional[List[Dict[str, np.ndarray]]]],
    ) -> Optional[List[Dict[str, np.ndarray]]]:
        if not extra_bioparam or not any(extra_bioparam):
            return None

        concatenated_result = []
        for bioparam_list in extra_bioparam:
            if bioparam_list is None:
                continue

            combined_dict = {}

            for param_dict in bioparam_list:
                for key, array in param_dict.items():
                    if key not in combined_dict:
                        combined_dict[key] = array
                    else:
                        combined_dict[key] = np.concatenate(
                            (combined_dict[key], array), axis=0
                        )

            concatenated_result.append(combined_dict)

        return concatenated_result

    def _read_bioparam(self, file):
        user_export = []
        b_node = file.get("biological_model")
        dp = file.get("probes", None)
        self.probes = None
        if dp is not None:
            self.probes = np.array(dp)
            self.probes = self.probes[self.probes != 0]

        def get_dict(name: str):
            node = file.get("biological_model")[name]
            d = {}
            for key in node:
                d[key] = np.array(node[key])
            user_export.append(d)

        # get_dict("initial")

        for i_ds in range( len(b_node)):
            get_dict(f"{i_ds}")
        # try:
        #     get_dict("final")
        # except Exception as _:
        #     user_export.append(user_export[-1])  # FIXME

        self.extra_bioparam = user_export

    @staticmethod
    def merge(*files: PartialResult) -> PartialResult:
        first = files[0]
        # last = files[-1]
        merged_result = PartialResult

        merged_result.particle_repartition = np.concatenate(
            [arg.particle_repartition for arg in files], axis=0
        )

        if first.extra_bioparam is not None:
            merged_result.extra_bioparam = []
            for i in range(0, len(files)):
                merged_result.extra_bioparam = PartialResult._concatenate_bioparam(
                    [file.extra_bioparam for file in files]
                )

        return merged_result


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
    def time(self) -> np.ndarray:
        return self._main.time

    @staticmethod
    def merge(*files: Results) -> Results:
        # first = files[0]
        # last = files[-1]
        main = MainExportResult.merge(*[f.main for f in files])
        total_repartion = np.concatenate([arg.total_repartion for arg in files], axis=0)
        # tot_par = PartialResult.merge([args.])

        return Results(main, files[0]._partial, total_repartion)


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
