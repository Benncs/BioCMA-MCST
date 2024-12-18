import h5py
import numpy as np


def get_np(n, file):
    if n is None:
        return np.array(file.get("/records/number_particle"))
    else:
        return n + np.array(file.get("/records/number_particle"))


def probe(file, probes):
    tmp_probe = np.array(file.get("probes"))
    if len(tmp_probe) != 0:
        if probes is not None:
            return np.concatenate((probes, tmp_probe))
        else:
            return tmp_probe


def get_dict(file, d, name: str):
    node = file.get("biological_model")[name]
    if d != {}:
        for key in node:
            d[key] = np.concatenate((d[key], np.array(node[key])))
    else:
        for key in node:
            d[key] = np.array(node[key])
    return d


def merge_hdf5_files(files):
    number_particle = None
    probes = None

    n_t_b = 0
    tgt_file = {}
    with h5py.File(files[0], "r") as f:
        model_group = f.get("biological_model")
        n_t_b = list(model_group.keys()).__len__()
        n_step_key = list(model_group.keys())
        model_key = list(model_group[n_step_key[0]].keys())
        for attr_name, attr_value in f.attrs.items():
            tgt_file[attr_name] = attr_value
    biological_model_data = [{} for _ in range(n_t_b)]
    for _, file_path in enumerate(files):
        with h5py.File(file_path, "r") as f:
            number_particle = get_np(number_particle, f)

            probes = probe(f, probes)

            model_group = f.get("biological_model")

            
            biological_model_data[0] = get_dict(f, biological_model_data[0], "initial")
            for i_ds in range(1, len(n_step_key) - 1):
                biological_model_data[i_ds] = get_dict(f, biological_model_data[i_ds], f"{i_ds}")
            biological_model_data[-1] = get_dict(f, biological_model_data[-1], "final")

  

    with h5py.File("tmp_merge.h5", "w") as f_out:
        f_out.create_dataset("/records/number_particle", data=number_particle,compression="gzip", compression_opts=9)
        if probes is not None:
            f_out.create_dataset("probes", data=probes,compression="gzip", compression_opts=9)

        model_group_out = f_out.create_group("biological_model")

        def mk_ts(ds_name, i):
            ds = model_group_out.create_group(ds_name)
            for key in model_key:
                if key != "spatial":
                    ds.create_dataset(key, data=biological_model_data[i][key],compression="gzip", compression_opts=9)

        mk_ts("initial", 0)

        for i in range(1, len(biological_model_data) - 1):
            mk_ts(f"{i}", i)

        mk_ts("final", len(biological_model_data) - 1)
        for attr_name, attr_value in tgt_file.items():
            f_out.attrs[attr_name] = attr_value


if __name__ == "__main__":
    files = [f"../results/tmp_merge/tmp_merge_partial_{i}.h5" for i in range(6)]
    merge_hdf5_files(files)
