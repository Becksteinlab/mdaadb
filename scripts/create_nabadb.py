import pathlib

from napalib.system.traj import trajectories
from mdaadb import Database, DBAnalysisRunner


def main():

    def get_topology(traj):
        return traj.topology

    def get_trajectory(traj):
        return traj.trajectory

    def get_repeat(traj):
        return int(traj.name().split("_")[2])

    def get_name(traj):
        repeat = get_repeat(traj)
        name = traj.name().split("_")
        name.remove(f"{repeat}")
        name = "_".join(name)
        return name

    def get_conformation(traj):
        if traj.is_inward:
            return "inward"
        if traj.is_outward:
            return "outward"
        if "_occ_" in traj.name():
            return "occluded"

    def get_temperature(traj):
        if traj.is_310:
            return 310
        if traj.is_358:
            return 358

    def get_protonation(traj):
        protonation = []
        if traj.has_s1:
            protonation.append("s1")
        if traj.has_s2:
            protonation.append("s2")
        if traj.has_s4:
            protonation.append("s4")

        return ",".join(protonation)


    def get_row(idx, traj):
        row = (
            idx,
            get_topology(traj),
            get_trajectory(traj),
            get_name(traj),
            get_repeat(traj),
            get_conformation(traj),
            get_temperature(traj),
            get_protonation(traj),
        )
        return row

    sims_schema = (
        "Simulations (simID INT PRIMARY KEY, topology TEXT, trajectory TEXT, name TEXT, repeat INT, conformation TEXT, temperature INT, protonation TEXT)"
    )
    obs_schema = (
        "Observables (name TEXT, progenitor TEXT)"
    )

    rows = [get_row(idx, traj) for (idx, traj) in enumerate(trajectories)]
    
    dbfile = pathlib.Path("~/projects/napadb/napa.sqlite")
    with Database(dbfile) as db:
        db.create_table(sims_schema)
        db.create_table(obs_schema)

        db.get_table("Simulations").insert_array(rows)


if __name__ == "__main__":
    main()
