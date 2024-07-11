from typing import List, NamedTuple
from pathlib import Path

from napalib.system.universe import NapAUniverse
import MDAnalysis as mda

from . import Database


def get_NapA_universe_by_simID(db: Database, simID: int) -> NapAUniverse:
    row = db.get_table("Simulations").get_row(simID)
    topology = row.topology
    trajectory = row.trajectory
    u = NapAUniverse(topology)
    u.load_new(trajectory)

    return u


def get_universe_by_simID(db: Database, simID: int) -> mda.Universe:
    row = db.get_table("Simulations").get_row(simID)
    topology = row.topology
    trajectory = row.trajectory

    return mda.Universe(topology, trajectory)


class DBAnalysisRunner:

    def __init__(self, db: Database, Analysis):

        self.db = db
        self.Analysis = Analysis
        self.name = self.Analysis.name
        self._analysis = None

        try:
            self.observables = self.db.get_table("Observables")
        except ValueError:
            self.observables = self.db.create_table(
                "Observables (name TEXT, progenitor TEXT)"
            )
        finally:
            self.observables.insert_array([
                (self.name, self.Analysis._path),
            ])

    def __enter__(self):
        self.db.open()
        return self

    def __exit__(self, *args):
        self.db.close()

    @property
    def results(self):
        if self._analysis is not None:
            return self._analysis.results

    def run_for_simID(self, simID: int, **kwargs) -> None:
        """"""
        universe = get_NapA_universe_by_simID(self.db, simID)
        self._analysis = self.Analysis(universe, **kwargs)
        self._analysis._simID = simID
        self._analysis.run()

    def save(self) -> None:
        if not self.results:
            raise ValueError("no results")

        if self.name not in self.db._get_table_names():
            self.db.create_table(self.Analysis.schema)

        rows = self.results[self.Analysis.results_key]

        table = self.db.get_table(self.name)
        table.insert_array(rows)

