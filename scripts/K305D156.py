from collections import namedtuple
import pathlib

import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array

from mdaadb.analysis import DBAnalysisRunner
from mdaadb import Database


class K305_D156(AnalysisBase):

    name = "K305D156"
    results_key = "distance"
    Row = namedtuple("Row", ["simID", "frame", "time", "dA", "dB"])
    schema = f"{name} (simID INT, frame INT, time REAL, dA REAL, dB REAL)"
    _path = str(pathlib.Path(__file__).resolve())

    def __init__(self, universe):
        self._simID = None

        super().__init__(universe.trajectory)
        self.u = universe
        self.n_frames = len(self.u.trajectory)

        self.OD_A = self.u.select_atoms(
            "resid 156 and name OD1 OD2 and segid A"
        )
        self.OD_B = self.u.select_atoms(
            "resid 156 and name OD1 OD2 and segid B"
        )
        self.NZ_A = self.u.select_atoms(
            "resid 305 and name NZ and segid A"
        )
        self.NZ_B = self.u.select_atoms(
            "resid 305 and name NZ and segid B"
        )

    def _prepare(self):
        self.results[self.results_key] = []

    def _single_frame(self):
        simID = self._simID
        ts = self.u.trajectory.ts
        frame = ts.frame
        time = ts.time

        d_A = np.min(
            distance_array(
                self.OD_A.positions,
                self.NZ_A.positions,
                box=self.u.dimensions
            )
        )
        d_B = np.min(
            distance_array(
                self.OD_B.positions,
                self.NZ_B.positions,
                box=self.u.dimensions
            )
        )

        row = self.Row(simID, frame, time, d_A, d_B)

        self.results[self.results_key].append(row)


def main():

    napa_dbfile = pathlib.Path("~/projects/napadb/napa.sqlite")
    napa = Database(napa_dbfile)

    analysis_runner = DBAnalysisRunner(napa, K305_D156)
    with analysis_runner as analysis:
        simulations = analysis.db.get_table("Simulations")
        simids = simulations._get_rowids()

        for simID in simids:
            analysis.run_for_simID(simID)
            analysis.save()


if __name__ == "__main__":
    main()