import inspect
from typing import List, NamedTuple, Optional, Callable

import MDAnalysis as mda
from MDAnalysis.analysis.base import Results

from .database import Database, Table


class DBAnalysisManager:
    """Class that connects database IO with analysis running.

    """

    def __init__(self, Analysis, dbfile, hooks=None):
        """

        Parameters
        ----------
        Analysis : mda.analysis.base.AnalysisBase
        dbfile : path-like
        hooks : dict

        """

        self.Analysis = Analysis
        self.db = Database(dbfile)
        self._analysis = None

        self.hooks = {
            "pre_run": None,
            "post_run": None,
            "get_universe": None,
            "post_save": None,
        }
        if hooks is not None:
            self.hooks.update(hooks)

        try:
            self._name = self.Analysis.name
        except AttributeError:
            self._name = self.Analysis.__name__
        try:
            self._notes = self.Analysis.notes
        except AttributeError:
            self._notes = None
        self._path = inspect.getfile(self.Analysis)

        try:
            self.obsv = self.db.get_table("Observables")
        except ValueError:
            self.obsv = self.db.create_table(
                """
                Observables (
                    obsName TEXT,
                    notes TEXT,
                    creator TEXT,
                    timestamp DATETIME DEFAULT (strftime('%m-%d-%Y %H:%M', 'now', 'localtime'))
                )
                """,
                STRICT=False
            )

        if self._name not in self.obsv.get_column("obsName").data:
            self.obsv.insert_row(
                (self._name, self._notes, self._path),
                columns=["obsName, notes, creator"],
            )

    @property
    def results(self) -> Results:
        """Analysis results"""

        if self._analysis is None:
            raise ValueError("Must call run() for results to exist.")
        return self._analysis.results

    def _get_universe(self, simID: int):

        if self.hooks["get_universe"] is not None:
        #if self.hooks["get_universe"]:
            return self.hooks["get_universe"](self.db, simID)

        row = self.db.get_table("Simulations").get_row(simID)
        return mda.Universe(row.topology, row.trajectory)

    def run(self, simID: int,  **kwargs: dict) -> None:
        """Run the analysis for a simulation given by `simID`.

        Parameters
        ----------
        simID : int
        **kwargs : dict
            additional keyword arguments to be passed to the Analysis class

        """

        univserse = self._get_universe(simID)

        self._analysis = self.Analysis(univserse, **kwargs)
        self._analysis._simID = simID

        if self.hooks["pre_run"] is not None:
            self.hooks["pre_run"](simID, self.db)

        self._analysis.run()

        if self.hooks["post_run"] is not None:
            self.hooks["post_run"](simID, self.db)

    def save(self) -> None:
        """Save the results of the analysis to the database."""

        assert self._analysis is not None

        if not self.results:
            raise ValueError("no results")

        try:
            analysis_table = self.db.get_table(self._name)
        except ValueError:
            analysis_table = self.db.create_table(self.Analysis.schema)
        else:
            assert analysis_table.schema == self.Analysis.schema

        simID = self._analysis._simID
        if simID in analysis_table._get_rowids():
            raise ValueError(
                f"{self._name} table already has data for simID {simID}"
            )

        rows = self.results[self.Analysis.results_key]
        analysis_table.insert_rows(rows)

        if self.hooks["post_save"] is not None:
            self.hooks["post_save"](simID, self.db)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.db.close()
