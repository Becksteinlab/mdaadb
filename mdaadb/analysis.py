import inspect
import pathlib
from typing import Optional, Callable

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase, Results

from .database import Database, Table


class DBAnalysisManager:
    """Class that connects database IO with analysis running.

    Parameters
    ----------
    Analysis :
    dbfile : path-like
    hooks : dict

    """

    def __init__(
        self,
        Analysis, # : class that inherits from AnalysisBase
        dbfile: pathlib.Path | str,
        hooks: Optional[dict] = None,
    ) -> None:
        self.Analysis = Analysis
        if isinstance(dbfile, Database):
            self.db = dbfile
        else:
            self.db = Database(dbfile)
        self._analysis = None

        self.hooks = {
            "pre_run": None,
            "post_run": None,
            "pre_save": None,
            "post_save": None,
            "get_universe": None,
        }
        if hooks is not None:
            self.hooks.update(hooks)

        try:
            self._name = self.Analysis.name
        except AttributeError:
            self._name = self.Analysis.__name__
        try:
            self._desc = self.Analysis.description
        except AttributeError:
            self._desc = "none"
        try:
            self._creator = inspect.getfile(self.Analysis)
        except OSError:
            self._creator = "unknown"

        try:
            self.observables = self.db.get_table("Observables")
        except ValueError:
            self.observables = self.db.create_table(
                """
                Observables (
                    name TEXT PRIMARY KEY,
                    description TEXT,
                    creator TEXT
                )
                """,
                STRICT=False,
            )

        if self._name not in self.observables.get_column("name").data:
            # print(self._name, self._desc, self._path)
            self.observables.insert_row(
                row=(self._name, self._desc, self._creator),
                columns=["name", "description", "creator"],
            )

    def _get_universe(self, simID: int) -> mda.Universe:

        if self.hooks["get_universe"] is not None:
            universe = self.hooks["get_universe"](self.db, simID)

            return universe

        row = self.db.get_table("Simulations").get_row(simID)
        universe = mda.Universe(row.Topology, row.Trajectory)

        return universe

    def run(self, simID: int, **kwargs: dict) -> None:
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
            self.hooks["pre_run"](self.db, simID)

        self._analysis.run()

        if self.hooks["post_run"] is not None:
            self.hooks["post_run"](self.db, simID)

    def save(self) -> None:
        """Save the results of the analysis to the database."""

        assert self._analysis is not None

        if not self.results:
            raise ValueError("no results")

        try:
            analysis_table = self.db.get_table(self._name)
        except ValueError:
            analysis_table = self.db.create_table(self.Analysis.schema)
        # else:
        #     assert analysis_table.schema == self.Analysis.schema

        simID = self._analysis._simID
        if simID in [row[0] for row in analysis_table.SELECT_DISTINCT("simID").execute().fetchall()]:
            raise ValueError(
                f"'{self._name}' table already has data for simID {simID}"
            )

        if self.hooks["pre_save"] is not None:
            self.hooks["pre_save"](self.db, simID)

        rows = self.results[self.Analysis.results_key]
        analysis_table.insert_rows(rows)

        if self.hooks["post_save"] is not None:
            self.hooks["post_save"](self.db, simID)

    @property
    def results(self) -> Results:
        if self._analysis is None:
            raise ValueError("Must call run() for results to exist.")
        return self._analysis.results

    def __enter__(self):
        self.db.open()
        return self

    def __exit__(self, *args):
        self.db.close()
