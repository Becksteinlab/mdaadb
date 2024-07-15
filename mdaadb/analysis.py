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
        Analysis :
        dbfile :
        hooks :

        """

        self.Analysis = Analysis
        self.db = Database(dbfile)

        try:
            self.analysis_name = self.Analysis.name
        except AttributeError:
            self.analysis_name = self.Analysis.__name__
        try:
            self.analysis_notes = self.Analysis.notes
        except AttributeError:
            self.analysis_notes = None
        self.analysis_path = inspect.getfile(self.Analysis)

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

        if self.analysis_name not in self.obsv.get_column("obsName").data:
            self.obsv.insert_row(
                (self.analysis_name, self.analysis_notes, self.analysis_path),
                columns=["obsName, notes, creator"],
            )

        self._analysis = None

    @property
    def results(self) -> Results:
        """Analysis results."""

        if self._analysis is None:
            raise ValueError("Must call run() for results to exist.")
        return self._analysis.results

    def _get_universe(self, simID: int, get_universe: Optional[Callable]):

        if get_universe is not None:
        #if self.hooks["get_universe"]:
            return get_universe(self.db, simID)

        row = self.db.get_table("Simulations").get_row(simID)
        return mda.Universe(row.topology, row.trajecory)

    def run(
        self,
        simID: int,
        get_universe: Optional[Callable] = None,
        **kwargs: dict,
        ) -> None:
        """
        Parameters
        ----------
        simID : int
        get_universe : Callable[Database, int]
        **kwargs : dict
            additional keyword arguments to be passed to the Analysis class

        """
        u = self._get_universe(simID, get_universe)

        self._analysis = self.Analysis(u, **kwargs)
        self._analysis._simID = simID
        self._analysis.run()

    def save(self) -> None:
        """Save the results of the analysis to the database."""

        if not self.results:
            raise ValueError("no results")

        analysis_table = Table(self.analysis_name, self.db)

        if analysis_table not in self.db:
            self.db.create_table(self.Analysis.schema)
        else:
            simID = self._analysis._simID
            if simID in analysis_table._get_rowids():
                raise ValueError(
                    f"{self.analysis_name} table already has data for simID {simID}"
                )

        rows = self.results[self.Analysis.results_key]

        self.db.get_table(self.analysis_name).insert_rows(rows)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.db.close()
