from functools import wraps
from collections import UserDict
from dataclasses import dataclass
from typing import Sequence
try:
    from typing import Self
except ImportError:
    from typing_extensions import Self
import pathlib
import sqlite3
import pandas as pd

from . import database


def register_command(dependency=None):
    """SQLite keyword decorator for `Query` methods.

    When a method is called, the keyword is bundled with input fields
    into a `Command` object, then registered into `Query.registered_commands`.

    """
    def decorator(func):

        @wraps(func)
        def wrapper(self, *args):
            keyword = wrapper.__name__.replace("_", " ")
            if args:
                fields = [arg for arg in args]
            else:
                fields = [""]
            self.registered_commands.append(
                Command(keyword, fields, dependency)
            )
            return func(self, *args)

        return wrapper

    return decorator


@dataclass
class Command:
    """SQLite command"""
    keyword: str
    fields: list[str]
    dependency: str | None

    def to_sql(self):
        return f"{self.keyword} {','.join(self.fields)}"


class Query:
    """SQLite query object.

    Implements the "builder" design pattern to build
    a full query statement from individual SQLite commands.

    Parameters
    ----------
    db : `Database` or path-like
        Database to execute queries on. If path-like, will
        call `Database(db)`. Use ":memory:" to execute queries
        on a temporary in memory database.

    """

    def __init__(self, db: database.Database | pathlib.Path | str) -> None:
        if isinstance(db, pathlib.Path) or isinstance(db, str):
            self.db = database.Database(db)
        else:
            assert isinstance(db, database.Database)
            self.db = db
        self.registered_commands: list[Command] = []

    def __repr__(self):
        return self.to_sql()

    def to_sql(self) -> str:
        """Convert the query from it's internal `Query` representation
        into a string that is the exact SQLite query/statement being executed.

        Returns
        -------
        statement : str
            The complete SQLite statement to be executed
            with proper SQL syntax.

        Raises
        ------
        ValueError
            If `self.registered_commands` is empty.
        ValueError
            If a command is missing a dependency command.

        """
        cmds = self.registered_commands

        if not cmds:
            raise ValueError("Query must consist of at least one command.")

        keywords = [cmd.keyword for cmd in cmds]
        for cmd in cmds:
            dep = cmd.dependency
            if dep is not None:
                if dep not in keywords:
                    raise ValueError("dependency not found in command")


        statement = "\n".join([cmd.to_sql() for cmd in cmds]) + ";"

        return statement

    def to_df(self) -> pd.DataFrame:
        """Execute the current Query and return result as pandas DataFrame."""
        return pd.read_sql(self.to_sql(), self.db.connection)

    def execute(self) -> sqlite3.Cursor:
        """Execute the current Query a single time.

        Returns
        -------
        `sqlite3.Cursor`
            Database cursor iterator that iterators over results of query

        """
        return self.db.connection.execute(self.to_sql())

    def executemany(self, rows) -> sqlite3.Cursor:
        """Execute the current parameterized Query for every item in `rows`.

        Parameters
        ----------
        rows : sequence

        Returns
        -------
        `sqlite3.Cursor`
            Database cursor iterator that iterates over results of query

        """
        return self.db.connection.executemany(self.to_sql(), rows)

    @register_command(dependency="ALTER TABLE")
    def ADD_COLUMN(self, column_def: str) -> Self:
        """Add ADD COLUMN command to the current Query.

        Parameters
        ----------
        column_def : str

        Returns
        -------
        self : `Query`
            ...

        """
        return self

    @register_command()
    def ALTER_TABLE(self, table_name: str) -> Self:
        """Add ALTER TABLE command to the current Query.

        Parameters
        ----------
        table_name : str

        Returns
        -------
        self : `Query`

        """
        return self

    @register_command()
    def CREATE_TABLE(self, table_schema: str) -> Self:
        """Add CREATE TABLE command to the current Query.

        Parameters
        ----------
        table_schema : str
            The table schema/structure.

        Returns
        -------
        self : `Query`
            ...

        """
        return self

    @register_command(dependency="FROM")
    def DELETE(self) -> Self:
        """Add DELETE command to the current Query.

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command()
    def FROM(self, *tables: str) -> Self:
        """Add FROM command to the current Query.

        Parameters
        ----------
        *tables : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="FROM")
    def INNER_JOIN(self, table: str) -> Self:
        """Add INNER JOIN command to the current Query.

        Parameters
        ----------
        table : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command()
    def INSERT_INTO(self, table: str) -> Self:
        """Add INSERT INTO command to the current Query.

        Parameters
        ----------
        table : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command()
    def LIMIT(self, limit: str) -> Self:
        """Add LIMIT command to the current Query.

        Parameters
        ----------
        limit : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="INNER JOIN")
    def ON(self, match: str) -> Self:
        """Add ON command to the current Query.

        Parameters
        ----------
        match : str
            ...

        Returns
        -------
        self : query
            ...

        """
        return self

    @register_command()
    def ORDER_BY(self, *columns: str) -> Self:
        """Add ORDER BY command to the current Query.

        Parameters
        ----------
        *columns : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command()
    def PRAGMA(self, schema_name: str) -> Self:
        """Add PRAGMA command to the current Query.

        Parameters
        ----------
        schema_name
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="FROM")
    def SELECT(self, *columns: str) -> Self:
        """Add SELECT command to the current Query.

        Parameters
        ----------
        *columns : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="FROM")
    def SELECT_COUNT(self, column: str) -> Self:
        """Add SELECT COUNT command to the current Query.

        Parameters
        ----------
        column : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="FROM")
    def SELECT_DISTINCT(self, *columns: str) -> Self:
        """Add SELECT DISTINCT command to the current Query.

        Parameters
        ----------
        *columns : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="UPDATE")
    def SET(self, *column_exprs) -> Self:
        """Add SET command to the current Query.

        Parameters
        ----------
        *column_exprs : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command()
    def UPDATE(self, table_name: str) -> Self:
        """Add UPDATE command to the current Query.

        Parameters
        ----------
        table_name : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="INSERT INTO")
    def VALUES(self, values: str) -> Self:
        """Add VALUES command to the current Query.

        Parameters
        ----------
        values : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command()
    def WHERE(self, condition: str) -> Self:
        """Add WHERE command to the current Query.

        Parameters
        ----------
        condition : str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self
