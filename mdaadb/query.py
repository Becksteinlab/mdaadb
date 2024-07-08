from functools import wraps
from collections import UserDict
from dataclasses import dataclass
from typing import List
try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

import pandas as pd


def register_command(dependency=None):
    """"""

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

            return func(self)

        return wrapper

    return decorator


@dataclass
class Command:
    """SQLite command"""
    keyword: str
    fields: List[str]
    dependency: str | None

    def to_sql(self):
        return f"{self.keyword} {','.join(self.fields)}"


class Query:
    """

    Parameters
    ----------
    db : Database or Path or str
        ...

    """

    def __init__(self, db):
        self.db = db
        self.registered_commands: List[Command] = []

    def __repr__(self):
        return self.to_sql()

    def to_sql(self) -> str:
        """Converts the query from it's internal `Query` representation
        into a string that is the exact SQLite query being executed.

        Returns
        -------
        sql_query : str
            ...

        Raises
        ------
        ValueError
            If `self.registered_commands` is empty.
        ValueError
            If a command is missing a dependency command.

        """

        cmds = self.registered_commands

        if not cmds:
            raise ValueError

        keywords = [cmd.keyword for cmd in cmds]
        for cmd in cmds:
            dep = cmd.dependency
            if dep is not None:
                if dep not in keywords:
                    raise ValueError("dependency not found in command")


        sql_query = "\n".join([cmd.to_sql() for cmd in cmds]) + ";"

        return sql_query

    def to_df(self) -> pd.DataFrame:
        """Execute the current Query and return result as pandas DataFrame."""
        if self.db is None:
            raise ValueError("No database to execute query on.")
        return pd.read_sql(self.to_sql(), self.db.connection)

    def execute(self):
        """Execute the current Query a single time.

        Returns
        -------
        `sqlite3.Cursor`
            Database cursor iterator that iterators over results of query

        Raises
        ------
        ValueError
            If `self.db is None`

        """
        if self.db is None:
            raise ValueError("No database to execute query on.")
        return self.db.cursor.execute(self.to_sql())

    def executemany(self, rows):
        """Execute the current parameterized Query for every item in `rows`.

        Parameters
        ----------
        rows : iterable of Row

        Returns
        -------
        `sqlite3.Cursor`
            Database cursor iterator that iterates over results of query

        Raises
        ------
        ValueError
            If `self.db is None`

        """
        if self.db is None:
            raise ValueError("No database to execute query on.")
        return self.db.cursor.executemany(self.to_sql(), rows)

    def commit(self):
        """"""
        self.db.connection.commit()

    @register_command()
    def CREATE_TABLE(self, name, schema) -> Self:
        """Add CREATE TABLE command to the current Query.

        Parameters
        ----------
        name : str
            ...
        schema : str
            ...

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
    def FROM(self, *tables) -> Self:
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
    def INNER_JOIN(self, table) -> Self:
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

    @register_command(dependency="INSERT")
    def INTO(self, table) -> Self:
        """Add (INSERT) INTO command to the current Query.

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

    @register_command(dependency="INTO")
    def INSERT(self) -> Self:
        """Add INSERT (INTO) command to the current Query.

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command()
    def LIMIT(self, limit: int) -> Self:
        """Add LIMIT command to the current Query.

        Parameters
        ----------
        limit : int
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="INNER JOIN")
    def ON(self, condition) -> Self:
        """Add ON command to the current Query.

        Parameters
        ----------
        condition
            ...

        Returns
        -------
        self : query
            ...

        """
        return self

    @register_command()
    def ORDER_BY(self, *columns) -> Self:
        """Add ORDER BY command to the current Query.

        Parameters
        ----------
        *columns
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command()
    def PRAGMA(self, *fields) -> Self:
        """Add PRAGMA command to the current Query.

        Parameters
        ----------
        *fields
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="FROM")
    def SELECT(self, *columns) -> Self:
        """Add SELECT command to the current Query.

        Parameters
        ----------
        columns : iterable of str
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self

    @register_command(dependency="INSERT")
    def VALUES(self, values) -> Self:
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
    def WHERE(self, condition) -> Self:
        """Add WHERE command to the current Query.

        Parameters
        ----------
        conditions
            ...

        Returns
        -------
        self : Query
            ...

        """
        return self
