from __future__ import annotations

import sqlite3
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from collections import namedtuple, UserDict
from typing import List, Iterable, Any

from .query import Query


def _namedtuple_factory(cursor, row):
    fields = [col[0] for col in cursor.description]
    Row = namedtuple("Row", fields)
    return Row(*row)


class Tables(UserDict):
    def __init__(self, db, *args, **kwargs):
        self.db = db


class Database:
    """

    Parameters
    ----------
    dbfile : pathlib.Path or str
        ...

    """

    def __init__(self, dbfile: Path | str):
        if isinstance(dbfile, str):
            self.dbfile = Path(dbfile)
        self.connection = sqlite3.connect(self.dbfile)
        self.connection.row_factory = _namedtuple_factory
        self.cursor = self.connection.cursor()

    def __iter__(self):
        return iter(self.tables)

    def __contains__(self, table):
        return table.db == self

    def table(self, name):
        """...

        Parameters
        ----------
        name : str
            table name

        Returns
        -------
        Table

        """
        if name not in self._get_table_names():
            raise ValueError("invalid table")
        return Table(name, self)

    @property
    def tables(self):
        """...

        Returns
        -------

        """
        return (Table(name, self) for name in self._get_table_names())

    def _get_table_names(self) -> List[str]:
        """...

        Returns
        -------

        """
        table_names = (
            self.SELECT("name")
            .FROM("sqlite_schema")
            .WHERE("type = 'table'")
            .execute()
            .fetchall()
        )
        return [name[0] for name in table_names]

    @property
    def schema(self) -> pd.DataFrame:
        """Get the database schema as a `pandas.DataFrame`"""
        return Query(db=self).SELECT("*").FROM("sqlite_schema").to_df()

    def insert_row_into_table(self, table, row) -> None:
        """...

        Parameters
        ----------
        table : Table or str
            ...
        row : Row
            ...

        Returns
        -------

        """
        if isinstance(table, str):
            table = Table(table, self)
        table.insert_row(row)

    def insert_array_into_table(self, table, array) -> None:
        """...

        Parameters
        ----------
        table : Table or str
            ...
        array : array_like
            ...

        Returns
        -------

        """
        if isinstance(table, str):
            table = Table(table, self)
        table.insert_array(array)

    def CREATE_TABLE(self, name, schema) -> Query:
        """SQLite CREATE_TABLE Statment entry point for ``Query``.

        Parameters
        ----------
        name : str
            ...
        schema : str
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).CREATE_TABLE(name, schema)

    def DELETE(self) -> Query:
        """SQLite DELETE Statment entry point for ``Query``.

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).DELETE()

    def FROM(self, *tables) -> Query:
        """SQLite FROM Statment entry point for ``Query``.

        Parameters
        ----------
        *tables
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).FROM(*tables)

    def INNER_JOIN(self, table) -> Query:
        """SQLite INNER_JOIN Statment entry point for ``Query``.

        Parameters
        ----------
        table : str
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).INNER_JOIN(table)

    def INSERT(self) -> Query:
        """SQLite INSERT Statment entry point for ``Query``.

        Returns
        -------
        Query

        """
        return Query(db=self).INSERT()

    def INTO(self, table) -> Query:
        """SQLite INTO Statment entry point for ``Query``.

        Parameters
        ----------
        table : Table
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).INTO(table)

    def LIMIT(self, limit: int) -> Query:
        """SQLite LIMIT Statment entry point for ``Query``.

        Parameters
        ----------
        limit : int
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).LIMIT(limit)

    def ON(self, condition) -> Query:
        """SQLite ON Statment entry point for ``Query``.

        Parameters
        ----------
        condition
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).ON(condition)

    def ORDER_BY(self, *fields) -> Query:
        """SQLite ORDER BY Statment entry point for ``Query``.

        Parameters
        ----------
        *fields
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).ORDER_BY(*fields)

    def PRAGMA(self, *fields) -> Query:
        """SQLite PRAGMA Statment entry point for ``Query``.

        Parameters
        ----------
        *fields
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).PRAGMA(*fields)

    def SELECT(self, *columns) -> Query:
        """SQLite SELECT Statment entry point for ``Query``.

        Parameters
        ----------
        *columns
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).SELECT(*columns)

    def VALUES(self, values) -> Query:
        """SQLite VALUES Statment entry point for ``Query``.

        Parameters
        ----------
        values
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).VALUES(*values)

    def WHERE(self, condition) -> Query:
        """SQLite WHERE Statment entry point for ``Query``.

        Parameters
        ----------
        condition
            ...

        Returns
        -------
        Query
            ...

        """
        return Query(db=self).WHERE(condition)


@dataclass
class Table:
    name: str
    db: Database

    @property
    def schema(self) -> str:
        """"""
        return (
            self.db
            .SELECT("sql")
            .FROM("sqlite_schema")
            .WHERE(f"name = '{self.name}'")
            .execute()
            .fetchone()[0]
        )

    @property
    def info(self) -> pd.DataFrame:
        """

        Returns
        -------
        pd.DataFrame

        """
        return self.db.PRAGMA(f"table_info('{self.name}')").to_df()

    def row(self, id:int) -> Row:
        """

        Parameters
        ----------
        id : int

        Returns
        -------
        Row

        Raises
        ------

        """
        if id not in self._get_row_ids():
            raise ValueError("id not in valid ids")
        return Row(id, self)

    def column(self, name:str) -> Column:
        """

        Parameters
        ----------
        name : str

        Returns
        -------
        Column

        Raises
        ------

        """
        if name not in self._get_column_names():
            raise ValueError("name not in column names")
        return Column(name, self)

    @property
    def n_rows(self) -> int:
        """Total number of rows in this table."""
        return len(self._get_row_ids())

    @property
    def n_cols(self) -> int:
        """Total number of columns in this table."""
        return len(self._get_column_names())

    @property
    def rows(self) -> Iterable[Row]:
        """An iterable of Rows contained in this table."""
        return (Row(id, self) for id in self._get_row_ids())

    @property
    def columns(self) -> Iterable[Column]:
        """An iterable of Columns contained in this table."""
        return (Column(name, self) for name in self._get_column_names())

    def _get_row_ids(self) -> List[int]:
        """Helper function that returns a list of integer indices
        that correspond to the primary key of this table."""
        row_ids = (
            self.SELECT(self.primary_key)
            .execute()
            .fetchall()
        )
        return [id[0] for id in row_ids]

    def _get_column_names(self) -> List[str]:
        """Helper function that returns a list of column names for this table.."""
        column_names = (
            self.db
            .SELECT("name")
            .FROM(f"pragma_table_info('simulations')")
            .execute()
            .fetchall()
        )
        return [name[0] for name in column_names]

    @property
    def primary_key(self) -> str:
        """Returns the primary key of this table as a string."""
        result = (
            self.db
            .SELECT("name")
            .FROM(f"pragma_table_info('{self.name}') as tblinfo")
            .WHERE("tblinfo.pk = 1")
            .execute()
            .fetchall()
        )
        assert len(result[0]) == 1
        return result[0][0]

    def insert_row(self, row) -> None:
        """

        Parameters
        ----------
        row : Row | tuple

        """
        return (
            self.db
            .INSERT()
            .INTO(self.name)
            .VALUES(row)
            .execute()
            .commit()
        )

    def insert_array(self, array) -> None:
        """

        Parameters
        ----------
        array : List[Row]

        """
        n_cols = self.n_cols
        values = "(" + ", ".join(["?" for _ in range(n_cols)]) + ")"
        return (
            self.db
            .INSERT()
            .INTO(self.name)
            .VALUES(values)
            .executemany(array)
            .commit()
        )

    def DELETE(self, *rows) -> Query:
        """

        Parameters
        ----------
        *rows

        Returns
        -------
        Query

        """
        return self.db.DELETE(*rows).FROM(self.name)

    def INSERT(self, *rows) -> Query:
        """

        Parameters
        ----------

        Returns
        -------
        Query

        """
        return self.db.INSERT(*rows).INTO(self.name)

    def SELECT(self, *columns) -> Query:
        """

        Parameters
        ----------

        Returns
        -------
        Query

        """
        return self.db.SELECT(*columns).FROM(self.name)

    def WHERE(self, *conditions) -> Query:
        """

        Parameters
        ----------

        Returns
        -------
        Query

        """
        return self.db.WHERE(*conditions).FROM(self.name)

    def to_sql(self) -> str:
        """Returns the sqlite query that generates this table."""
        return self.SELECT("*").to_sql()

    def to_df(self) -> pd.DataFrame:
        """Returns a pandas DataFrame of this table."""
        return pd.read_sql(self.to_sql(), self.db.connection)

    def __len__(self) -> int:
        return self.n_rows

    def __iter__(self) -> Iterable[Row]:
        return iter(self.rows)


@dataclass
class Row:
    id: int
    table: Table

    @property
    def db(self) -> Database:
        """Returns the database that this row belongs to."""
        return self.table.db

    @property
    def data(self):
        """"""
        return (
            self.table
            .SELECT("*")
            .WHERE(f"{self.table.primary_key}={self.id}")
            .execute()
            .fetchall()
        )

    def to_sql(self) -> str:
        """Returns the sqlite query that generates this row."""
        return (
            self.table
            .SELECT("*")
            .WHERE(f"{self.table.primary_key}={self.id}")
            .to_sql()
        )

    def to_df(self) -> pd.DataFrame:
        """Returns a pandas DataFrame of this row."""
        return pd.read_sql(self.to_sql(), self.db.connection)


@dataclass
class Column:
    name: str
    table: Table

    @property
    def db(self) -> Database:
        """"""
        return self.table.db

    @property
    def type_(self) -> str:
        """"""
        _type = (
            Query(db=self.db)
            .SELECT("type")
            .FROM(f"pragma_table_info('{self.table.name}') as t_info")
            .WHERE(f"t_info.name='{self.name}'")
            .execute()
            .fetchall()
        )
        return _type[0][0]

    @property
    def data(self) -> List[Any]:
        """"""
        _data = self.table.SELECT(self.name).execute().fetchall()
        return [data[0] for data in _data]

    def to_sql(self) -> str:
        """"""
        return (
            self.table
            .SELECT(self.name)
            .to_sql()
        )

    def to_df(self) -> pd.DataFrame:
        """"""
        return pd.read_sql(self.to_sql(), self.db.connection)

