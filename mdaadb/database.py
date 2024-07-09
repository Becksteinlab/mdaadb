from __future__ import annotations

import sqlite3
import pathlib
from dataclasses import dataclass
from collections import namedtuple, UserDict
from typing import Optional, List, NamedTuple, Iterator, Any
from numpy.typing import ArrayLike
import pandas as pd

from . import query


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
    dbfile : path-like
        Path to the database file.
        Use ':memory:' for a temporary in-memory database.

    """

    def __init__(self, dbfile: pathlib.Path | str):
        self.dbfile = pathlib.Path(dbfile) if isinstance(dbfile, str) else dbfile
        self.connection = sqlite3.connect(self.dbfile)
        self.connection.row_factory = _namedtuple_factory
        self.cursor = self.connection.cursor()

    def __iter__(self):
        return iter(self.tables)

    def __contains__(self, table: Table) -> bool:
        return table.db == self

    @property
    def schema(self) -> pd.DataFrame:
        """Get the database schema as a `pandas.DataFrame`."""
        return self._schema.to_df()

    @property
    def _schema(self) -> Table:
        return Table("sqlite_schema", self)

    @property
    def table_list(self) -> pd.DataFrame:
        """Get a list of tables as a `pandas.DataFrame`."""
        return (
            self._table_list
            .SELECT("*")
            .WHERE("schema='main'")
            .to_df()
        )

    @property
    def _table_list(self) -> Table:
        return Table("pragma_table_list", self)

    @property
    def tables(self):
        """...

        Returns
        -------

        """
        return (Table(name, self) for name in self._get_table_names())

    def _get_table_names(self) -> List[str]:
        """Helper function to get table names in this database."""
        return list(self.table_list["name"])

    def get_table(self, tbl_name) -> Table:
        """...

        Parameters
        ----------
        tbl_name : str
            Table name

        Returns
        -------
        Table

        """
        if tbl_name not in self._get_table_names():
            raise ValueError("invalid table")
        return Table(tbl_name, self)

    def create_table(self, tbl_schema: str, STRICT: bool = True) -> Table:
        """Create a table inside this database.

        Parameters
        ----------
        tbl_schema : str
            Table schema that defines the table and column names.
        STRICT : bool

        Note
        ----
        `self.connection` is used a context manager so that
        `self.connection.commit()` is automatically called.

        """
        if STRICT:
            tbl_schema += " STRICT"

        with self.connection:
            self.CREATE_TABLE(tbl_schema).execute()

        tbl_name = tbl_schema.split()[0]

        return self.get_table(tbl_name)

    def insert_row_into_table(self, table: Table | str, row: tuple) -> None:
        """...

        Parameters
        ----------
        table : Table or str
            ...
        row : tuple
            ...

        """
        if isinstance(table, str):
            table = Table(table, self)
        table.insert_row(row)

    def insert_array_into_table(self, table: Table | str, array: ArrayLike) -> None:
        """...

        Parameters
        ----------
        table : Table or str
            ...
        array : ArrayLike
            ...

        """
        if isinstance(table, str):
            table = Table(table, self)
        table.insert_array(array)

    def CREATE_TABLE(self, tbl_schema: str) -> query.Query:
        """SQLite CREATE_TABLE Statment entry point for `query.Query`.

        Parameters
        ----------
        tbl_schema : str
            ...

        Returns
        -------
        query.Query
            ...

        """
        return query.Query(db=self).CREATE_TABLE(tbl_schema)

    def DELETE(self) -> query.Query:
        """SQLite DELETE Statment entry point for `query.Query`.

        Returns
        -------
        query.Query
            ...

        """
        return query.Query(db=self).DELETE()

    def INSERT_INTO(self, table: str) -> query.Query:
        """SQLite INSERT INTO Statment entry point for `query.Query`.

        Parameters
        ----------
        table : str
            ...

        Returns
        -------
        query.Query
            ...

        """
        return query.Query(db=self).INSERT_INTO(table)

    def PRAGMA(self, schema_name: str) -> query.Query:
        """SQLite PRAGMA Statment entry point for `query.Query`.

        Parameters
        ----------
        schema_name : str
            ...

        Returns
        -------
        query.Query
            ...

        """
        return query.Query(db=self).PRAGMA(schema_name)

    def SELECT(self, *columns: str) -> query.Query:
        """SQLite SELECT Statment entry point for `query.Query`.

        Parameters
        ----------
        *columns : str
            ...

        Returns
        -------
        query.Query
            ...

        """
        return query.Query(db=self).SELECT(*columns)

    def SELECT_COUNT(self, column: str) -> query.Query:
        """SQLite SELECT COUNT Statment entry point for `query.Query`.

        Parameters
        ----------
        column: str
            ...

        Returns
        -------
        query.Query
            ...

        """
        return query.Query(db=self).SELECT_COUNT(column)


@dataclass
class Table:
    name: str
    db: Database

    @property
    def info(self) -> pd.DataFrame:
        """Get the PRAGMA info table as a `pandas.DataFrame`."""
        return (
            self.db
            .PRAGMA(f"table_info('{self.name}')")
            .to_df()
        )

    @property
    def _info(self) -> Table:
        return Table(f"pragma_table_info('{self.name}')", self.db)

    @property
    def schema(self) -> str:
        """Get the schema/structure of this table."""
        _schema = (
            self.db._schema
            .SELECT("sql")
            .WHERE(f"name='{self.name}'")
            .execute()
            .fetchone()
            .sql
        )
        return _schema.replace("CREATE TABLE ", "")

    @property
    def n_rows(self) -> int:
        """Total number of rows in this table."""
        return len(self._get_rowids())

    @property
    def n_cols(self) -> int:
        """Total number of columns in this table."""
        return (
            self.db.
            _table_list
            .SELECT("ncol")
            .WHERE(f"name='{self.name}'")
            .execute()
            .fetchone()
            .ncol
        )

    @property
    def rows(self) -> Iterator[Row]:
        """An iterable of Rows contained in this table."""
        return (Row(id, self) for id in self._get_rowids())

    @property
    def columns(self) -> Iterator[Column]:
        """An iterable of Columns contained in this table."""
        return (Column(name, self) for name in self._get_column_names())

    def get_row(self, id: int) -> Row:
        """

        Parameters
        ----------
        id : int
            INT PRIMARY KEY rowid

        Returns
        -------
        Row

        Raises
        ------
        ValueError

        """
        if id not in self._get_rowids():
            raise ValueError(f"rowid {id} not in rowids")
        return Row(id, self)

    def get_column(self, name: str) -> Column:
        """

        Parameters
        ----------
        name : str
            Name of the column

        Returns
        -------
        Column

        Raises
        ------
        ValueError

        """
        if name not in self._get_column_names():
            raise ValueError("name not in column names")
        return Column(name, self)

    def _get_rowids(self) -> List[int]:
        """Helper function that returns a list of integer indices
        that correspond to the primary key of this table."""
        rowids = (
            self.SELECT(self.pk)
            .execute()
            .fetchall()
        )
        return [id[0] for id in rowids]

    def _get_column_names(self) -> List[str]:
        return list(self.info["name"])

    @property
    def pk(self) -> str:
        """Get the primary key of this table."""
        result = (
            self._info
            .SELECT("name")
            .WHERE("pk=1")
            .execute()
            .fetchone()
        )

        if result is None:
            # if the table has no INT PRIMARY KEY, 'rowid' is the default
            return "rowid"

        return result.name

    def insert_row(self, row: tuple):
        """

        Parameters
        ----------
        row : tuple

        """
        with self.db.connection:
            (
                self.INSERT_INTO(self.name)
                .VALUES(str(row))
                .execute()
            )

    def insert_array(self, array: ArrayLike) -> None:
        """

        Parameters
        ----------
        array : ArrayLike

        """
        n_cols = self.n_cols
        values = "(" + ", ".join(["?" for _ in range(n_cols)]) + ")"
        with self.db.connection:
            (
                self.INSERT_INTO(self.name)
                .VALUES(values)
                .executemany(array)
            )

    def DELETE(self) -> query.Query:
        """

        Returns
        -------
        query.Query

        """
        return self.db.DELETE().FROM(self.name)

    def INSERT_INTO(self, table: str) -> query.Query:
        """

        Parameters
        ----------
        table : str

        Returns
        -------
        query.Query

        """
        return self.db.INSERT_INTO(self.name)

    def SELECT(self, *columns: str) -> query.Query:
        """

        Parameters
        ----------
        *columns : str

        Returns
        -------
        query.Query

        """
        return self.db.SELECT(*columns).FROM(self.name)

    def SELECT_COUNT(self, column: str) -> query.Query:
        """

        Parameters
        ----------
        column : str

        Returns
        -------
        query.Query

        """
        return self.db.SELECT_COUNT(column).FROM(self.name)

    def to_sql(self) -> str:
        """Returns the sqlite query that generates this table."""
        return self.SELECT("*").to_sql()

    def to_df(self) -> pd.DataFrame:
        """Returns a pandas DataFrame of this table."""
        return pd.read_sql(self.to_sql(), self.db.connection)

    def __len__(self) -> int:
        return self.n_rows

    def __iter__(self) -> Iterator[Row]:
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
            .WHERE(f"{self.table.pk}={self.id}")
            .execute()
            .fetchone()
        )

    def to_sql(self) -> str:
        """Returns the sqlite query that generates this row."""
        return (
            self.table
            .SELECT("*")
            .WHERE(f"{self.table.pk}={self.id}")
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
        return (
            self.table._info
            .SELECT("type")
            .WHERE(f"name='{self.name}'")
            .execute()
            .fetchone()
            .type
        )

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

