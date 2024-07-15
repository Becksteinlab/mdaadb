from __future__ import annotations

import sqlite3
import pathlib
from dataclasses import dataclass
from collections import namedtuple, UserDict
from typing import (
    Optional,
    List,
    Tuple,
    NamedTuple,
    Iterator,
    Sequence,
    Any,
)
import pandas as pd

from . import query


def _namedtuple_factory(cursor, row):
    fields = [col[0] for col in cursor.description]
    Row = namedtuple("Row", fields)
    return Row(*row)


class Database:

    def __init__(self, dbfile: pathlib.Path | str):
        """Database constructor.

        Parameters
        ----------
        dbfile : path-like
            Path to the database file.
            Use ':memory:' for a temporary in-memory database.

        """
        if isinstance(dbfile, str):
            self.dbfile = pathlib.Path(dbfile)
        else:
            self.dbfile = dbfile
        self.connection = None
        self.cursor = None
        self.open()

    def open(self):
        """Open a connection to `self.dbfile` if it is not already open."""
        if self.connection is None:
            self.connection = sqlite3.connect(self.dbfile)
            self.connection.row_factory = _namedtuple_factory
            self.cursor = self.connection.cursor()

    def close(self):
        """Close the connection to `self.dbfile` if it is open."""
        if self.connection is not None:
            self.connection.close()
            self.connection = None
            self.cursor = None

    def get_table(self, table_name: str) -> Table:
        """Get a table from this database by name.

        Parameters
        ----------
        table_name : str

        Returns
        -------
        Table

        """
        if table_name not in self._get_table_names():
            raise ValueError(f"'{table_name}' not in database")
        return Table(table_name, self)

    def create_table(
        self,
        schema: Schema | str,
        STRICT: Optional[bool] = True,
        ) -> Table:
        """Create a table inside this database.

        Parameters
        ----------
        schema : Schema or str
            Table schema that defines the table and column names.
        STRICT : bool (default True)
            Set the SQLite 'STRICT' flag for CREATE TABLE command.

        Note
        ----
        `self.connection` is used a context manager so that
        `self.connection.commit()` is automatically called.

        """

        if self.connection is None:
            raise ValueError("Database connection is closed.")

        if isinstance(schema, str):
            schema = Schema.from_str(schema)

        table_name = schema.table_name
        schema_str = schema.to_sql()

        if STRICT:
            schema_str += "\nSTRICT"

        with self.connection:
            self.CREATE_TABLE(schema_str).execute()

        return self.get_table(table_name)

    def insert_row_into_table(
        self,
        table: Table | str,
        row: Tuple,
        columns: Optional[Sequence[str]] = None,
        ) -> None:
        """Insert a single row into a table within this database.

        Parameters
        ----------
        table : Table or str
        row : tuple
        columns : sequence of strings
            Column names to insert rows into.
            Must be defined if not inserting an element into every column.

        """
        if self.connection is None:
            raise ValueError("Database connection is closed.")

        if isinstance(table, str):
            table = Table(table, self)

        if columns:
            column_names = ", ".join(columns)
            tbl_with_cols = f"{table.name} ({column_names})"
            with self.connection:
                (
                    self.INSERT_INTO(tbl_with_cols)
                    .VALUES(str(row))
                    .execute()
                )
            return

        with self.connection:
            (
                self.INSERT_INTO(table.name)
                .VALUES(str(row))
                .execute()
            )

    def insert_rows_into_table(
        self,
        table: Table | str,
        rows: List[Tuple],
        columns: Optional[Sequence[str]] = None,
        ) -> None:
        """Insert multiple rows into a table within this database.

        Parameters
        ----------
        table : Table or str
        rows : list of tuples
            Sequence of rows to be inserted.
        columns : sequence of strings
            Column names to insert rows into.
            Must be defined if not inserting an element into every column.

        """
        if self.connection is None:
            raise ValueError("Database connection is closed.")

        if isinstance(table, str):
            table = Table(table, self)

        if columns:
            place_holders = "(" + ", ".join(len(columns)*"?".split()) + ")"
            column_names = ", ".join(columns)
            tbl_with_cols = f"{table.name} ({column_names})"
            with self.connection:
                (
                    self.INSERT_INTO(tbl_with_cols)
                    .VALUES(place_holders)
                    .executemany(rows)
                )

            return

        place_holders = "(" + ", ".join(table.n_cols*"?".split()) + ")"
        with self.connection:
            (
                self.INSERT_INTO(table.name)
                .VALUES(place_holders)
                .executemany(rows)
            )

    def add_column_to_table(
        self,
        table: Table | str,
        column_def: str,
        ) -> None:
        """Add a column to a table within this database.

        Parameters
        ----------
        table : Table or str
        column_def : str
            Column definition with fields '<column-name> <TYPE>'.
            <TYPE> must be a valid SQLite type.

        """
        if self.connection is None:
            raise ValueError("Database connection is closed.")

        if isinstance(table, str):
            table = Table(table, self)

        with self.connection:
            (
                self.ALTER_TABLE(table.name)
                .ADD_COLUMN(column_def)
                .execute()
            )

    def update_column_in_table(
        self,
        table: Table | str,
        column_name: str,
        values: Sequence,
        ids: Optional[Sequence[int]] = None,
        ) -> None:
        """Update a column of a table within this database.

        Parameters
        ----------
        table : Table or str
        column_name : str
        values : sequence of values
        ids : sequence of int

        """
        if self.connection is None:
            raise ValueError("Database connection is closed.")

        if isinstance(table, str):
            table = Table(table, self)

        if ids is None:
            assert len(values) == len(table)
            # better to use generator or list comprehension here ??
            data = (
                (val, id)
                for val, id in zip(values, table._get_rowids())
            )
        else:
            assert len(values) == len(ids)
            data = (
                (val, id)
                for val, id in zip(values, ids)
            )

        with self.connection:
            (
                self.UPDATE(table.name)
                .SET(f"{column_name} = ?")
                .WHERE(f"{table.pk} = ?")
                .executemany(data)
            )

    def insert_column_into_table(
        self,
        table: Table | str,
        column_def: str,
        values: Sequence,
        ) -> None:
        """Insert a column into a table within this database.

        Parameters
        ----------
        table : Table or str
        column_def : str
        values : sequence of values

        """
        if self.connection is None:
            raise ValueError("Database connection is closed.")

        if isinstance(table, str):
            table = Table(table, self)

        col_name, col_type = column_def.strip().split()
        if col_name in table._get_column_names():
            raise ValueError(
                f"Column '{col_name}' is already in table '{table.name}'."
            )

        self.add_column_to_table(table, column_def)
        self.update_column_in_table(table, col_name, values)

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
        """An iterable of `Table` in this database

        TODO: implement this as a UserDict container?

        """
        return (Table(name, self) for name in self._get_table_names())

    def _get_table_names(self) -> List[str]:
        """Helper function to get table names in this database."""
        return list(self.table_list["name"])

    def ALTER_TABLE(self, table_name: str) -> query.Query:
        """SQLite ALTER TABLE statement entry point via database.

        Parameters
        ----------
        table_name : str

        Returns
        -------
        query.Query

        """
        return query.Query(db=self).ALTER_TABLE(table_name)

    def CREATE_TABLE(self, tbl_schema: str) -> query.Query:
        """SQLite CREATE TABLE statement entry point via database.

        Parameters
        ----------
        tbl_schema : str

        Returns
        -------
        query.Query

        """
        return query.Query(db=self).CREATE_TABLE(tbl_schema)

    def DELETE(self) -> query.Query:
        """SQLite DELETE statement entry point via database.

        Returns
        -------
        query.Query

        """
        return query.Query(db=self).DELETE()

    def INSERT_INTO(self, table: str) -> query.Query:
        """SQLite INSERT INTO statement entry point via database.

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
        """SQLite PRAGMA statement entry point via database.

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
        """SQLite SELECT statement entry point via database.

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
        """SQLite SELECT COUNT statement entry point via database.

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

    def SELECT_DISTINCT(self, *columns: str) -> query.Query:
        """SQLite SELECT DISTINCT statement entry point via database.

        Parameters
        ----------
        *columns: str
            ...

        Returns
        -------
        query.Query
            ...

        """
        return query.Query(db=self).SELECT_DISTINCT(*columns)

    def UPDATE(self, table_name: str) -> query.Query:
        """SQlite UPDATE statement entry point via database.

        Parameters
        ----------
        table_name : str

        Returns
        -------
        query.Query

        """
        return query.Query(db=self).UPDATE(table_name)

    def __contains__(self, table: Table | str) -> bool:
        if isinstance(table, str):
            table = Table(table, self)
        return (table.name in self._get_table_names()) and (table.db == self)

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()

    def __iter__(self):
        return iter(self.tables)


@dataclass
class Table:
    name: str
    db: Database

    def get_row(self, id: int) -> Row:
        """Get a row from this table by integer index.

        Parameters
        ----------
        id : int

        Returns
        -------
        Row

        Raises
        ------
        ValueError

        """
        if id not in self._get_rowids():
            raise ValueError(f"{self.pk} {id} not found in '{self.name}' table")

        return Row(id, self)

    def get_column(self, name: str) -> Column:
        """Get a column from this table by name.

        Parameters
        ----------
        name : str

        Returns
        -------
        Column

        Raises
        ------
        ValueError

        """
        if name not in self._get_column_names():
            raise ValueError(f"{name} not a column in '{self.name}' table")

        return Column(name, self)

    def insert_row(
        self,
        row: Tuple,
        columns: Optional[Sequence[str]] = None,
        ) -> None:
        """Insert a single row into this table.

        Parameters
        ----------
        row : tuple
        columns : sequence of strings
            Column names to insert rows into.
            Must be defined if not inserting an element into every column.

        """
        self.db.insert_row_into_table(self, row, columns=columns)

    def insert_rows(
        self,
        rows: List[Tuple],
        columns: Optional[Sequence[str]] = None,
        ) -> None:
        """Insert multiple rows into this table.

        Parameters
        ----------
        rows : list of tuples
            Sequence of rows to be inserted.
        columns : sequence of strings
            Column names to insert rows into.
            Must be defined if not inserting an element into every column.

        """
        self.db.insert_rows_into_table(self, rows, columns=columns)

    def add_column(
        self,
        column_def: str
        ) -> None:
        """Add a column to this table.

        Parameters
        ----------
        column_def : str
            Column definition with fields '<column-name> <TYPE>'.
            <TYPE> must be a valid SQLite type.

        """
        self.db.add_column_to_table(self, column_def)

    def update_column(
        self,
        column: str,
        values: Sequence,
        ids: Optional[Sequence[int]] = None,
        ) -> None:
        """Update a column in this table.

        Parameters
        ----------
        column : str
        values :
        ids :

        """
        self.db.update_column_in_table(self, column, values, ids=ids)

    def insert_column(
        self,
        column_name: str,
        values: Sequence,
        ) -> None:
        """Insert a column into this table.

        Parameters
        ----------
        column_name : str
        values : sequence of values

        """
        self.db.insert_column_into_table(self, column_name, values)

    def to_sql(self) -> str:
        """Return the SQLite query that generates this table."""
        return self.SELECT("*").to_sql()

    def to_df(self) -> pd.DataFrame:
        """Return a `pandas.DataFrame` of this table."""
        return pd.read_sql(self.to_sql(), self.db.connection)

    @property
    def info(self) -> pd.DataFrame:
        """The PRAGMA info table as a `pandas.DataFrame`."""
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
        """The schema/structure of this table."""
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
            self.db._table_list
            .SELECT("ncol")
            .WHERE(f"name='{self.name}'")
            .execute()
            .fetchone()
            .ncol
        )

    @property
    def rows(self) -> Iterator[Row]:
        """An iterable of Rows contained in this table.

        TODO: implement this as a UserDict container?

        """
        return (Row(id, self) for id in self._get_rowids())

    @property
    def columns(self) -> Iterator[Column]:
        """An iterable of Columns contained in this table.

        TODO: implement this as a UserDict container?

        """
        return (Column(name, self) for name in self._get_column_names())

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
        """The primary key of this table."""
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

    def DELETE(self) -> query.Query:
        """SQLite DELETE statement entry point via table.

        Returns
        -------
        query.Query

        """
        return self.db.DELETE().FROM(self.name)

    def SELECT(self, *columns: str) -> query.Query:
        """SQlite SELECT statement entry point via table.

        Parameters
        ----------
        *columns : str

        Returns
        -------
        query.Query

        """
        return self.db.SELECT(*columns).FROM(self.name)

    def SELECT_COUNT(self, column: str) -> query.Query:
        """SQLite SELECT COUNT statement entry point via table.

        Parameters
        ----------
        column : str

        Returns
        -------
        query.Query

        """
        return self.db.SELECT_COUNT(column).FROM(self.name)

    def SELECT_DISTINCT(self, *columns: str) -> query.Query:
        """SQlite SELECT DISTINCT statement entry point via table.

        Parameters
        ----------
        *columns : str

        Returns
        -------
        query.Query

        """
        return self.db.SELECT_DISTINCT(*columns).FROM(self.name)

    def SET(self, *column_exprs) -> query.Query:
        """SQLite SET statement entry point via table.

        Parameters
        ----------
        *column_exprs : str

        Returns
        -------
        query.Query

        """
        return self.db.UPDATE(self.name).SET(*column_exprs)

    def __len__(self) -> int:
        return self.n_rows

    # def __iter__(self) -> Iterator[NamedTuple]:
    #     return iter(self.rows)


class Schema:
    def __init__(self, table_name: str, *column_defs: str):
        self.table_name = table_name
        self.column_defs = column_defs

    def __repr__(self):
        return self.to_sql()

    def to_sql(self):
        col_defs = ",\n    ".join(self.column_defs)
        statement = f"{self.table_name} (\n    {col_defs}\n)"
        return statement

    @classmethod
    def from_str(cls, string):
        table_name = string[:string.find("(")].strip()
        column_defs = [
            s.strip()
            for s in string[string.find("(")+1:string.rfind(")")].replace("\n", "").split(",")
        ]
        return cls(table_name, *column_defs)


@dataclass
class Row:
    id: int
    table: Table

    @property
    def db(self) -> Database:
        return self.table.db

    @property
    def data(self):
        return self._q.execute().fetchone()

    def to_sql(self) -> str:
        return self._q.to_sql()

    def to_df(self) -> pd.DataFrame:
        return pd.read_sql(self.to_sql(), self.db.connection)

    @property
    def _q(self) -> query.Query:
        return (
            self.table
            .SELECT("*")
            .WHERE(f"{self.table.pk}={self.id}")
        )

    def __getattr__(self, name):
        return getattr(self.data, name)


@dataclass
class Column:
    name: str
    table: Table

    @property
    def db(self) -> Database:
        return self.table.db

    @property
    def type_(self) -> str:
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
        return [row[0] for row in self._q.execute().fetchall()]

    def to_sql(self) -> str:
        return self._q.to_sql()

    def to_df(self) -> pd.DataFrame:
        return pd.read_sql(self.to_sql(), self.db.connection)

    @property
    def _q(self) -> query.Query:
        return self.table.SELECT(self.name)


def touch(dbfile: pathlib.Path | str) -> None:
    """Create a minimal database file.

    A legal database includes a 'Simulations' and 'Observables' table
    at a minimum.

    Parameters
    ----------
    dbfile : path-like
        Path to database file.

    Raises
    ------
    ValueError
        If database already exists.

    """
    if isinstance(dbfile, str):
        dbfile = pathlib.Path(dbfile)

    if dbfile.exists():
        raise ValueError(f"{dbfile} already exists")

    sims_schema = """
    Simulations (
        simID INT PRIMARY KEY,
        topology TEXT,
        trajectory TEXT
    )
    """
    obsv_schema = """
    Observables (
        obsName TEXT,
        notes TEXT,
        creator TEXT,
        timestamp DATETIME DEFAULT (strftime('%m-%d-%Y %H:%M', 'now', 'localtime'))
    )
    """

    with Database(dbfile) as db:
        db.create_table(sims_schema)
        db.create_table(obsv_schema, STRICT=False)


# class Tables(UserDict):
#     def __init__(self, db, *args, **kwargs):
#         self.db = db
