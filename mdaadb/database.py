from __future__ import annotations

import sqlite3
import pandas as pd
from dataclasses import dataclass
from collections import namedtuple
from typing import List, Iterable, Any
from query import Query


def _namedtuple_factory(cursor, row):
    fields = [col[0] for col in cursor.description]
    Row = namedtuple("Row", fields)
    return Row(*row)


class Database:

    def __init__(self, database):
        self.connection = sqlite3.connect(database)
        self.connection.row_factory = _namedtuple_factory
        self.cursor = self.connection.cursor()

    def table(self, name):
        """"""
        if name not in self._get_table_names():
            raise ValueError("invalid table")
        return Table(name, self)

    @property
    def tables(self):
        """"""
        return (Table(name, self) for name in self._get_table_names())

    def _get_table_names(self) -> List[str]:
        """"""
        table_names = (
            self.select("name")
            .from_("sqlite_schema")
            .where("type = 'table'")
            .execute()
            .fetchall()
        )
        return [name[0] for name in table_names]

    @property
    def schema(self) -> pd.DataFrame:
        return Query(db=self).select("*").from_("sqlite_schema").to_df()

    def select(self, *fields):
        """"""
        return Query(db=self).select(*fields)

    def pragma(self, *fields):
        """"""
        return Query(db=self).pragma(*fields)

    def insert_into(self, *tables):
        """"""
        return Query(db=self).insert_into(*tables)

    def create_table(self, *fields):
        """"""
        return Query(db=self).create_table(*fields)

    def insert_row_into_table(self, table, row) -> None:
        if isinstance(table, str):
            table = Table(table, self)
        table.insert_row(row)

    def insert_array_into_table(self, table, array) -> None:
        if isinstance(table, str):
            table = Table(table, self)
        table.insert_array(array)

    def __iter__(self):
        return iter(self.tables)

    def __contains__(self, table):
        return table.database == self


@dataclass
class Table:
    name: str
    db: Database

    @property
    def schema(self) -> str:
        """"""
        return self.db.select("sql").from_("sqlite_schema").where(f"name = '{self.name}'").execute().fetchone()[0]

    @property
    def info(self) -> pd.DataFrame:
        """"""
        return self.db.pragma(f"table_info('{self.name}')").to_df()

    def row(self, id:int) -> Row:
        """"""
        if id not in self._get_row_ids():
            raise ValueError("id not in valid ids")
        return Row(id, self)

    def column(self, name:str) -> Column:
        """"""
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
        row_ids = self.select(self.primary_key).execute().fetchall()
        return [id[0] for id in row_ids]

    def _get_column_names(self) -> List[str]:
        """Helper function that returns a list of column names for this table.."""
        column_names = (
            self.db.select("name")
            .from_(f"pragma_table_info('simulations')")
            .execute()
            .fetchall()
        )
        return [name[0] for name in column_names]

    @property
    def primary_key(self) -> str:
        """Returns the primary key of this table as a string."""
        result = (
            self.db.select("name")
            .from_(f"pragma_table_info('{self.name}') as ti")
            .where("ti.pk = 1")
            .execute()
            .fetchall()
        )
        assert len(result[0]) == 1
        return result[0][0]

    def select(self, *columns) -> Query:
        """"""
        return self.db.select(*columns).from_(self.name)

    def insert_row(self, row) -> None:
        """"""
        return self.db.insert_into(self.name).values(row).execute().commit()

    def insert_array(self, array) -> None:
        """"""
        n_cols = self.n_cols
        values = "(" + ", ".join(["?" for _ in range(n_cols)]) + ")"
        return (
            self.db
            .insert_into(self.name)
            .values(values)
            .executemany(array)
            .commit()
        )

    def delete_row(self, row) -> None:
        ...

    def delete_row_by_id(self, id) -> None:
        ...

    def to_sql(self) -> str:
        """Returns the sqlite query that generates this table."""
        return self.select("*").to_sql()

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
            .select("*")
            .where(f"{self.table.primary_key}={self.id}")
            .execute()
            .fetchall()
        )

    def to_sql(self) -> str:
        """Returns the sqlite query that generates this row."""
        return (
            self.table
            .select("*")
            .where(f"{self.table.primary_key}={self.id}")
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
            .select("type")
            .from_(f"pragma_table_info('{self.table.name}') as t_info")
            .where(f"t_info.name='{self.name}'")
            .execute()
            .fetchall()
        )
        return _type[0][0]

    @property
    def data(self) -> List[Any]:
        """"""
        _data = self.table.select(self.name).execute().fetchall()
        return [data[0] for data in _data]

    def to_sql(self) -> str:
        """"""
        return self.table.select(self.name).to_sql()

    def to_df(self) -> pd.DataFrame:
        """"""
        return pd.read_sql(self.to_sql(), self.db.connection)

