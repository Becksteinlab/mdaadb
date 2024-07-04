from functools import wraps
from collections import UserDict
import pandas as pd
try:
    from typing import Self
except ImportError:
    from typing_extensions import Self


# def validate_fields(expected_type):
#
#     print(expected_type)
#
#     def decorator(func):
#
#         @wraps(func)
#         def wrapper(self, *fields):
#             print(f"validating fields {fields}")
#             return func(self)
#
#         return wrapper
#
#     return decorator


def register_fields(func):

    @wraps(func)
    def wrapper(self, *fields):
        prefix = wrapper.__name__
        print(prefix)
        print(f"registering {prefix} fields: {fields}")
        self._registered_fields[prefix] = []
        for field in fields:
            if not isinstance(field, str):
                raise ValueError("only strings work atm")
            self._registered_fields[prefix].append(field)

        return func(self)

    return wrapper


class Fields(UserDict):
    ...


class Query:

    def __init__(self, db):
        self.db = db
        self._registered_fields = {}

    @register_fields
    def CREATE_TABLE(self, *fields) -> Self:
        return self

    @register_fields
    def DELETE(self, *fields) -> Self:
        return self

    @register_fields
    def FROM(self, *fields) -> Self:
        return self

    @register_fields
    def INNER_JOIN(self, *fields) -> Self:
        return self

    @register_fields
    def INTO(self, *fields) -> Self:
        return self

    @register_fields
    def INSERT(self, *fields) -> Self:
        return self

    @register_fields
    def LIMIT(self, *fields) -> Self:
        return self

    @register_fields
    def ON(self, *fields) -> Self:
        return self

    @register_fields
    def ORDER_BY(self, *fields) -> Self:
        return self

    @register_fields
    def PRAGMA(self, *fields) -> Self:
        return self

    @register_fields
    def SELECT(self, *fields) -> Self:
        return self

    @register_fields
    def VALUES(self, *fields) -> Self:
        return self

    @register_fields
    def WHERE(self, *fields) -> Self:
        return self

    def to_sql(self) -> str:
        """Converts the query from it's internal `Query` representation
        into a string that is the exact SQL query being executed."""

        statements = []
        for prefix, fields in self._registered_fields.items():
            statement = prefix + " " + ",".join(fields)
            statements.append(statement)

        sort_by = (
            "SELECT",
            "PRAGMA",
            "INSERT",
            "DELETE",
            "CREATE_TABLE",
            "VALUES",
            "FROM",
            "INTO",
            "INNER_JOIN",
            "ON",
            "WHERE",
            "ORDER_BY",
            "LIMIT",
        )

        statements = sorted(statements, key=lambda x: sort_by.index(x.split(" ")[0]))
        statements = [statement.replace("_", " ") for statement in statements]
        sql_query = "\n".join(statements) + ";"

        return sql_query

    def to_df(self) -> pd.DataFrame:
        """Executes the current query and outputs result as pandas DataFrame."""
        return pd.read_sql(self.to_sql(), self.db.connection)

    def execute(self):
        """Executes the current query and returns Cursor iterator."""
        return self.db.cursor.execute(self.to_sql())

    def executemany(self, rows):
        """"""
        return self.db.cursor.executemany(self.to_sql(), rows)

    def commit(self):
        """"""
        return self.db.connection.commit()

    def __repr__(self):
        return self.to_sql()
