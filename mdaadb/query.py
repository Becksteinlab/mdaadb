from functools import wraps
from typing import Self
import pandas as pd


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

def register_fields(prefix):

    print(prefix)
    def decorator(func):

        @wraps(func)
        def wrapper(self, *fields):
            print(f"registering {prefix} fields: {fields}")
            self._registered_fields[prefix] = []
            for field in fields:
                if not isinstance(field, str):
                    raise ValueError("only strings work atm")
                self._registered_fields[prefix].append(field)
            return func(self)

        return wrapper

    return decorator


class Query:

    def __init__(self, db):
        self.db = db
        self._registered_fields = {}

    @register_fields("SELECT")
    def select(self, *columns) -> Self:
        return self

    # @register_fields("SELECT")
    # def select_table(self, *tables) -> Self:
    #     return self

    @register_fields("PRAGMA")
    def pragma(self, *fields) -> Self:
        return self

    @register_fields("FROM")
    def from_(self, *tables) -> Self:
        return self

    @register_fields("WHERE")
    def where(self, *conditions) -> Self:
        return self

    @register_fields("INSERT_INTO")
    def insert_into(self, *tables) -> Self:
        return self

    @register_fields("VALUES")
    def values(self, *values) -> Self:
        return self

    @register_fields("INNER_JOIN")
    def inner_join(self, *joins) -> Self:
        return self

    @register_fields("DELETE")
    def delete(self, *statements) -> Self:
        return self

    @register_fields("CREATE_TABLE")
    def create_table(self, *fields) -> Self:
        return self

    @register_fields("ORDER_BY")
    def order_by(self, *orders) -> Self:
        return self

    def to_sql(self) -> str:
        """"""
        statements = []
        for prefix, fields in self._registered_fields.items():
            statement = prefix + " " + ",".join(fields)
            statements.append(statement)

        sort_by = (
            "SELECT",
            "PRAGMA",
            "FROM",
            #"AS",
            "WHERE",
            "INNER_JOIN",
            "INSERT_INTO",
            "VALUES",
            "DELETE",
            "CREATE_TABLE",
            "ORDER_BY"
        )
        statements = sorted(statements, key=lambda x: sort_by.index(x.split(" ")[0]))
        statements = [statement.replace("_", " ") for statement in statements]
        sql_query = "\n".join(statements) + ";"

        return sql_query

    def to_df(self) -> pd.DataFrame:
        """"""
        return pd.read_sql(self.to_sql(), self.db.connection)

    def execute(self):
        """"""
        return self.db.cursor.execute(self.to_sql())

    def executemany(self, rows):
        """"""
        return self.db.cursor.executemany(self.to_sql(), rows)

    def commit(self):
        """"""
        return self.db.connection.commit()
