import pytest
from numpy.testing import assert_equal
from mdaadb.database import Database, Table, Schema


class TestMemoryDatabase:

    @pytest.fixture()
    def memdb(self):
        db = Database(":memory:")
        tbl1 = db.create_table(
            """
            Table1 (
                ID INT PRIMARY KEY,
                text_col TEXT,
                real_col REAL
            )
            """
        )
        tbl2 = db.create_table(
            """
            Table2 (
                ID INT PRIMARY KEY,
                text_col TEXT,
                real_col REAL
            )
            """
        )
        rows = [(i, f"text{i}", i*10.0) for i in range(1,11)]
        tbl1.insert_rows(rows)
        for row in rows:
            tbl2.insert_row(row)

        return db

    def test_get_table(self, memdb):
        tbl1 = memdb.get_table("Table1")
        tbl2 = memdb.get_table("Table2")
        assert tbl1 in memdb
        assert tbl2 in memdb

    def test_nrows(self, memdb):
        assert memdb.get_table("Table1").n_rows == 10
        assert memdb.get_table("Table2").n_rows == 10

    def test_ncols(self, memdb):
        assert memdb.get_table("Table1").n_cols == 3
        assert memdb.get_table("Table2").n_cols == 3

    def test_table_names(self, memdb):
        actual = memdb._get_table_names()
        desired = ["Table2", "Table1", "sqlite_schema"]
        assert_equal(actual, desired)

    def test_insert_rows(self, memdb):
        rows = [(i, f"text{i}", float(i*10)) for i in range(1,11)]
        rows_from_db = [
            tuple(row)
            for row in memdb.get_table("Table1").SELECT("*").execute().fetchall()
        ]

        assert_equal(rows_from_db, rows)

    def test_insert_row_vs_rows(self, memdb):
        tbl1 = memdb.get_table("Table1")
        tbl2 = memdb.get_table("Table2")
        result1 = tbl1.SELECT("*").execute().fetchall()
        result2 = tbl2.SELECT("*").execute().fetchall()

        assert_equal(result1, result2)

    def test_add_column(self, memdb):
        ...

    def test_update_column(self, memdb):
        ...

    def test_insert_column(self, memdb):
        ...

    def test_Rows(self, memdb):
        tbl1 = memdb.get_table("Table1")
        tbl2 = memdb.get_table("Table2")

        for i, row in enumerate(tbl1.rows):
            assert row.ID == i + 1
            assert row.text_col == f"text{i+1}"
            assert row.real_col == (i+1) * 10.0
        for i, row in enumerate(tbl2.rows):
            assert row.ID == i + 1
            assert row.text_col == f"text{i+1}"
            assert row.real_col == (i+1) * 10.0

    def test_Columns(self, memdb):
        tbl1 = memdb.get_table("Table1")
        tbl2 = memdb.get_table("Table2")

        for name, type_, col in zip(
            ["ID", "text_col", "real_col"],
            ["INT", "TEXT", "REAL"],
            tbl1.columns,
        ):
            assert col.type_ == tbl1.get_column(name).type_
        for name, type_, col in zip(
            ["ID", "text_col", "real_col"],
            ["INT", "TEXT", "REAL"],
            tbl2.columns,
        ):
            assert col.type_ == tbl1.get_column(name).type_

