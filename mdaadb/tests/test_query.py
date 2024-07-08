import pytest
from ..query import Query


class TestQuery:



    def test_SELECT(self):
        correct_q = "SELECT simID\nFROM simulations;"
        q = Query(db="dummy").SELECT("simID").FROM("simulations").to_sql()

        assert q == correct_q

    def test_SELECT_missing_FROM(self, q):
        with pytest.raises(ValueError):
            q = Query(db="dummy").SELECT("simID").to_sql()

    # def test_FROM_missing_select(self, q):
    #     with pytest.raises(ValueError):
    #         q = Query(db="dummy").FROM("simulations").to_sql()
