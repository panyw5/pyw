import sqlite3
import tempfile
from pathlib import Path

import pytest


@pytest.mark.sage
class TestSQLiteCacheInit:
    def test_init_creates_db_and_table(self):
        from sage.all import WeylGroup
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        with tempfile.TemporaryDirectory() as tmpdir:
            W = WeylGroup(["A", 2])
            kl = KazhdanLusztigPolynomials(W, cache_dir=Path(tmpdir), persistent_cache=True)

            assert kl._db_conn is not None
            db_path = str(kl._db_path)
            assert Path(db_path).exists()

            db = sqlite3.connect(db_path)
            assert db.execute("PRAGMA journal_mode").fetchone()[0] == "wal"
            tables = db.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='q_at_one'"
            ).fetchall()
            assert len(tables) == 1
            db.close()

    def test_init_loads_existing_entries(self):
        from sage.all import WeylGroup
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        with tempfile.TemporaryDirectory() as tmpdir:
            W = WeylGroup(["A", 2])
            identity_word = tuple()
            s1_word = (1,)

            kl1 = KazhdanLusztigPolynomials(
                W, cache_dir=Path(tmpdir), persistent_cache=True, auto_save_min_new_entries=0
            )
            identity = W.one()
            s1 = W.simple_reflection(1)
            q_val = kl1.Q(identity, s1, at_one=True)
            kl1._db_conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            kl1._db_conn.close()

            kl2 = KazhdanLusztigPolynomials(W, cache_dir=Path(tmpdir), persistent_cache=True)
            cached = kl2._Q_at_one_cache.get((identity_word, s1_word))
            assert cached is not None
            assert cached == q_val


@pytest.mark.sage
class TestSQLiteCachePersistence:
    def test_persist_q_at_one_writes_immediately(self):
        from sage.all import WeylGroup
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        with tempfile.TemporaryDirectory() as tmpdir:
            W = WeylGroup(["A", 2])
            kl = KazhdanLusztigPolynomials(
                W, cache_dir=Path(tmpdir), persistent_cache=True, auto_save_min_new_entries=99999
            )

            identity = W.one()
            s1 = W.simple_reflection(1)
            q_val = kl.Q(identity, s1, at_one=True)

            rows = kl._db_conn.execute("SELECT count(*) FROM q_at_one").fetchone()
            assert rows[0] >= 1

    def test_cache_hit_on_second_lookup(self):
        from sage.all import WeylGroup
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        with tempfile.TemporaryDirectory() as tmpdir:
            W = WeylGroup(["A", 2])
            kl = KazhdanLusztigPolynomials(
                W, cache_dir=Path(tmpdir), persistent_cache=True, auto_save_min_new_entries=99999
            )
            kl.set_profiling(True)

            identity = W.one()
            s1 = W.simple_reflection(1)

            kl.Q(identity, s1, at_one=True)
            kl.Q(identity, s1, at_one=True)

            stats = kl.profile_stats()
            assert stats["Q_cache_hits_at_one"] >= 1

    def test_insert_or_replace_handles_duplicates(self):
        from sage.all import WeylGroup
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        with tempfile.TemporaryDirectory() as tmpdir:
            W = WeylGroup(["A", 2])
            kl = KazhdanLusztigPolynomials(W, cache_dir=Path(tmpdir), persistent_cache=True)

            cache_key = ((1, 2), (1,))
            kl._persist_q_at_one(cache_key, 42)
            kl._persist_q_at_one(cache_key, 99)

            rows = kl._db_conn.execute("SELECT count(*) FROM q_at_one").fetchone()
            assert rows[0] == 1


@pytest.mark.sage
class TestSQLiteCacheConcurrency:
    def test_two_instances_share_cache(self):
        from sage.all import WeylGroup
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        with tempfile.TemporaryDirectory() as tmpdir:
            W1 = WeylGroup(["A", 2])
            W2 = WeylGroup(["A", 2])

            kl1 = KazhdanLusztigPolynomials(
                W1, cache_dir=Path(tmpdir), persistent_cache=True, auto_save_min_new_entries=0
            )
            identity1 = W1.one()
            s1_1 = W1.simple_reflection(1)
            q_val = kl1.Q(identity1, s1_1, at_one=True)
            kl1._db_conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            kl1._db_conn.close()

            kl2 = KazhdanLusztigPolynomials(
                W2, cache_dir=Path(tmpdir), persistent_cache=True, auto_load_cache=True
            )
            assert len(kl2._Q_at_one_cache) >= 1


@pytest.mark.sage
class TestSQLiteDisabled:
    def test_persistent_cache_disabled_no_db(self):
        from sage.all import WeylGroup
        from pyw.core.kazhdan_lusztig import KazhdanLusztigPolynomials

        with tempfile.TemporaryDirectory() as tmpdir:
            W = WeylGroup(["A", 2])
            kl = KazhdanLusztigPolynomials(
                W, cache_dir=Path(tmpdir), persistent_cache=False, auto_load_cache=False
            )

            assert kl._db_conn is None


@pytest.mark.sage
class TestCharacterCaching:
    def test_character_reuses_cache_across_instances(self):
        from pyw.core.affine_lie_algebra import AffineLieAlgebra
        from pyw.core.affine_weight import AffineWeight
        from pyw.core.character import KazhdanLusztigCharacter

        ala1 = AffineLieAlgebra(["A", 2, 1])
        lam1 = AffineWeight.affine_fundamental_weight(ala1, 1)
        klc1 = KazhdanLusztigCharacter(ala1)
        klc1.kl_polynomial.set_profiling(True)
        klc1.character(lam1, order=1)
        stats_first = klc1.kl_polynomial.profile_stats()

        ala2 = AffineLieAlgebra(["A", 2, 1])
        lam2 = AffineWeight.affine_fundamental_weight(ala2, 1)
        klc2 = KazhdanLusztigCharacter(ala2)
        klc2.kl_polynomial.set_profiling(True)
        klc2.character(lam2, order=1)
        stats_second = klc2.kl_polynomial.profile_stats()

        assert stats_second["Q_cache_hits_at_one"] > 0
        assert stats_second["Q_cache_hits_at_one"] == stats_second["Q_calls"]
