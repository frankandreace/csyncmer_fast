import csyncmer_fast

SMALL_SEQUENCE = "ACGTAACGTATACGTA"
SMALL_K = 7
SMALL_S = 4

class TestFunctionality:
    """Testing that the library correctly returns closed syncmers."""

    def test_import(self):
        """Testing that the library is imported correctly."""
        
        assert hasattr(csyncmer_fast, '__version__')

    def test_returning_object_get_all_syncmer(self):
        """Testing that the library correctly returns a list of tuples."""
        
        m_closed_syncmer_iterator = csyncmer_fast.SyncmerIterator(SMALL_SEQUENCE,SMALL_K,SMALL_S)
        result_list: list = m_closed_syncmer_iterator.get_all_syncmers()

        assert isinstance(result_list, list)
        assert result_list != None
        assert len(result_list) > 0
        assert isinstance(result_list[0], tuple)

    def test_returing_object_next(self):
        """Testing that the __next__ function correctly return the tuple of (hash_value: int, k_mer_positon: int, is_forward: bool, s_mer_position_is_start: bool)."""

        m_closed_syncmer_iterator = csyncmer_fast.SyncmerIterator(SMALL_SEQUENCE,SMALL_K,SMALL_S)
        next_syncmer = m_closed_syncmer_iterator.__next__

        assert isinstance(next_syncmer, tuple)
        assert isinstance(next_syncmer[0], int)
        assert isinstance(next_syncmer[1], int)
        assert isinstance(next_syncmer[2], bool)
        assert isinstance(next_syncmer[3], bool)

        assert next_syncmer != None
        for value in next_syncmer:
            assert value != None
        

    def test_returning_object_iter(self):
        """Testing that the library correctly returns a an iterator."""
        m_closed_syncmer_iterator = csyncmer_fast.SyncmerIterator(SMALL_SEQUENCE,SMALL_K,SMALL_S)
        syncmer_iterator = m_closed_syncmer_iterator.__iter__

        assert isinstance(syncmer_iterator, iter) 

    def test_closed_syncmer_correctness_next(self):
        """Testing that the __next__ function correctly return the correct values for (hash_value: int, k_mer_positon: int, is_forward: bool, s_mer_position_is_start: bool)."""
        m_closed_syncmer_iterator = csyncmer_fast.SyncmerIterator(SMALL_SEQUENCE,SMALL_K,SMALL_S)
        next_syncmer = m_closed_syncmer_iterator.__next__

        assert next_syncmer[0] == 5
        assert next_syncmer[1] == 1
        assert next_syncmer[2] == False
        assert next_syncmer[3] == False

        next_syncmer = m_closed_syncmer_iterator.__next__

        assert next_syncmer[0] == 5
        assert next_syncmer[1] == 4
        assert next_syncmer[2] == False
        assert next_syncmer[3] == True

        next_syncmer = m_closed_syncmer_iterator.__next__

        assert next_syncmer[0] == 52
        assert next_syncmer[1] == 6
        assert next_syncmer[2] == False
        assert next_syncmer[3] == True

        next_syncmer = m_closed_syncmer_iterator.__next__

        assert next_syncmer[0] == 52
        assert next_syncmer[1] == 7
        assert next_syncmer[2] == False
        assert next_syncmer[3] == False


    def test_closed_syncmer_correctness_get_all_syncmer(self):
        """Testing that the result returned"""
        m_closed_syncmer_iterator_1 = csyncmer_fast.SyncmerIterator(SMALL_SEQUENCE,SMALL_K,SMALL_S)
        syncmer_list_1: list = m_closed_syncmer_iterator_1.get_all_syncmers()

        syncmer_list_2: list = list(csyncmer_fast.SyncmerIterator(SMALL_SEQUENCE,SMALL_K,SMALL_S))

        for syncmer_1,syncmer_2 in zip(syncmer_list_1, syncmer_list_2):
            for value1, value2 in zip(syncmer_1,syncmer_2):
                assert value1 == value2

