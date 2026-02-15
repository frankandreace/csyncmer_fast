"Fast closed syncmer detection using ntHash rolling hash."

from ._bindings import (
    SyncmerIterator,
    CanonicalSyncmerIterator,
    count_syncmers,
    count_syncmers_canonical,
)

__version__ = "0.2.0"

__all__ = [
    'SyncmerIterator',
    'CanonicalSyncmerIterator',
    'count_syncmers',
    'count_syncmers_canonical',
]
