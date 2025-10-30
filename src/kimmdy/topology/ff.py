from __future__ import annotations

import logging
from gmx_top4py.topology.ff import FF as BasicFF

from kimmdy.parsing import read_edissoc


logger = logging.getLogger(__name__)


class FF(BasicFF):
    """Container for parsed forcefield data.

    Also see <https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#topology-file>
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.default_edissoc: dict[str, dict[tuple[str, str], float]] = {}
        if self.gmxdir is not None:
            gmx_builtin_ffs = self.gmxdir / "top"
            self.default_edissoc = read_edissoc(gmx_builtin_ffs / "edissoc.dat")
