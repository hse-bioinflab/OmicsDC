from typing import Union

from pybedtools import Interval, BedTool

BedLike = Union[BedTool, list[Interval]]
