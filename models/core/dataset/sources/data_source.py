from abc import ABC, abstractmethod
from typing import Literal

import numpy as np


class DataSource(ABC):
    @abstractmethod
    def fetch(self, contig: str, strand: Literal["+", "-", "."], start: int, end: int) -> np.ndarray:
        raise NotImplementedError()
