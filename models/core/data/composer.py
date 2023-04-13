from typing import List, Callable, Optional

from pybedtools import Interval
from torch.utils.data import Dataset

from .sources import DataSource


class Composer(Dataset):
    def __init__(
        self, features: dict[str, DataSource], targets: dict[str, DataSource], windows: List[Interval],
        window_transform: Optional[Callable[[Interval], Interval]] = None,
        features_transform: Optional[Callable] = None,
        target_transform: Optional[Callable] = None,
    ):
        if not features or not targets or not windows:
            raise ValueError('There must be at least 1 data source for features/targets and at least 1 genomic window.')
        self.features = features
        self.targets = targets
        self.windows = windows
        self.window_transform = window_transform
        self.features_transfrom = features_transform
        self.target_transform = target_transform

    def __getitem__(self, idx):
        win = self.windows[idx]
        if self.window_transform is not None:
            win = self.window_transform(win)

        features = {key: ds.fetch(win.chrom, win.strand, win.start, win.end) for key, ds in self.features.items()}
        if self.features_transfrom is not None:
            features = self.features_transfrom(features)

        targets = {key: ds.fetch(win.chrom, win.strand, win.start, win.end) for key, ds in self.targets.items()}
        if self.target_transform is not None:
            targets = self.target_transform(targets)

        return {'features': features, 'targets': targets}

    def __len__(self):
        return len(self.windows)
