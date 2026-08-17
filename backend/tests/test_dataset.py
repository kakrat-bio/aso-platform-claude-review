import os

import pytest

from backend.datasets.huesken import HueskenDataset
from backend.datasets.torch_dataset import ASODataset, ASOEmbeddingDataset
from backend.features.embed_cache import CACHE_FILE

from torch.utils.data import DataLoader


def test_embedding_dataset():
    if not os.path.exists(CACHE_FILE):
        pytest.skip(
            f"no RNA-FM embedding cache at {CACHE_FILE} — run "
            "features.embed_cache.precompute_embeddings first. Returning "
            "silently here would report PASSED for a test that ran nothing."
        )

    dataset = HueskenDataset("OligoFormer/data/Hu.csv")
    emb_dataset = ASOEmbeddingDataset(dataset, cache_path=CACHE_FILE)

    x, y = emb_dataset[0]

    print("x.shape:", x.shape)
    print("y:", y)

    assert x.shape == (1291,)
    assert y.shape == ()


def test_embedding_dataloader():
    if not os.path.exists(CACHE_FILE):
        pytest.skip(
            f"no RNA-FM embedding cache at {CACHE_FILE} — run "
            "features.embed_cache.precompute_embeddings first. Returning "
            "silently here would report PASSED for a test that ran nothing."
        )

    dataset = HueskenDataset("OligoFormer/data/Hu.csv")
    emb_dataset = ASOEmbeddingDataset(dataset, cache_path=CACHE_FILE)

    loader = DataLoader(
        emb_dataset,
        batch_size=32,
        shuffle=True,
    )

    for X, y in loader:
        print("X.shape:", X.shape)
        print("y.shape:", y.shape)

        assert X.shape == (32, 1291)
        assert y.shape == (32,)
        break


def test_torch_dataset():
    dataset = HueskenDataset("OligoFormer/data/Hu.csv")
    torch_dataset = ASODataset(dataset)

    x, y = torch_dataset[0]

    print("x:", x)
    print("x.shape:", x.shape)
    print("y:", y)

    assert x.shape == (9,)
    assert y.shape == ()


def test_dataloader():
    dataset = HueskenDataset("OligoFormer/data/Hu.csv")
    torch_dataset = ASODataset(dataset)

    loader = DataLoader(
        torch_dataset,
        batch_size=32,
        shuffle=True,
    )

    for X, y in loader:
        print("X.shape:", X.shape)
        print("y.shape:", y.shape)

        assert X.shape == (32, 9)
        assert y.shape == (32,)
        break


if __name__ == "__main__":
    test_embedding_dataset()
    test_embedding_dataloader()
    test_torch_dataset()
    test_dataloader()
    print("All dataset tests passed!")
