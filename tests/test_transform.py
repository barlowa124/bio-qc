import numpy as np
import pandas as pd

from cytof_qc import transform


def test_arcsinh_roundtrip():
    df = pd.DataFrame({"CD45": [0.0, 1.0, 50.0, -2.0], "CD3": [10.0, 0.0, -1.0, 4.0]})
    rt = transform.inverse_arcsinh(transform.arcsinh(df))
    np.testing.assert_allclose(rt.to_numpy(), df.to_numpy(), rtol=1e-9)


def test_arcsinh_handles_negatives():
    df = pd.DataFrame({"x": [-100.0, -1.0, 0.0, 1.0]})
    out = transform.arcsinh(df)
    assert np.isfinite(out["x"]).all()
    assert out["x"].iloc[0] < 0


def test_arcsinh_rejects_nonpositive_cofactor():
    import pytest

    with pytest.raises(ValueError, match="cofactor"):
        transform.arcsinh(pd.DataFrame({"x": [1.0]}), cofactor=0)
