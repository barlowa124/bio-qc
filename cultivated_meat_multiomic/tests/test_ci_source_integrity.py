import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def generate(root):
    shutil.copy(ROOT / 'ci_fixtures.py', root / 'ci_fixtures.py')
    result = subprocess.run([sys.executable, str(root / 'ci_fixtures.py')], cwd=root,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_generator_does_not_create_api_source(tmp_path):
    generate(tmp_path)
    assert not (tmp_path / 'api' / 'app.py').exists()


def test_generator_preserves_existing_api_source(tmp_path):
    api = tmp_path / 'api'
    api.mkdir()
    source = api / 'app.py'
    source.write_bytes(b'SENTINEL = True\n')
    generate(tmp_path)
    assert source.read_bytes() == b'SENTINEL = True\n'


def test_production_api_uses_fixture_data_without_replacement(tmp_path):
    pytest.importorskip('flask')
    api = tmp_path / 'api'
    api.mkdir()
    source = api / 'app.py'
    shutil.copy(ROOT / 'api' / 'app.py', source)
    before = source.read_bytes()
    generate(tmp_path)
    assert source.read_bytes() == before
    spec = importlib.util.spec_from_file_location('isolated_fixture_api', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    client = module.app.test_client()
    assert client.get('/health').get_json()['model_accuracy'] is None
    assert client.get('/health').get_json()['data_source'] == 'ci_fixture_synthetic'
    valid = client.post('/predict', json={'expression': [1.0] * 30})
    assert valid.status_code == 200
    assert valid.get_json()['data_source'] == 'ci_fixture_synthetic'
    assert client.post('/predict', json={'expression': [float('nan')] * 30}).status_code == 400
