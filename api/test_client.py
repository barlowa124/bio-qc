"""Test client for the prediction API."""
import json

import numpy as np
import requests

# Example: random expression values for 30 genes
test_expr = np.random.lognormal(0, 1, 30).tolist()

# Local test
response = requests.post("http://localhost:5000/predict", json={"expression": test_expr})
print(json.dumps(response.json(), indent=2))
