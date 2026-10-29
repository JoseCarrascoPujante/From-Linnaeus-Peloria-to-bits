import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from epibit import entropy_binary, kl_binary, corrected_probability, jsd_binary

def test_entropy():
    assert abs(entropy_binary(0.5) - 1.0) < 1e-12
    assert entropy_binary(0.0) == 0.0
    assert entropy_binary(1.0) == 0.0

def test_probability():
    assert abs(corrected_probability(0, 10, 0.5) - (0.5/11.0)) < 1e-12

def test_divergence_zero():
    p = 0.37
    assert abs(kl_binary(p, p)) < 1e-12
    assert abs(jsd_binary(p, p)) < 1e-12

if __name__ == "__main__":
    test_entropy()
    test_probability()
    test_divergence_zero()
    print("All tests passed.")
