"""R01: sanity check that adversarial_opt.n3_objective is non-trivial (non-Einstein points > bound)."""
import sys, math, time
import numpy as np
sys.argv = [sys.argv[0], sys.argv[1], "/dev/null"]
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import adversarial_opt as ao
for x in ([0.0, 0.0], [math.log(3), 0.0], [math.log(8), 0.0], [-math.log(8), 1.0], [math.log(50), -1.0]):
    t = time.time()
    print(x, ao.n3_objective(np.array(x), 2, 3.0, 0.0, math.inf), f"{time.time()-t:.2f}s", flush=True)
