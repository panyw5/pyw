from sage.all import *
import pathlib
import sys
import time

sys.path.insert(0, "/Users/lelouch/pyw/refs/Kazhdan-Lusztig")
import Algebra as Alg


log = pathlib.Path("/Users/lelouch/pyw/.trellis/tasks/04-25-D6-level-(-4)/legacy_d6_omega0.log")


def write_line(handle, *args):
    msg = " ".join(str(a) for a in args)
    print(msg, flush=True)
    handle.write(msg + "\n")
    handle.flush()


with log.open("w") as f:
    write_line(f, "start")
    alg = Alg.Alg(["D", 6, 1], QLoad=False, WLoad=False)
    llambda = -4 * alg.omega[0]
    order = 2
    write_line(f, "llambda", llambda)

    t0 = time.time()
    numerator = alg.Kazhdan_Lusztig_numerator(llambda, order)
    write_line(f, "numerator_done_s", time.time() - t0)

    denominator = alg.Kazhdan_Lusztig_denominator(order)
    write_line(f, "denominator_done")

    q = var("q")
    character = simplify((numerator / denominator).taylor(q, 0, order))
    write_line(f, "character", character)
    write_line(f, "done")
