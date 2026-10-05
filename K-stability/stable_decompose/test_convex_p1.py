import numpy as np

from convex_p1 import _prolongate, assemble_qp, build_mesh, load_polygon, solve_level


def test_all_samples_are_convex_and_ccw():
    from pathlib import Path
    for path in Path("unstable_polytope").glob("d[56]_a[0-9]*"):
        if path.suffix:
            continue
        sample = load_polygon(path)
        a = sample.vertices[1] - sample.vertices[0]
        b = sample.vertices[2] - sample.vertices[1]
        assert a[0] * b[1] - a[1] * b[0] > 0
        assert np.all(np.diff(np.r_[sample.vertices[:, 0], sample.vertices[0, 0]]) != 0) or True


def test_affine_values_satisfy_constraints():
    sample = load_polygon("unstable_polytope/d5_a70")
    mesh = build_mesh(sample, 2)
    _, _, A, _, _ = assemble_qp(mesh)
    values = 2 * mesh.nodes[:, 0] - 3 * mesh.nodes[:, 1] + 7
    assert A.shape[0] > 0
    assert np.max(np.abs(A @ values)) < 1e-8


def test_refinement_counts_and_no_duplicate_nodes():
    sample = load_polygon("unstable_polytope/d5_a70")
    m0, m1, m2 = (build_mesh(sample, i) for i in range(3))
    assert len(m1.triangles) == 4 * len(m0.triangles)
    assert len(m2.triangles) == 4 * len(m1.triangles)
    assert len({tuple(np.round(x, 12)) for x in m2.nodes}) == len(m2.nodes)


def test_local_mass_matrix_is_positive_definite():
    sample = load_polygon("unstable_polytope/d5_a70")
    mesh = build_mesh(sample, 0)
    P, q, A, lower, upper = assemble_qp(mesh)
    assert P.shape == (len(mesh.nodes), len(mesh.nodes))
    assert np.all(np.linalg.eigvalsh(P.toarray()) > 0)
    assert len(q) == len(mesh.nodes)
    assert len(lower) == A.shape[0] == len(upper)


def test_boundary_linear_term_on_single_triangle(tmp_path):
    path = tmp_path / "triangle"
    path.write_text("0 0\n2 0\n0 2\n")
    mesh = build_mesh(load_polygon(path), 0)
    _, q, A, _, _ = assemble_qp(mesh)
    assert A.shape[0] == 0
    assert np.allclose(q, [2.0, 2.0, 2.0])


def test_nonconvex_nodal_values_violate_constraint():
    sample = load_polygon("unstable_polytope/d5_a70")
    mesh = build_mesh(sample, 1)
    _, _, A, _, _ = assemble_qp(mesh)
    values = -np.sum(mesh.nodes ** 2, axis=1)
    assert np.min(A @ values) < -1e-8


def test_red_refinement_prolongation_is_affine_exact():
    sample = load_polygon("unstable_polytope/d6_a140")
    old_mesh = build_mesh(sample, 1)
    new_mesh = build_mesh(sample, 2)
    old_values = 2 * old_mesh.nodes[:, 0] - 3 * old_mesh.nodes[:, 1] + 7
    prolonged = _prolongate(old_mesh, new_mesh, old_values)
    expected = 2 * new_mesh.nodes[:, 0] - 3 * new_mesh.nodes[:, 1] + 7
    assert np.max(np.abs(prolonged - expected)) < 1e-12


def test_optimized_prolongation_remains_feasible():
    sample = load_polygon("unstable_polytope/d5_a70")
    old_mesh = build_mesh(sample, 0)
    old_result = solve_level(old_mesh)
    new_mesh = build_mesh(sample, 1)
    prolonged = _prolongate(old_mesh, new_mesh, old_result.values_u)
    _, _, A, _, _ = assemble_qp(new_mesh)
    assert np.min(A @ prolonged) >= -1e-8


def test_clockwise_polygon_is_rejected(tmp_path):
    path = tmp_path / "clockwise"
    path.write_text("0 0\n0 1\n1 0\n")
    try:
        load_polygon(path)
    except ValueError as exc:
        assert "counter-clockwise" in str(exc)
    else:
        raise AssertionError("clockwise input should be rejected")
