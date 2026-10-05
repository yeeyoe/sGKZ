"""Create a self-contained interactive 3D HTML view from a solver NPZ file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _boundary_edges(triangles: np.ndarray) -> list[list[int]]:
    owners: dict[tuple[int, int], int] = {}
    for tri in triangles:
        for a, b in ((int(tri[0]), int(tri[1])),
                     (int(tri[1]), int(tri[2])),
                     (int(tri[2]), int(tri[0]))):
            edge = tuple(sorted((a, b)))
            owners[edge] = owners.get(edge, 0) + 1
    return [list(edge) for edge, count in owners.items() if count == 1]


def _crease_edges(triangles: np.ndarray, nodes: np.ndarray, gradients: np.ndarray,
                  threshold: float) -> list[list[int]]:
    """Return interior edges whose adjacent triangle gradients jump strongly."""
    owners: dict[tuple[int, int], list[int]] = {}
    for ti, tri in enumerate(triangles):
        for a, b in ((int(tri[0]), int(tri[1])),
                     (int(tri[1]), int(tri[2])),
                     (int(tri[2]), int(tri[0]))):
            owners.setdefault(tuple(sorted((a, b))), []).append(ti)
    marked = []
    for (a, b), adjacent in owners.items():
        if len(adjacent) == 2:
            jump = float(np.linalg.norm(gradients[adjacent[0]] - gradients[adjacent[1]]))
            if jump >= threshold:
                marked.append([a, b])
    return marked


def make_html(npz_path: str | Path, output_path: str | Path,
              title: str | None = None, z_scale: float = 8.0,
              crease_threshold: float | None = None) -> Path:
    """Write an interactive 3D surface view and return its path."""
    npz_path, output_path = Path(npz_path), Path(output_path)
    with np.load(npz_path) as data:
        nodes = np.asarray(data["nodes"], dtype=float)
        triangles = np.asarray(data["triangles"], dtype=np.int64)
        values = np.asarray(data["values_u"], dtype=float)
        gradients = np.asarray(data["triangle_gradients"], dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("NPZ nodes must have shape (N, 2)")
    if len(values) != len(nodes):
        raise ValueError("values_u length does not match nodes")
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError("NPZ triangles must have shape (T, 3)")
    if title is None:
        title = npz_path.stem.replace("_", " ")
    boundary = _boundary_edges(triangles)
    if crease_threshold is None:
        scale = float(np.max(np.linalg.norm(gradients, axis=1)))
        crease_threshold = max(1e-9, 1e-3 * scale)
    creases = _crease_edges(triangles, nodes, gradients, crease_threshold)
    payload = {
        "nodes": nodes.tolist(),
        "triangles": triangles.tolist(),
        "values": values.tolist(),
        "boundary": boundary,
        "creases": creases,
        "crease_threshold": float(crease_threshold),
        "title": title,
        "z_scale": float(z_scale),
        "psi_note": "values_u from convex_p1.py",
    }
    raw = json.dumps(payload, separators=(",", ":"))
    html = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  :root {{ color-scheme: light dark; }}
  body {{ margin: 0; font: 15px system-ui, sans-serif; background: Canvas; color: CanvasText; }}
  header {{ padding: 12px 16px 4px; }}
  h1 {{ font-size: 18px; font-weight: 500; margin: 0 0 4px; }}
  p {{ margin: 0; color: GrayText; }}
  #view {{ width: 100vw; height: calc(100vh - 60px); min-height: 360px; touch-action: none; }}
  canvas {{ display: block; width: 100%; height: 100%; }}
</style>
</head>
<body>
<header><h1 id="title"></h1><p>Drag to rotate. Scroll or pinch to zoom. Height is scaled for visibility. Axes: X red, Y green, Z blue.</p></header>
<div id="view" role="img" aria-label="Interactive three-dimensional triangular surface of the convex P1 minimizer"></div>
<script type="module">
import * as THREE from "https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js";
const data = {raw};
document.getElementById("title").textContent = data.title;
const host = document.getElementById("view");
const renderer = new THREE.WebGLRenderer({{antialias: true, alpha: true}});
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
host.appendChild(renderer.domElement);
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(38, 1, 0.1, 2000);
const points = data.nodes.map((p, i) => new THREE.Vector3(p[0], p[1], data.values[i] * data.z_scale));
const center = points.reduce((s, p) => s.add(p), new THREE.Vector3()).multiplyScalar(1 / points.length);
const model = new THREE.Group();
scene.add(model);
// Keep a world-space coordinate frame visible while the surface is rotated.
const axisLength = Math.max(
  ...points.map(p => Math.max(Math.abs(p.x), Math.abs(p.y), Math.abs(p.z))), 1
) * 1.15;
scene.add(new THREE.AxesHelper(axisLength));
const geometry = new THREE.BufferGeometry().setFromPoints(points);
geometry.setIndex(data.triangles.flat());
geometry.computeVertexNormals();
model.add(new THREE.Mesh(geometry, new THREE.MeshPhongMaterial({{color: 0x3478b5, side: THREE.DoubleSide, transparent: true, opacity: 0.85, shininess: 24}})));
model.add(new THREE.LineSegments(new THREE.EdgesGeometry(geometry), new THREE.LineBasicMaterial({{color: 0x667085, transparent: true, opacity: 0.5}})));
for (const [a, b] of data.boundary) {{
  const line = new THREE.BufferGeometry().setFromPoints([points[a], points[b]]);
  model.add(new THREE.Line(line, new THREE.LineBasicMaterial({{color: 0xd14b4b}})));
}}
for (const [a, b] of data.creases) {{
  const line = new THREE.BufferGeometry().setFromPoints([points[a], points[b]]);
  model.add(new THREE.Line(line, new THREE.LineBasicMaterial({{color: 0xff2414, linewidth: 3}})));
}}
scene.add(new THREE.HemisphereLight(0xffffff, 0x556070, 1.5));
const light = new THREE.DirectionalLight(0xffffff, 1.2);
light.position.set(20, -20, 40);
scene.add(light);
const span = Math.max(...points.map(p => p.distanceTo(center)), 1);
camera.position.copy(center).add(new THREE.Vector3(span * 1.2, -span * 1.6, span * 1.3));
camera.lookAt(center);
let dragging = false, lastX = 0, lastY = 0, distance = camera.position.distanceTo(center);
renderer.domElement.addEventListener("pointerdown", e => {{ dragging = true; lastX = e.clientX; lastY = e.clientY; renderer.domElement.setPointerCapture(e.pointerId); }});
renderer.domElement.addEventListener("pointermove", e => {{
  if (!dragging) return;
  model.rotation.z += (e.clientX - lastX) * 0.008;
  model.rotation.x += (e.clientY - lastY) * 0.008;
  lastX = e.clientX; lastY = e.clientY;
}});
renderer.domElement.addEventListener("pointerup", e => {{ dragging = false; renderer.domElement.releasePointerCapture(e.pointerId); }});
renderer.domElement.addEventListener("wheel", e => {{ e.preventDefault(); distance *= Math.exp(e.deltaY * 0.001); distance = Math.max(span * 0.25, Math.min(span * 8, distance)); camera.position.copy(center).add(camera.position.clone().sub(center).normalize().multiplyScalar(distance)); camera.lookAt(center); }}, {{passive: false}});
function resize() {{ const w = host.clientWidth, h = host.clientHeight; renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix(); }}
new ResizeObserver(resize).observe(host); resize();
function animate() {{ renderer.render(scene, camera); requestAnimationFrame(animate); }} animate();
</script>
</body>
</html>
'''
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="solver NPZ file")
    parser.add_argument("--output", required=True, help="HTML output path")
    parser.add_argument("--title", help="optional figure title")
    parser.add_argument("--z-scale", type=float, default=8.0,
                        help="vertical display scale, default: 8")
    parser.add_argument("--crease-threshold", type=float,
                        help="mark interior edges with gradient jump at least this value")
    args = parser.parse_args()
    make_html(args.input, args.output, args.title, args.z_scale, args.crease_threshold)
    print(args.output)


if __name__ == "__main__":
    main()
