//! Orthographic camera math shared by GPU uniforms and source-indexed picking.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    pub yaw: f64,
    pub pitch: f64,
    pub zoom: f64,
    pub pan: [f64; 2],
}
impl Default for Camera {
    fn default() -> Self {
        Self {
            yaw: 0.0,
            pitch: 0.0,
            zoom: 1.0,
            pan: [0.0; 2],
        }
    }
}
impl Camera {
    pub fn drag(&mut self, dx: f64, dy: f64, pan: bool, width: f64, height: f64) {
        if ![dx, dy, width, height].into_iter().all(f64::is_finite) || width <= 0.0 || height <= 0.0
        {
            return;
        }
        if pan {
            self.pan[0] = (self.pan[0] + 2.0 * dx / width).clamp(-4.0, 4.0);
            self.pan[1] = (self.pan[1] - 2.0 * dy / height).clamp(-4.0, 4.0);
        } else {
            self.yaw = (self.yaw + dx * 0.01).rem_euclid(std::f64::consts::TAU);
            self.pitch = (self.pitch + dy * 0.01)
                .clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        }
    }
    pub fn scroll(&mut self, delta: f64) {
        if delta.is_finite() {
            self.zoom =
                (self.zoom * (-delta.clamp(-1000.0, 1000.0) * 0.001).exp()).clamp(0.1, 20.0);
        }
    }
    pub fn project(&self, point: [f32; 3], width: f64, height: f64) -> [f64; 3] {
        let [x, y, z] = point.map(f64::from);
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        let rx = cy * x + sy * z;
        let rz = -sy * x + cy * z;
        let ry = cp * y - sp * rz;
        let z = sp * y + cp * rz;
        let scale = 0.82 * self.zoom;
        [
            (rx * scale / (width / height) + self.pan[0] + 1.0) * width * 0.5,
            (1.0 - ry * scale - self.pan[1]) * height * 0.5,
            z,
        ]
    }
    pub fn pick(
        &self,
        scene: &collagene_core::structure_view::coordinates::scene::Scene,
        x: f64,
        y: f64,
        width: f64,
        height: f64,
    ) -> Option<u32> {
        if ![x, y, width, height].into_iter().all(f64::is_finite) || width <= 0.0 || height <= 0.0 {
            return None;
        }
        scene
            .points
            .iter()
            .filter_map(|p| {
                let q = self.project(p.position, width, height);
                let d = (q[0] - x).powi(2) + (q[1] - y).powi(2);
                (d <= 36.0).then_some((p.atom_index, q[2], d))
            })
            .max_by(|a, b| a.1.total_cmp(&b.1).then_with(|| b.2.total_cmp(&a.2)))
            .map(|p| p.0)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn picking_returns_frontmost_original_index_and_changes_with_rotation() {
        use collagene_core::structure_view::coordinates::scene::{
            Alternates, Chains, Point, Scene, Selection,
        };
        let scene = Scene {
            source_hash: "a".repeat(64),
            selection: Selection {
                model: 1,
                chains: Chains::All,
                alternates: Alternates::All,
                hydrogens: true,
            },
            points: vec![
                Point {
                    atom_index: 9,
                    position: [0.0, 0.0, -0.5],
                    color: [0.5; 3],
                },
                Point {
                    atom_index: 77,
                    position: [0.0, 0.0, 0.5],
                    color: [0.5; 3],
                },
            ],
            edges: vec![],
            scaled_center: [0.0; 3],
            coordinate_magnitude: 1.0,
            scaled_radius: 1.0,
            unresolved_connections: 0,
            withheld_symmetry_or_relation: 0,
        };
        let c = Camera::default();
        assert_eq!(c.pick(&scene, 400.0, 200.0, 800.0, 400.0), Some(77));
        assert_eq!(c.pick(&scene, 10.0, 10.0, 800.0, 400.0), None);
        assert_eq!(c.pick(&scene, 0.0, 0.0, 0.0, 400.0), None);
        let c = Camera {
            yaw: std::f64::consts::PI,
            ..c
        };
        assert_eq!(c.pick(&scene, 400.0, 200.0, 800.0, 400.0), Some(9));
    }
    #[test]
    fn orthographic_basis_rotation_and_aspect_match_screen_coordinates() {
        let c = Camera::default();
        assert_eq!(c.project([0.0; 3], 800.0, 400.0), [400.0, 200.0, 0.0]);
        let x = c.project([1.0, 0.0, 0.0], 800.0, 400.0);
        assert!((x[0] - 564.0).abs() < 1e-9);
        let y = c.project([0.0, 1.0, 0.0], 800.0, 400.0);
        assert!((y[1] - 36.0).abs() < 1e-9);
        let c = Camera {
            yaw: std::f64::consts::FRAC_PI_2,
            ..c
        };
        let p = c.project([1.0, 0.0, 0.0], 800.0, 400.0);
        assert!((p[0] - 400.0).abs() < 1e-9);
        assert!((p[2] + 1.0).abs() < 1e-9);
    }
    #[test]
    fn pointer_pan_zoom_and_nonfinite_events_remain_bounded() {
        let mut c = Camera::default();
        c.drag(40.0, 20.0, true, 800.0, 400.0);
        let p = c.project([0.0; 3], 800.0, 400.0);
        assert!((p[0] - 440.0).abs() < 1e-9);
        assert!((p[1] - 220.0).abs() < 1e-9);
        c.scroll(-1e308);
        assert!(c.zoom <= 20.0);
        let saved = c;
        c.scroll(f64::NAN);
        c.drag(f64::INFINITY, 1.0, false, 800.0, 400.0);
        assert_eq!(c, saved);
    }
}
